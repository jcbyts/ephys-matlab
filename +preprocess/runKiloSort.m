function thisSession = runKiloSort(thisSession, varargin)
% runKiloSort runs Kilosort on data files
% It can either take in an ops struct (as specified in the Kilosort
% github) or it can take session meta table
%
% e.g.,
% ops = io.loadOps(fpath)
% preprocess.runKiloSort(ops(1))
%
% meta = io.getExperimentsAnd();
% preprocess.runKilosort(meta(1,:));
%
% Optional Arguments:
%   'merge'    [logical] - run post sorting routine to automatically merge
%                          units based on feature similarity and ISI 
%                          distribution

% parse relevant optional arguments
ip = inputParser();
ip.KeepUnmatched = true; % important. arguments that are not matched can be passed on
ip.addOptional('overwrite', false)
ip.parse(varargin{:})
        
% detect if this is a meta object
ops = thisSession;

if isa(ops, 'table') % Ops is a meta data session
    disp('Session Meta detected.')
    
    if isnan(ops.oe2dat) || ~ops.oe2dat
        disp('oe2dat has not been run. aborting')
        return
    end
    
    % Check if KiloSort has already been run
    if ~any(strcmp(ops.SpikeSorting, 'Kilo'))
        disp('Kilosort has never been run')
        disp('Running Kilosort now')
        
        % Here the meta data is converted into a Kilosort Ops struct
        ops = io.loadOps(ops);
    else
        disp('Kilosort has already been run')
        
        if ip.Results.overwrite
            disp('overwrite is true')
            ops = io.loadOps(ops);
            disp('Running Kilosort now')
        else % exit
            disp('skipping')
            return
        end
    end
        
end

assert(isstruct(ops), 'ops is not a struct for some reason')
nOps = numel(ops);
fprintf('Found %d ops structs\n', nOps)

% Loop over ops structs and run KiloSort
success = false(1,nOps);
kpaths  = cell(nOps,1);
for i = 1:nOps
    if ops(i).Nchan==1
        continue
    else
        [success(i), kpaths{i}] = run_KiloSort_helper(ops(i), varargin{:});
    end
end

if istable(thisSession)
    newThisSession = thisSession;
    
    ssList = thisSession.SpikeSorting{1};
    if ~isnan(ssList)
        ssList = regexp(ssList, ',', 'split');
        ssList = union(ssList, {'Kilo'});
        ssList = sprintf('%s,', ssList{:});
        newThisSession.SpikeSorting{1} = ssList(1:end-1);
    else
        ssList = {'Kilo'};
        ssList = sprintf('%s,', ssList{:});
        newThisSession.SpikeSorting{1} = ssList(1:end-1);
    end
    
    io.writeMeta(newThisSession, 2)
    thisSession = newThisSession;
end

if any(success)
    commandwindow
    fprintf('KiloSort has been Successfully run\n')
    
    fprintf('Open the START menu and select Anaconda prompt\n')

    fprintf('Copy and paste in the following\n')

    fprintf('activate phy\n')
    for i = 1:nOps
        fprintf('cd %s\n', kpaths{i})
        fprintf('phy template-gui params.py\n\n\n')
    end
        
    fprintf('Make sure to copy cluster_groups.csv to the server directory')
else
    fprintf('Uh oh. Kilosort failed for some reason.\n')
end



% -------------------------------------------------------------------------
% --- Helper function: This actually calls KiloSort
function [success, phypath] = run_KiloSort_helper(ops, varargin)

ip = inputParser();
ip.addOptional('merge',     false)
ip.addOptional('GPU',       true)
ip.addOptional('parfor',    true)
ip.addOptional('verbose',   true)
ip.addOptional('showfigures', false)
ip.addOptional('overwrite', false)

ip.parse(varargin{:});

% check if Kilosort has already been run
if exist(fullfile(ops.root,  'rez.mat'), 'file')
    fprintf('Kilosort has already been run. Skipping\n')
    success = true;
    phypath = [];
    return
end

fprintf('-----------------------------------------------------------------\n')
fprintf('-----------------------------------------------------------------\n')
fprintf('-----------------------------------------------------------------\n')
fprintf('Running Kilosort\n')
OLD_DIR = ops.root;

% --- Copy the data to the local disk before running.

% 1) create a local directory
LOCAL_DIR = fullfile(getpref('EPHYS', 'LOCAL_DATA'), ['kilo-tmp-' datestr(now, 'yyyymmddHHMM')]);
mkdir(LOCAL_DIR)

% 2) convert the ops struct to point to the local directory
ops = io.convertOpsToNewDirectory(ops, LOCAL_DIR);
phypath = ops.root; % save for easy output

% 3) create subdirectory for the individual electrode
[~,~,mess] = mkdir(ops.root);
assert(~strcmp(mess, 'MATLAB:MKDIR:DirectoryExists'), 'tmp directory exists. clean it up.')

% 4) copy the necessary files over
copyfile(fullfile(OLD_DIR, 'ephys.dat'), ops.fbinary)
copyfile(fullfile(OLD_DIR, 'ephys_info.mat'), fullfile(ops.root, 'ephys_info.mat'))
copyfile(fullfile(OLD_DIR, 'chanMap.mat'), ops.chanMap)

% 5) correct the parfor_ naming issue
if isfield(ops, 'parfor_')
    ops.parfor = ops.parfor_;
    ops = rmfield(ops, 'parfor_');
end

% 6) set optional arguments for running kilosort
merge           = ip.Results.merge;
ops.GPU         = ip.Results.GPU;
ops.parfor      = ip.Results.parfor;
ops.verbose     = ip.Results.verbose;
ops.showfigures = ip.Results.showfigures;

% 7) prepare to run
if ops.GPU
    fprintf('Using GPU for faster processing\n')
    fprintf('Opening gpuDevice. If this is the first time this was run during \nthis matlab session, this can be very slow\n')
    gpu = gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
	ops.ForceMaxRAMforDat = min(gpu.AvailableMemory, ops.ForceMaxRAMforDat);
end

% 8) main spike-sorting routine
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization

pause(.5) % pause between functions to keep KiloSort from crashing

rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively

pause(.5) % pause between functions to keep KiloSort from crashing

rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% 9) Cleanup GPU
if ops.GPU
    reset(gpu)
end

% 10) Set Ops to point back to the server directory
ops = io.convertOpsToNewDirectory(ops, OLD_DIR);

% 11) Save to output for phy (Save to server)
fprintf('saving matlab results file\n')
save(fullfile(ops.root,  'rez.mat'), 'rez', 'ops', '-v7.3');

fprintf('saving python files for Phy\n')
rezToPhy(rez, ops.root);

% 12) Save the high-pass filtered data for visualization in phy
fprintf('Saving high-pass filtered data for phy\n')
fprintf('Copying whitened datafile to server\n')

% % build filter
% if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
%     [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
% else
%     [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
% end
% 
% % load Raw
% dataRAW = io.loadRaw(ops, [], false);
% dataRAW = double(dataRAW');
% 
% % common average reference
% dataRAW = bsxfun(@minus, dataRAW, mean(dataRAW,2));
%     
% % filter
% datr = filter(b1, a1, dataRAW);
% datr = flipud(datr);
% datr = filter(b1, a1, datr);
% datr = flipud(datr);
% 
% % convert to a 16-bit integer
% datr = int16(datr');
% 
% % open file for writing
% fid = fopen(ops.fproc, 'w');
% 
% % write to disk
% fwrite(fid, datr, 'int16');
%     
% % cleaup
% fclose(fid);
% fprintf('Done\n')
preprocess.save_filtered_data(ops);

% 13) Run the automatic merge
if merge
    fprintf('Attemptying Automerge\n')
    rez                = merge_posthoc3(rez);
    fprintf('saving python files for Phy\n')
    % save python results file for Phy
    rezToPhy(rez, ops.root);
end

fprintf('Kilo tmp folder is [%s]\n', LOCAL_DIR)

% 14) Copy files back over to the local disk for running phy
% Yes, we probably could've done this in a different order (feels like
% we're copying back and forth), but this was slapped on at the end and it
% works. It can be revisited later.

% copy files that are required for phy back to the local disk
phy_files = {'params.py', 'amplitudes.npy', 'channel_map.npy', 'channel_positions.npy', 'pc_features.npy', ...
    'pc_feature_ind.npy', 'similar_templates.npy', 'spike_clusters.npy', 'spike_templates.npy', ...
    'spike_times.npy', 'templates.npy', 'templates_ind.npy', 'template_features.npy', ...
    'template_feature_ind.npy', 'whitening_mat.npy', 'whitening_mat_inv.npy', 'rez.mat', 'tempWh.dat'};

newPath = regexp(ops.root, '\\', 'split');

fprintf('Copying files for phy over to the local disk\n')
for iFile = 1:numel(phy_files)
    server_file = fullfile(ops.root, phy_files{iFile});
    if exist(server_file, 'file')
        fprintf('%[s]\n', phy_files{iFile})
        copyfile(server_file, fullfile(LOCAL_DIR, newPath{end}, phy_files{iFile}))
    end
end

fprintf('Done')

success = true;
