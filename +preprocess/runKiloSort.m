function success = runKiloSort(ops, varargin)
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

% detect if this is a meta object
if isa(ops, 'table')
    disp('Session Meta detected.')
    if isnan(ops.oe2dat) || ~ops.oe2dat
        disp('oe2dat has not been run. aborting')
        return
    end
    
    if ~any(strcmp(ops.SpikeSorting, 'Kilo'))
        disp('Kilosort has never been run')
        disp('Running Kilosort now')
        
        ops = io.loadOps(ops);
    else
        disp('Kilosort has already been run')
        ip = inputParser();
        ip.KeepUnmatched = true;
        ip.addOptional('overwrite', false)
        ip.parse(varargin{:})
        if ip.Results.overwrite
            disp('overwrite is true')
            ops = io.loadOps(ops);
            disp('Running Kilosort now')
        else
            success = true;
            disp('skipping')
            return
        end
    end
        
end

assert(isstruct(ops), 'ops is not a struct for some reason')
fprintf('found %d ops structs\n', numel(ops))
success = false(1,numel(ops));
for i = 1:numel(ops)
    if ops(i).Nchan==1
        continue
    else
        success(i) = run_KiloSort_helper(ops(i), varargin{:});
    end
end


function success = run_KiloSort_helper(ops, varargin)

ip = inputParser();
ip.addOptional('merge',     true)
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
    return
end

fprintf('-----------------------------------------------------------------\n')
fprintf('-----------------------------------------------------------------\n')
fprintf('-----------------------------------------------------------------\n')
fprintf('Running Kilosort\n')
OLD_DIR = ops.root;
LOCAL_DIR = fullfile(getpref('EPHYS', 'LOCAL_DATA'), ['kilo-tmp-' datestr(now, 'yyyymmddHHMM')]);
mkdir(LOCAL_DIR)
% ops_old = ops;
ops = io.convertOpsToNewDirectory(ops, LOCAL_DIR);

[~,~,mess] = mkdir(ops.root);
assert(~strcmp(mess, 'MATLAB:MKDIR:DirectoryExists'), 'tmp directory exists. clean it up.')
    
copyfile(fullfile(OLD_DIR, 'ephys.dat'), ops.fbinary)
copyfile(fullfile(OLD_DIR, 'chanMap.mat'), ops.chanMap)

if isfield(ops, 'parfor_')
    ops.parfor = ops.parfor_;
    ops = rmfield(ops, 'parfor_');
end

merge           = ip.Results.merge;
ops.GPU         = ip.Results.GPU;
ops.parfor      = ip.Results.parfor;
ops.verbose     = ip.Results.verbose;
ops.showfigures = ip.Results.showfigures;

if ops.GPU
    fprintf('Using GPU for faster processing\n')
    fprintf('Opening gpuDevice. If this is the first time this was run during \nthis matlab session, this can be very slow\n')
    gpu = gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
    %     gpu = parallel.gpu.GPUDevice.getDevice(1);
	ops.ForceMaxRAMforDat = min(gpu.AvailableMemory, ops.ForceMaxRAMforDat);
end

% main spike-sorting routine
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
pause(.5) % sometimes kilosort crashes. try pausing
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
pause(.5) % sometimes kilosort crashes. try pausing
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

if ops.GPU
    reset(gpu)
end

ops = io.convertOpsToNewDirectory(ops, OLD_DIR);

% saving
fprintf('saving matlab results file\n')
save(fullfile(ops.root,  'rez.mat'), 'rez', 'ops', '-v7.3');

if merge
    fprintf('Attemptying Automerge\n')
    rez                = merge_posthoc3(rez);
end

fprintf('saving python files for Phy\n')

% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
% delete(ops.fproc);
rmdir(LOCAL_DIR, 's')

success = true;

