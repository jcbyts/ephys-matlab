function importSession(thisSession)
% this is a bottom up import of a session. Some things may already have
% been run. 
% call by passing in an existing session meta data or by starting from
% scratch
% example calls:
% meta = io.getMetaTable();
% 
% io.importSession(meta(10,:)); % run import on session 10 from the meta
%                               table
% 
% or 
% io.importSession(); % no arguments lets you select things


SERVER_DATA_DIR = getpref('EPHYS', 'SERVER_DATA');

if nargin == 0 % no session passed in
    % use gui to select session
    directoryname = uigetdir(SERVER_DATA_DIR, 'Pick session to import');
    meta = io.getMetaTable();
    directory = strrep(directoryname, SERVER_DATA_DIR, ''); % path relative to server data dir
    sessionix = find(strcmp(meta.Directory, directory));
    if isempty(sessionix) % session does not exist yet
        pat = '(?<subject>\D+)\_(?<date>[\d-]+)\_(?<time>[\d-]+)\_(?<note>[^]+)';
        s = regexp(directory, pat, 'names'); % build sesion from filename
        assert(~isempty(s.subject), 'This is not a session') % TODO: this should be a more robust check
        
        % get the current variable names in the meta table
        tableColumns = meta.Properties.VariableNames;
        tableValues  = repmat({nan}, 1, numel(tableColumns)); % populate empty values
        
        % --- populate specifics we can know from the file name
        
        % subject name
        ix = strcmp(tableColumns, 'Subject');
        tableValues{ix} = s.subject;
        
        % date
        ix = strcmp(tableColumns, 'Date');
     
        dstr = cellfun(@str2double, regexp(s.date, '-', 'split'), 'uni', false);
        tableValues{ix} = datestr(datenum(dstr{:}), 'mm/dd/yyyy');
        
        % directory
        ix = strcmp(tableColumns, 'Directory');
        tableValues{ix} = directory;
        
        % time
        ix = strcmp(tableColumns, 'Time');
        tableValues{ix} = strrep(s.time, '-', ':');
        
        % tag
        ix = strcmp(tableColumns, 'Tag');
        tableValues{ix} = s.note;
        
        
        thisSession = cell2table(tableValues, 'VariableNames', tableColumns);
        
        
    else % session has already been added to the meta file. Load what is already there
        thisSession = meta(sessionix,:);
    end
        
keyboard
end


disp(thisSession)

electrodeName = thisSession.Electrode;
if iscell(electrodeName)
    electrodeName = electrodeName{1};
end

if isnan(electrodeName)
    disp('No electrode data. Select from list of available.')
    shank = hardware.electrodeFactory();
else
    shank = hardware.electrodeFactory(electrodeName);
end


oepath = fullfile(SERVER_DATA_DIR, thisSession.Directory{1});

% check if import has already been run
derived = dir(fullfile(oepath, '_*'));
try
    for iDir = 1:numel(derived)
        ops_old = io.loadOps(oepath);
        
        for i = 1:numel(ops_old)
            o = ops_old(i);
            
            o = io.convertOpsToNewDirectory(o, oepath);
            
            save(fullfile(o.root, 'ops.mat'), '-v7.3', '-struct', 'o')
        end
    end
end

io.oe2dat(oepath, shank, 'overwrite', false, 'verbose', true);

thisSession.oe2dat = true;
% --- load session: This is central to all operations
sess = io.loadSession(oepath);



overwrite = false; % if it breaks, run again with true
PDS = io.getPds(sess, overwrite);


io.getEdf(sess, PDS, overwrite);


%% Step 4: Import LFP
% Loop over sessions and import the local field potential
% This is really slow because of the line-noise fitting and correction
ops = io.loadOps(oepath);

plotIt = false;
for i = 1:numel(ops)
    [~, ~, info] = io.getLFP(ops(i),overwrite,plotIt);
end

if info.phaseCorrection
    thisSession.LfpPhaseCorrection = true;
end

io.writeMeta(thisSession)
