function PDS = getPds(sessionInfo, overwrite, forceAll)
% GET PDS loads PLDAPS files and synchronizes with the OE clock
% Inputs:
%   SessionInfo@struct - session info struct from io.loadSession(oepath)
%   overwrite@logical  - flag
%   includeBadSyncs@logical - keep sessions that aren't matched with ephys
% Outputs:
%   PDS@cell    - array of PDS structs
% Example call:
%   PDS = io.getPds(sessionInfo, overwrite, includeBadSyncs)

% 2017.08.14    jly     wrote it

if nargin < 3
    forceAll = false;
    if nargin < 2
        overwrite = false;
    end
end

if isa(sessionInfo, 'table') % it's meta data -- load the struct
    sessionInfo = io.loadSession(sessionInfo);
end

behaviorDir = fullfile(sessionInfo.path, '_behavior');
try
    if exist(behaviorDir, 'dir') && ~overwrite
        tmp = load(fullfile(behaviorDir, 'PDS.mat'));
        PDS = tmp.PDS;
        return
    end
catch
    warning('Failed to load PDS. Going to try re-importing')
end

pdsList = io.getPdsFileList(sessionInfo);

evList = io.getEventsFileList(sessionInfo);

nPdsFiles = numel(pdsList);

PDS = cell(nPdsFiles,1);
fprintf('Found [%d] PDS files\n', nPdsFiles)


reconError = nan(nPdsFiles,1);

for kPdsFile = 1:nPdsFiles
    
    tmp = load(pdsList{kPdsFile}, '-mat');
    
    % check that you have the proper PLDAPS and PEP branches currently in
    % the path to ensure that all objects are recreated properly
    
    if tmp.PDS.initialParametersMerged.git.use
        branch = regexp(tmp.PDS.initialParametersMerged.git.pep.status, '(?<=branch\s)(\w+)', 'match');
        data_branch = branch{1};
        [pep_path, ~] = fileparts(which('calibrationGUI.m'));
        curr_path = pwd;
        
        cd(pep_path);
        stat = pds.git.git('status');
        cd(curr_path)
        branch = regexp(stat, '(?<=branch\s)(\w+)', 'match');
        curr_branch = branch{1};
        assert(strcmp(data_branch, curr_branch), 'You are on the wrong branch of PEP')
    else
        reposDir = getpref('ephysmatlab', 'repos');
        old_pep_path = fullfile(reposDir, 'pds-stimuli-pldapsGUI\');
        addpath(old_pep_path)
        % reload PDS file
        tmp = load(pdsList{kPdsFile}, '-mat');
    end
    
   
    
    if ~isempty(tmp.PDS.functionHandles)
        fhlist = fieldnames(tmp.PDS.functionHandles);
        for i = 1:numel(fhlist)
            if isa(tmp.PDS.functionHandles.(fhlist{i}), 'matlab.ui.Figure')
                close(tmp.PDS.functionHandles.(fhlist{i}))
            end
        end
    end
    
    if ~forceAll
        % --- Synchronize the clocks ------------------------------------------
        [OE2PTBfit, ~,PTB2OE, maxreconstructionerror] = io.sync2OeClock(tmp.PDS, evList);
        
        tmp.PDS.PTB2OE    = PTB2OE;
        tmp.PDS.OE2PTBfit = OE2PTBfit;
        tmp.PDS.maxreconstructionerror = maxreconstructionerror;
        
        if isempty(maxreconstructionerror)
            fprintf('%d) No Ephys Data\n', kPdsFile);
        else
            fprintf('%d) Strobe times aligned. Max reconstruction error is %2.3f ms\n', kPdsFile, maxreconstructionerror*1e3)
            reconError(kPdsFile) = maxreconstructionerror*1e3;
        end
        
    end
    
    PDS{kPdsFile} = tmp.PDS;
end

fprintf('Aligning with OE clock now\n')

if forceAll
    
    trial = pds.getPdsTrialData(PDS);
    % --- Synchronize the clocks ------------------------------------------
    fprintf('Forcing all PDS files to be synced together\n')
    
    if isempty(evList)
        fprintf('no ephys this session\n')
        OE2PTBfit = @(x) x;
        PTB2OE = @(x) x;
        maxreconstructionerror = inf;
    else
        [OE2PTBfit, ~,PTB2OE, maxreconstructionerror] = io.sync2OeClock_trial(trial, evList);
        fprintf('Strobe times aligned. Max reconstruction error is %2.3f ms\n', maxreconstructionerror*1e3)
    end
    
    for kPdsFile = 1:numel(PDS)
        PDS{kPdsFile}.OE2PTBfit = OE2PTBfit;
        PDS{kPdsFile}.PTB2OE = PTB2OE;
        PDS{kPdsFile}.maxreconstructionerror = maxreconstructionerror;
    end
else
    badSync = reconError > .1;
    noEphys = cellfun(@(x) isempty(x.maxreconstructionerror), PDS);
    PDS(noEphys | badSync) = [];
end

dt = cellfun(@(x) x.initialParametersMerged.session.initTime, PDS);
[~, id] = sort(dt);
PDS = PDS(id);

nT = cellfun(@(x) numel(x.data), PDS);
PDS(nT==1) = [];

if ~exist(behaviorDir, 'dir')
    mkdir(behaviorDir)
end


% remove figures embedded within the PDS files (they are big!)
for iFile = 1:numel(PDS)
    if isempty(PDS{iFile}.functionHandles)
        continue
    else
        fn_ = fieldnames(PDS{iFile}.functionHandles);
        for k = 1:numel(fn_)
%             disp(class(fn_{k}))
            if isa(PDS{iFile}.functionHandles.(fn_{k}), 'matlab.ui.Figure')
                PDS{iFile}.functionHandles.(fn_{k}) = [];
                disp('Removing figure from file')
            end
        end
    end
end

fprintf('Saving PDS files out as a mat file\n')
save(fullfile(behaviorDir, 'PDS.mat'), 'PDS', '-v7.3')

