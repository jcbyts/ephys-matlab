function PDS = getPds(sessionInfo, overwrite, includeBadSyncs)
% GET PDS loads PLDAPS files and synchronizes with the OE clock
% Inputs:
%   SessionInfo@struct - session info struct from io.loadSession(oepath)
% Outputs:
%   PDS@cell    - array of PDS structs
% Example call:
%   PDS = io.getPds(sessionInfo)

% 2017.08.14    jly     wrote it

if nargin < 3
    includeBadSyncs = false;
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
fprintf('Aligning with OE clock now\n')

reconError = nan(nPdsFiles,1);

for kPdsFile = 1:nPdsFiles
    
    tmp = load(pdsList{kPdsFile}, '-mat');
    
    % check that you have the proper PLDAPS and PEP branches currently in
    % the path to ensure that all objects are recreated properly
    
    % I hate regexp! why is it so confusing!?
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
    
    
    
    if ~isempty(tmp.PDS.functionHandles)
        fhlist = fieldnames(tmp.PDS.functionHandles);
        for i = 1:numel(fhlist)
            if isa(tmp.PDS.functionHandles.(fhlist{i}), 'matlab.ui.Figure')
                close(tmp.PDS.functionHandles.(fhlist{i}))
            end
        end
    end
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
    
    
    PDS{kPdsFile} = tmp.PDS;
end


badSync = reconError > .1;
if any(badSync)
    keyboard
end
noEphys = cellfun(@(x) isempty(x.maxreconstructionerror), PDS);
if includeBadSyncs
    PDS(noEphys) = [];
else
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

