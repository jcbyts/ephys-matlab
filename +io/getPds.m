function PDS = getPds(sessionInfo)
% GET PDS loads PLDAPS files and synchronizes with the OE clock
% Inputs:
%   SessionInfo@struct - session info struct from io.loadSession(oepath)
% Outputs:
%   PDS@cell    - array of PDS structs
% Example call:
%   PDS = io.getPds(sessionInfo)

% 2017.08.14    jly     wrote it

behaviorDir = fullfile(sessionInfo.path, '_behavior');
if exist(behaviorDir, 'dir')
    load(fullfile(behaviorDir, 'PDS.mat'));
    return
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
noEphys = cellfun(@(x) isempty(x.maxreconstructionerror), PDS);

PDS(noEphys | badSync) = [];

dt = cellfun(@(x) x.initialParametersMerged.session.initTime, PDS);
[~, id] = sort(dt);
PDS = PDS(id);

nT = cellfun(@(x) numel(x.data), PDS);
PDS(nT==1) = [];

if ~exist(behaviorDir, 'dir')
    mkdir(behaviorDir)
end

save(fullfile(behaviorDir, 'PDS.mat'), 'PDS')

