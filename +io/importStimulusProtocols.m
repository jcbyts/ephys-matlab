function newThisSession = importStimulusProtocols(thisSession)

newThisSession = thisSession;

% --- import the stimulus
disp('Loading PDS files')
PDS = io.getPds(thisSession);

disp('Checking for available stimulus protocols')
ephys_root = fileparts(which('EphysSession'));
stimulusProtocols = arrayfun(@(x) x.name(1:end-2), dir(fullfile(ephys_root, '+session', '*.m')), 'uni', false);

nProt = numel(stimulusProtocols);
hasStim = false(nProt,1);
for kStim = 1:nProt
    tmp = feval(str2func(sprintf('session.%s', stimulusProtocols{kStim})), PDS);
    tmp
    if tmp.numTrials > 1
        hasStim(kStim) = true;
    end
end

stimList = thisSession.StimulusProtocols{1};
if ~isnan(stimList)
    stimList = regexp(stimList, ',', 'split');
    stimList = union(stimList, stimulusProtocols(hasStim));
    stimList = sprintf('%s,', stimList{:});
    newThisSession.StimulusProtocols{1} = stimList(1:end-1);
else
    stimList = stimulusProtocols(hasStim);
    stimList = sprintf('%s,', stimList{:});
    newThisSession.StimulusProtocols{1} = stimList(1:end-1);
end