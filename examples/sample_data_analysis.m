
meta = io.getExperimentsAnd('Subject', 'Ellie');

figure(1); clf
plot(datenum(meta.Date, 'mm/dd/yyyy'), meta.Weight, '-o')
datetick(gca)
ylabel('Weight (g)')


%% get all V1 sessions

meta = io.getExperimentsAnd('Subject', 'Ellie', 'Chamber', 'V1', 'Lens', 1);
new_meta = meta;
nSessions = size(meta,1);
fprintf('Found %d sessions that match the criteria\n', nSessions)

% loop through and 
for iSession = 1:nSessions
    disp(iSession)
    oepath = fullfile(getpref('EPHYS', 'SERVER_DATA'), meta.Directory{iSession});
    
    try
        [sess, ops, info] = io.loadSession(oepath);
    catch me
        disp('Oops. Looks like this session hasn''t been imported yet.')
        io.importSession(meta(iSession,:))
    end
        
    % --- important check which stimulus protocols have been run
    PDS = io.getPds(sess);
    if isempty(PDS)
        continue
    end
    
    stimulusProtocols = {'hartleyFF', 'csdFlash', 'natImgBackground', 'squareFlash', 'psaForage'};
    nProt = numel(stimulusProtocols);
    hasStim = false(nProt,1);
    for kStim = 1:nProt
        tmp = feval(str2func(sprintf('session.%s', stimulusProtocols{kStim})), PDS);
        if tmp.numTrials > 1
            hasStim(kStim) = true;
        end
    end
    
    thisSession = meta(iSession,:);
    stimList = thisSession.StimulusProtocols{1};
    if ~isnan(stimList)
        stimList = regexp(stimList, ',', 'split');
        stimList = union(stimList, stimulusProtocols(hasStim));
        stimList = sprintf('%s,', stimList{:});
        new_meta.StimulusProtocols{iSession} = stimList(1:end-1);
    else
        stimList = stimulusProtocols(hasStim);
        stimList = sprintf('%s,', stimList{:});
        new_meta.StimulusProtocols{iSession} = stimList(1:end-1);
    end
    
    
    spList = thisSession.SpikeSorting{1};
    if all(~isnan(spList)) && any(strfind(spList, 'Kilo'))
        continue
    end
        
    
    for i = 1:numel(ops)
        preprocess.runKiloSort(ops(i), 'GPU', false)
        preprocess.KiloAutomerge(ops(i))
        io.getSpikesFromKilo(ops(i))
    end
    
    
    
    if ~isnan(spList)
        spList = regexp(spList, ',', 'split');
        spList = union(spList, {'Kilo'});
        spList = sprintf('%s,', spList{:});
        new_meta.SpikeSorting{iSession} = spList(1:end-1);
    else
        stimList = {'Kilo'};
        stimList = sprintf('%s,', stimList{:});
        new_meta.SpikeSorting{iSession} = stimList(1:end-1);
    end
%     sp = io.getSpikes(sess);
    
    
end

%%
io.writeMeta(new_meta)

%% get all sessions with hartley FF 
meta = io.getExperimentsAnd('Subject', 'Ellie', 'Chamber', 'V1', 'Lens', 1, 'StimulusProtocols', 'hartleyFF', 'SpikeSorting', 'Kilo');
