%% get all sessions with hartley FF
SERVER_DIR = getpref('EPHYS', 'SERVER_DATA');
LOCAL_DIR  = getpref('EPHYS', 'LOCAL_DATA');

datadir = fullfile(LOCAL_DIR, 'hartleyFFdata');
mkdir(datadir)

meta = io.getExperimentsAnd('Subject', 'Ellie', 'Chamber', 'V1', 'Lens', 1, 'StimulusProtocols', 'hartleyFF', 'SpikeSorting', 'Kilo');

nSessions = size(meta,1);
%%
for iSession = 10:nSessions
    
    thisSession = meta(iSession,:);
    
    
    
    [sess, ops, info] = io.loadSession(thisSession);
    warning off
    PDS = io.getPds(sess);
    warning on 
    
    hart = session.hartleyFF(PDS);
    
    frameTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), hart.trial, 'UniformOutput', false));
    kx         = cell2mat(arrayfun(@(x) x.kx(:), hart.trial, 'UniformOutput', false));
    ky         = cell2mat(arrayfun(@(x) x.ky(:), hart.trial, 'UniformOutput', false));
    
    % frozenSequence starts, duration
    frozenTrials = arrayfun(@(x) x.frozenSequence, hart.trial);
    if ~any(frozenTrials)
        frozen_seq_starts = [];
        frozen_seq_dur = 0;
    else
        frozen_seq_dur = arrayfun(@(x) x.frozenSequenceLength, hart.trial(frozenTrials));
        assert(all(frozen_seq_dur == mode(frozen_seq_dur)), 'frozen sequences are of different length')
        frozen_seq_dur = mode(frozen_seq_dur);
        trialDuration = arrayfun(@(x) numel(x.frameTimes), hart.trial);
        trialStartFrames = 1 + [0; cumsum(trialDuration((1:end-1)))];
        frozen_seq_starts = [];
        for iTrial = find(frozenTrials(:))'
            numRepeats = floor(trialDuration(iTrial)/frozen_seq_dur);
            trial_seq_starts = trialStartFrames(iTrial) + (0:(numRepeats-1))*frozen_seq_dur;
            frozen_seq_starts = [frozen_seq_starts trial_seq_starts]; %#ok<AGROW>
        end
    end
    
    
    % get eye data
    [eyeDeg, timestamps, elInfo] = io.getEdf(sess);
    
    % detect saccades
    ix = any(eyeDeg(1:2,:)==elInfo.bitDeg(2));
    eyeDeg(1,ix) = nan;
    eyeDeg(2,ix) = nan;
    
    % --- detect those saccades
    ix = ~any(isnan(eyeDeg));
    [saccades] = pdsa.detectSaccades(timestamps(ix), eyeDeg(1:2,ix), 'verbose', false);
    
    eyeData = [timestamps(:) eyeDeg'];
    
    
    sp = io.getSpikes(sess);
    
    iShank = 1; % TODO: how to handle multiple shanks
    sp = sp{iShank};
    
    if numel(sp.cids) ~= numel(unique(sp.clu))
        sp = io.getSpikesFromKilo(ops(1));
    end
   
    
    thresh = 0;
    ix = find(sp.uQ > thresh);
    [cdepth, id] = sort(sp.clusterDepths(ix));
    ix = ix(id);
    
    nUnits = numel(ix);
    
    sp_new = struct('st', [], 'clu', [], 'ss', []);
    for iUnit = 1:nUnits
        ii = sp.clu == sp.cids(ix(iUnit));
        stu = sp.st(ii);
        sp_new.st   = [sp_new.st; stu];
        sp_new.ss   = [sp_new.ss; sp.ss(ii)];
        sp_new.clu  = [sp_new.clu; iUnit*ones(numel(stu),1)];
    end
    
    sp_new.cids = unique(sp_new.clu);
    sp_new.yc = sp.yc;
    sp_new.xc = sp.xc;
    sp_new.clusterDepths = cdepth;
    sp_new.uQ = sp.uQ(ix);
	sp_new.cgs = sp.cgs(ix);
    
    % --- calculate cluster depths with CSD 
    csdF = session.csdFlash(PDS);

    figure; clf
    csdstats = csdF.computeCsd(ops(iShank), 'plot', true);
    title(thisSession.Directory{1})
    drawnow
    
    sp = sp_new;
    fname = fullfile(datadir, [thisSession.Directory{1} '.mat']);
    save(fname, '-v7', 'frameTimes', 'saccades', 'kx', 'ky', 'sp', 'frozen_seq_dur', 'frozen_seq_starts', 'eyeData', 'csdstats');
end

%%
tmp = dir(fullfile(datadir, '*.mat'));

flist = arrayfun(@(x) fullfile(x.folder, x.name), tmp, 'UniformOutput', false);
tar(fullfile(datadir, 'data.tgz'), flist)

b = untar(fullfile(datadir, 'data.tgz'))