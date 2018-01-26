%% get all sessions with hartley FF
SERVER_DIR = getpref('EPHYS', 'SERVER_DATA');
LOCAL_DIR  = getpref('EPHYS', 'LOCAL_DATA');

datadir = fullfile(LOCAL_DIR, 'hartleyFFdata');
mkdir(datadir)

meta = io.getExperimentsAnd('Subject', 'Ellie', 'Chamber', 'V1', 'Lens', 1, 'StimulusProtocols', 'hartleyFF', 'SpikeSorting', 'Kilo');

nSessions = size(meta,1);

for iSession = 1:nSessions
    
    thisSession = meta(iSession,:);
    
    
    
    [sess, ops, info] = io.loadSession(thisSession);
    
    PDS = io.getPds(sess);
    
    hart = session.hartleyFF(PDS);
    
    frameTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), hart.trial, 'UniformOutput', false));
    kx         = cell2mat(arrayfun(@(x) x.kx(:), hart.trial, 'UniformOutput', false));
    ky         = cell2mat(arrayfun(@(x) x.ky(:), hart.trial, 'UniformOutput', false));
    
    % get eye data
    [eyeDeg, timestamps, elInfo] = io.getEdf(sess);
    
    % detect saccades
    ix = any(eyeDeg(1:2,:)==elInfo.bitDeg(2));
    eyeDeg(1,ix) = nan;
    eyeDeg(2,ix) = nan;
    
    % --- detect those saccades
    ix = ~any(isnan(eyeDeg));
    [saccades] = pdsa.detectSaccades(timestamps(ix), eyeDeg(1:2,ix), 'verbose', false);
    
    sp = io.getSpikes(sess);
    
    sp = sp{1};
    
    sp.cids = unique(sp.clu);
    
    thresh = 10;
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
    
%     figure(1); clf
%     plot.spikeWaveformsFromOps(ops, sp_new, 'numWaveforms', 15)
    
    sp = sp_new;
    fname = fullfile(datadir, [thisSession.Directory{1} '.mat']);
    save(fname, '-v7', 'frameTimes', 'saccades', 'kx', 'ky', 'sp');
end

%%
tmp = dir(fullfile(datadir, '*.mat'));

flist = arrayfun(@(x) fullfile(x.folder, x.name), tmp, 'UniformOutput', false);
tar(fullfile(datadir, 'data.tgz'), flist)

b = untar(fullfile(datadir, 'data.tgz'))