function sp = runSingleChannelSpikeSortMog(ops)

if istable(ops)
    thisSession = ops;
    ops = io.loadOps(ops);
else
    thisSession = [];
end

info = load(fullfile(ops.root, 'ephys_info.mat'));
data = io.loadRaw(ops, [], true); % load data in mV

% highpass filter
%[b,a] = butter(3, 300/30e3*2, 'high');

[b,a] = butter(3, [(300/30e3*2), (6000/30e3*2)]);

size(data)

% detect artifacts
%[data, bad] = preprocess.removeChannelArtifacts(data, 1000, 1, 30, false);
%data(bad,:) = 0;

data = data';
data = filter(b,a,data);
data = flipud(data);
data = filter(b,a,data);
data = flipud(data);
% SD_data = std(data);

sp = struct('ss', [], 'clu', []);
nChannels = size(data,2);
Fs = info.sampleRate; % sampling rate
threshhold = -3.5; % in SD
clustOffset = 0;
for iCh = 1:nChannels
    
    fprintf('Spike sorting channel %d\n', iCh)
    % get threshold crossings
    ss = ephys.spikeSorting.detectSpikes(data(:,iCh), Fs, threshhold);
    
    wf = ephys.spikeSorting.extractWaveforms(data(:,iCh), ss);
    
    b = ephys.spikeSorting.extractFeatures(wf);
    try
        df = 3; % degrees of freedom
        model = ephys.spikeSorting.MixtureModel.fit(b, df, 'verbose', true);
        results.s = ss;
        results.w = wf;
        results.b = b;
        results.model = model;
    catch
        fprintf('Mixture model failed. Trying Gaussian model\n')
        df = inf; % degrees of freedom
        model = ephys.spikeSorting.MixtureModel.fit(b, df, 'verbose', true);
        results.s = ss;
        results.w = wf;
        results.b = b;
        results.model = model;
    end
    X0 = ephys.spikeSorting.keepMaxClusters(results, size(data,1), .6);
    
    nUnits = size(X0,2);
    [ss, clu] = find(X0);
    
    figure(1); clf
    cmap = lines;
    for iUnit = 1:nUnits
        subplot(2,nUnits,iUnit)
        spikeSampleTimes = ss(clu==iUnit);
        nSpikes = numel(spikeSampleTimes);
        wf = ephys.spikeSorting.extractWaveforms(data(:,iCh), spikeSampleTimes);
        plot(wf(:,1:ceil(nSpikes/500):end), 'Color', cmap(iUnit,:));
        bins = 0:(Fs/1e3):ceil(.15*Fs);
        cnt = histc(diff(spikeSampleTimes), bins);
        cnt = cnt./sum(cnt);
        subplot(2,nUnits, iUnit+nUnits)
        bar([-bins(2:end), bins]/Fs, [cnt(2:end)' cnt'], 'FaceColor', cmap(iUnit,:));
        xlim([-1 1]*max(bins)/Fs);
        xlabel('sec')
    end
    drawnow
    
    sp.ss = [sp.ss; ss(:)];
    sp.clu = [sp.clu; clu(:)+clustOffset];
    
    clustOffset = max(sp.clu);
end

%% save spikes
sp.st = io.convertSamplesToTime(sp.ss, Fs, info.timestamps(:), info.fragments(:));
sp.cids = unique(sp.clu)';

load(ops.chanMap)
sp.yc = ycoords;
sp.ycoords = ycoords;
sp.xc = xcoords;
sp.xcoords = xcoords;
n = numel(sp.cids);
sp.cgs = zeros(n,1);
sp.cgs2= zeros(n,1);
sp.cR  = zeros(n,1);
sp.clusterAmps   = zeros(1,n);
sp.clusterDepths = zeros(1,n);
sp.firingRates   = zeros(1,n);
sp.isiV = zeros(1,n);


save(fullfile(ops.root, 'sp-MOG.mat'), '-v7.3', '-struct', 'sp')


if istable(thisSession)
    newThisSession = thisSession;
    
    ssList = thisSession.SpikeSorting{1};
    if ~isnan(ssList)
        ssList = regexp(ssList, ',', 'split');
        ssList = union(ssList, {'MOG'});
        ssList = sprintf('%s,', ssList{:});
        newThisSession.SpikeSorting{1} = ssList(1:end-1);
    else
        ssList = {'MOG'};
        ssList = sprintf('%s,', ssList{:});
        newThisSession.SpikeSorting{1} = ssList(1:end-1);
    end
    
    io.writeMeta(newThisSession, 2)
    
end