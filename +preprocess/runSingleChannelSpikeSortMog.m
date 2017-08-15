function sp = runSingleChannelSpikeSortMog(ops)

info = load(fullfile(ops.root, 'ephys_info.mat'));
data = io.loadRaw(ops, [], true); % load data in mV

% highpass filter
[b,a] = butter(3, 300/30e3*2, 'high');

data = data';
data = filter(b,a,data);
data = flipud(data);
data = filter(b,a,data);
data = flipud(data);

% detect artifacts
[data, bad] = preprocess.removeChannelArtifacts(data, 1000, 1, 30, false);
data(bad,:) = 0;

figure(1); clf
plot(data)
ylabel('mV')
xlabel('samples')

sp = struct('ss', [], 'clu', []);

nChannels = size(data,2);
Fs = info.sampleRate; % sampling rate
threshhold = -4; % in SD
clustOffset = 0;
for iCh = 1:nChannels
    
    fprintf('Spike sorting channel %d\n', iCh)
    % get threshold crossings
    ss = ephys.spikeSorting.detectSpikes(data(:,iCh), Fs, threshhold);
    
    wf = ephys.spikeSorting.extractWaveforms(data(:,iCh), ss);
    
    b = ephys.spikeSorting.extractFeatures(wf);
    
    df = 3; % degrees of freedom
    model = ephys.spikeSorting.MixtureModel.fit(b, df, 'verbose', true);
    results.s = ss;
    results.w = wf;
    results.b = b;
    results.model = model;
    
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
        plot(wf(:,1:(nSpikes/1e3):end), 'Color', cmap(iUnit,:));
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
sp.st = io.convertSamplesToTime(sp.ss, Fs, info.timestamps, info.fragments);
sp.cids = unique(sp.clu)';

save(fullfile(ops.root, 'sp.mat'), '-v7.3', '-struct', 'sp')