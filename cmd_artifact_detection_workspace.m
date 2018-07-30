
%%
[sessionInfo, ops, info] = io.loadSession(thisSession);

%%

nSampleTotal = sum(info.fragments);

blockSize = info.sampleRate*220; % one minute of data at a time
% blockSize = nSampleTotal; % run on full dataset

nBlocks = ceil(nSampleTotal/blockSize);

wfs  = [];
ssix = [];

thresh = 2e4;

% build filter
[b,a] = butter(3, [(300/30e3*2), (5000/30e3*2)]);

for iBlock = 1:nBlocks
    fprintf('Loading block %d/%d\n', iBlock, nBlocks)
    % blockSize = nSampleTotal;
    
    win = (iBlock - 1) * blockSize + [1 blockSize];
    win(2) = min(win(2), nSampleTotal);
    
    dataRAW = io.loadRaw(ops, win, true);
        
    
    fprintf('Filtering\n')
    dataRAW = dataRAW';
    dataHP = filter(b,a,dataRAW);
    dataHP = flipud(dataHP);
    dataHP = filter(b,a,dataHP);
    dataHP = flipud(dataHP);
    
    % common average reference
    dataHP = bsxfun(@minus, dataHP, mean(dataHP,2));
    
    % threshold based on the energy
    En = max(dataHP.^2,[],2); % maximum across channels
    
    figure(1); clf
    plot(En)
    axis tight
    xlabel('Sample #')
    ylabel('Energy')
    hold on
    
    
    if iBlock==1
        title('Set threshold for clipping artifacts')
        drawnow
        
        [x, thresh] = ginput(1);
    end
    
    plot(xlim, thresh*[1 1], 'r--')
	drawnow
    
    % detect threshold crossings
    fprintf('Thresholding\n')
    ss = find(any(En > thresh,2));
    ss(diff(ss) < 32) = [];
    spWin = -10:32;
    
    ss((ss + min(spWin)) < 1)= []; 
    ss((ss + max(spWin)) > size(dataHP,1)) = [];
    
    fprintf('Clipping Waveforms\n')
    wf = ephys.spikeSorting.extractWaveforms(dataHP, ss, spWin);
    X = permute(wf, [1 3 2]);
    sz = size(X);
    X = reshape(X, [], sz(3));
    
    wfs = [wfs X];
    ssix = [ssix; (iBlock - 1) * blockSize + ss];
end

%% Select waveforms to delete
app = WaveformSelector(ssix, wfs);

%% delete waveforms
removeInd = app.SS;

%% load on the whole thing right now
[dataRAW, timestamps] = io.loadRaw(ops, [], true);
dataRAW = dataRAW';
n = size(dataRAW,1);
nSegments = numel(removeInd);
fprintf('%d segments set to be removed\n', nSegments)

nCh = size(dataRAW,2);
dataNEW = dataRAW;

% TODO: remove from each block -- then resave
blockStart = ((1:nBlocks) - 1) * blockSize + 1;

for iSeg = 1:nSegments
    % clip out a whole second
    widx = (-16e3:32e3) + removeInd(iSeg); 
    repix = widx(widx > 0 & widx <= n);
    for i = 1:nCh
        dataNEW(repix, i) = nan;
    end
end

for i = 1:nCh
    fprintf('channel %d\n', i)
    dataNEW(:,i) = repnan(dataNEW(:,i), 'nearest');
end
    
%     figure(1); clf
%     plotidx = (-100e3:100e3) + ixRemove(iSeg); 
%     plotIx = plotidx(plotidx > 0 & plotidx <= n);
%     plot(plotIx, bsxfun(@plus, dataRAW(plotIx,:), (1:nCh)*200), 'k'); hold on
%     plot(plotIx, bsxfun(@plus, dataNEW(plotIx,:), (1:nCh)*200), 'b');
%    	drawnow
% end
i = 0;
%%
i = i + 1;
figure(1); clf
n = size(dataRAW,1);
ix = 1:n;
times = io.convertSamplesToTime(ix, 32e3, info.timestamps, info.fragments);
plot(times, dataRAW(:, i)); hold on
% plot(times, dataNEW(:, i));
plot([1; 1]*ts', ylim'*ones(1, numel(ts)), 'k')

plot([1; 1]*[nats.trial.start], ylim'*ones(1, numel(nats.trial)), 'r')
plot([1; 1]*[hart.trial.start], ylim'*ones(1, numel(hart.trial)), 'g')
plot([1; 1]*[csd.trial.start], ylim'*ones(1, numel(csd.trial)), 'b')
%%
figure(2); clf
x = dataRAW(10e4+ (1:10e4),2);
fs = 32e3;

nwin = 5e3;
wind = hanning(nwin);
nlap = nwin-10;
freq = 0:2:256;

[S, f, t] = spectrogram(x, nwin, nlap, freq, fs, 'yaxis'); %, 8, freq, fs,  'psd', 'yaxis');
% view(0,0)
% spectrogram(x2,[],8,freq,fs)
% figure(2); clf
% imagesc(log(abs(S))')

% x = dataRAW(1:10e4,1);
%%
clf
imagesc(log(abs(S(1:128,:))))
axis xy
%%

PDS = io.getPds(thisSession);

ts = cell2mat(cellfun(@(x) cellfun(@(y) x.PTB2OE(y.timing.flipTimes(1)), x.data)', PDS, 'uniformOutput', false));

% plot([1; 1]*removeInd', ylim'*ones(1,nSegments), 'k')
%%
iix = io.findPDScontainingStimModule(PDS, 'natImgBackground');

nats = session.natImgBackground(PDS);
hart = session.hartleyFF(PDS);
csd  = session.csdFlash(PDS);

%%

nlEn = abs(dataRAW(2:end-1,:).* (dataRAW(2:end-1,:)-dataRAW(1:end-2,:)) .*dataRAW(3:end,:));
en = dataRAW(2:end-1,:).^2;

ch = 0;
%%
ch = ch + 1;
clf

plot(zscore(nlEn(:,ch))+15)

hold on
plot(zscore(en(:,ch))+4)

plot(zscore(dataRAW(:,ch)));

%%
clf
ch = ch+1;
plot(en(:,ch), nlEn(:,ch), '.')
hold on
mu = mean([en(:,ch), nlEn(:,ch)]);
C = cov([en(:,ch), nlEn(:,ch)]);
plotellipse(mu, C, 40, 'Linewidth', 4)


