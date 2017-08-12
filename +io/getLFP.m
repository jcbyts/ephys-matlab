function [data, timestamps, info] = getLFP(ops, plotIt)
% [data, timestamps, info] = getLFP(ops, showOutput)

if nargin < 2
    plotIt = true;
end

[~, ~, info] = io.loadSession(ops.root);

flfp     = fullfile(ops.root, 'lfp.dat');
flfpInfo = fullfile(ops.root, 'lfp_info.mat');

if exist(flfp, 'file')
    fid = fopen(flfp, 'r');
    
    fseek(fid, 0, 'eof');
    filesize = ftell(fid);

    Nchan = ops.Nchan;

    fseek(fid, 0, 'bof');
    
    nTotSamp = filesize/Nchan/2;
    bufferSize = [Nchan nTotSamp];

    data = fread(fid, bufferSize, '*int16');
    
    info = load(flfpInfo);
%     info.timestamps = info.timestamps/30e3;
    
    data = double(data')*info.bitVolts;
    
    load(ops.chanMap)
    data = data(:,chanMap(connected));
    timestamps = io.convertSamplesToTime((1:nTotSamp)', info.sampleRate, info.timestamps(:), info.fragments(:));

    return
end

% --- Declare constants
LOW_CUTOFF  = 300;
ARTIFACT_SD = 5; % reject anything outside
LINE_NOISE_FREQ = [60 120 180];
NEW_FS = 1e3;
OE_HIGHPASS = 1;

% --- begin analysis

% build low-pass filter
[b, a] = butter(1, LOW_CUTOFF/info.sampleRate*2, 'low');

% --- load raw data
tic
fprintf('Computing LFP from Raw\n')
fprintf('Loading raw data...\t')

fid  = fopen(ops.fbinary);
fout = fopen(flfp, 'w');

fseek(fid, 0, 'eof');
filesize = ftell(fid);

downStep = info.sampleRate/NEW_FS;

Nchan = ops.Nchan;

fseek(fid, 0, 'bof');

nTotSamp = filesize/Nchan/2;
bufferSize = [Nchan nTotSamp];

dataRAW = fread(fid, bufferSize, '*int16');

fprintf('[%02.2fs]\n', toc)

dataRAW = dataRAW';
% --- correct phase from Open Ephys headstage (e.g., Okun 2016)
% fprintf('Correcting phase from OE headstage @ %dHz...\t', OE_HIGHPASS)
% dataRAW = double(dataRAW');
% 
% for ch = 1:Nchan
%     dataRAW(:,ch) = preprocess.correctOeLfpPhase(dataRAW(:,ch), info.sampleRate, OE_HIGHPASS);
% end
% fprintf('[%02.2fs]\n', toc)

% --- low pass filter
fprintf('Lowpass filter at %dHz...\t', LOW_CUTOFF)

dataRAW = filtfilt(b,a,double(dataRAW));

fprintf('[%02.2fs]\n', toc)

dataRAW = dataRAW*info.bitVolts;

% --- downsample
fprintf('Downsample to %dHz...\t', NEW_FS)
data = dataRAW(downStep:downStep:nTotSamp,:);
fprintf('[%02.2fs]\n', toc)

% --- remove artifacts
fprintf('Removing artifacts at %dSD...\t', ARTIFACT_SD)

% step through, calculate running standard deviation for calculating
% threshold
stepSize = 10*NEW_FS;
nSteps = ceil(size(data,1)/stepSize);
s = nan(nSteps,Nchan);
for i = 1:(nSteps-1)
    iix = (i-1)*stepSize + (1:stepSize);
    s(i,:) = nanstd(data( iix ,:));
end

thresh = ARTIFACT_SD * nanmedian(s(:));

[~, bad] = preprocess.removeChannelArtifacts(data, thresh, 1, 50, true);
data(bad,:) = nan;
fprintf('[%02.2fs]\n', toc)

% --- fit and remove line noise
fprintf('Removing line noise...\t')
iigood = ~isnan(data(:,1));

if plotIt
    f = figure('Position', [100 100 500 800]); clf
    f(2) = figure('Position', [600 100 500 800]); clf
end

for ch = 1:Nchan
    if plotIt
        figure(f(1)); clf
    end
      
    newdata = preprocess.removeLineNoiseChunkwise(data(iigood,ch), NEW_FS, LINE_NOISE_FREQ, 2, 30, plotIt);
    
    
    if plotIt
        figure(f(2)); clf
        subplot(211)
        [PxxMean1, fr, ~, ~, Pxx1] = ephys.estimatePSD_welch(zscore(data(iigood,ch)), [], [], 5e3, NEW_FS, []);
        plot(fr, 10 * log10(PxxMean1), '-k'); hold on
        [PxxMean, fr, ~, ~, Pxx] = ephys.estimatePSD_welch(zscore(newdata), [], [], 5e3, NEW_FS, []);
        plot(fr, 10 * log10(PxxMean), '-r');
        subplot(212)
        plot(fr, 10 * log10(Pxx1), '-k'); hold on
        plot(fr, 10 * log10(Pxx), '-r');
        drawnow
    end
    
    data(iigood,ch) = newdata;
    
end
fprintf('[%02.2fs]\n', toc)

% --- save output
newInfo = info;
newInfo.sampleRate = NEW_FS;
newInfo.fragments  = info.fragments(:)/downStep;
newInfo.phaseCorrection = OE_HIGHPASS;
info = newInfo;

dataInt = int16(data'/info.bitVolts);

fwrite(fout, dataInt, 'int16');

save(flfpInfo, '-v7.3', '-struct', 'newInfo')

timestamps = io.convertSamplesToTime((1:nTotSamp)', info.sampleRate, info.timestamps(:), info.fragments(:));
fclose(fid);
fclose(fout);