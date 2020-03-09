function [data, timestamps, info] = getLFP(ops, overwrite, plotIt)
% GET LFP gets the lfp from a specific recording shank specified by an ops
% It will import the lfp and remove line noise if this is the first time it
% is called.
% Inputs:
%   ops@struct     - the ops struct from a single recording shank
%   overwrite@logical - whether to overwrite the import
%   plotIt@logical  - boolean. plots line noise estimation during import.
% Output:
%   data@double       - nSamples x nChannels data (in mV)
%   timestamps@double - nSamples x 1 time (in OE time)
%   info@str
% [data, timestamps, info] = getLFP(ops, showOutput)

% 2017.08.14    jly     wrote it

if nargin < 3
    plotIt = true;
    if nargin <2
        overwrite = false;
    end
end

info = io.loadEphysInfo(ops);

flfp     = fullfile(ops.root, 'lfp.dat');
flfpInfo = fullfile(ops.root, 'lfp_info.mat');

if exist(flfp, 'file') && ~overwrite
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
OE_HIGHPASS = 0.1;
CHUNK_SIZE = 1000; % seconds

useGpu = false;
% try
%     if gpuDeviceCount > 0
%         useGpu = true;
%     end
% end

% --- begin analysis

% build low-pass filter
% [b, a] = butter(1, LOW_CUTOFF/info.sampleRate*2, 'low'); don't use
% butterworh
loratio = LOW_CUTOFF/info.sampleRate*2;

% --- load raw data
tic
fprintf('Computing LFP from Raw\n')
fprintf('Loading raw data...\t')

fid  = fopen(ops.fbinary);
fout = fopen(flfp, 'w');

fseek(fid, 0, 'eof');
filesize = ftell(fid);

downStep = info.sampleRate/NEW_FS;

assert(mod(downStep,1)==0, 'Downsampling must be in integer values. Check your sampling rates')

Nchan = ops.Nchan;

fseek(fid, 0, 'bof');

nTotSamp = filesize/Nchan/2;

% Correct the distortion of the low LFP frequencies by the high pass filter
% of the amplifier (see "Artifactual origin of biphasic cortical spike-LFP 
% correlation", bioRxiv, http://dx.doi.org/10.1101/051029 for details).

% The values below come from direct measurement of the equipment, as described in the bioRxiv note
switch OE_HIGHPASS
  case 0.1
    freq  = [0.0615    0.1214    0.1845    0.2247    0.2914    0.3732    0.8186    1.1102    2.0432    3.0051 11.1815   20.7900   30.1811];
    phase = [3.0659    2.7502    2.4215    2.1768    2.0019    1.7454    1.1285    0.8774    0.5578    0.4007 0.1024    0.0389    0.0145];
  case 1
    freq  = [0.1067    0.1775    0.2308    0.2950    0.4092    0.8221    1.1241    2.0472    2.9354   11.2952  20.3804   29.6150];
    phase = [3.6443    3.2698    2.8039    2.5902    2.2015    1.4175    1.1069    0.6644    0.4840   0.1121   0.0500    0.0213];
  otherwise
    error('correctOeLfpPhase: this cutoff value isn''t supported')
end
freq(end+1) = 50;
phase(end+1) = 0;

Nsample = 2^nextpow2(min(nTotSamp, info.sampleRate*CHUNK_SIZE)); % ~500s chunks

freqX  = info.sampleRate * (1:Nsample/2-1)/Nsample; % frequency of each FFT component (to be computed below)
freqXI = freqX >= freq(1) & freqX <= freq(end); % frequencies for which phase distortion was measured
phaseDistort = zeros(Nsample/2-1,1);
phaseDistort(freqXI)  = interp1(log(freq), phase, log(freqX(freqXI)), 'pchip');
phaseDistort(~freqXI) = 0; % The rationale is that for frequencies below freq(1) the power is so low it doesn't matter if we don't correct them.

% --- low pass filter
fprintf('Preprocessing the data in %d second chunks\n', CHUNK_SIZE)
fprintf('Preprocessing Steps:\n')
fprintf('1) Correct phase distortion from OE headstage highpass at %dHz...\n', OE_HIGHPASS)
fprintf('2) Lowpass filter at %dHz...\n', LOW_CUTOFF)
fprintf('3) Downsample to %dHz...\n', NEW_FS)
fprintf('4) Reject artifacts bigger than %d SD...\n', ARTIFACT_SD)
fprintf('5) Remove line noise at %d,%d,%d Hz...\n', LINE_NOISE_FREQ)


%% Do this chunkwise
% while 1
nBlocks = ceil(nTotSamp/Nsample);
block = 1;

bufferSize = [Nchan Nsample];

while 1
fprintf('Block %d / %d\t', block, nBlocks)

% --- read in raw
dataRAW = fread(fid, bufferSize, '*int16');

% --- convert from integer to double
dataRAW = double(dataRAW');

if useGpu
    dataRAW = gpuArray(dataRAW);
end

if block == nBlocks
    Nsample = 2^nextpow2(min(size(dataRAW,1), info.sampleRate*CHUNK_SIZE)); % ~500s chunks
    
    freqX  = info.sampleRate * (1:Nsample/2-1)/Nsample; % frequency of each FFT component (to be computed below)
    freqXI = freqX >= freq(1) & freqX <= freq(end); % frequencies for which phase distortion was measured
    phaseDistort = zeros(Nsample/2-1,1);
    phaseDistort(freqXI)  = interp1(log(freq), phase, log(freqX(freqXI)), 'pchip');
    phaseDistort(~freqXI) = 0;
end

% --- loop over channels and undo the phase distortion from the OE highpass
for ch = 1:Nchan
    X = fft(dataRAW(:,ch), Nsample); % do we need to do this on the whole thing?
    X = X(2:Nsample/2);  
    X = abs(X).*exp(1i.*(angle(X) - phaseDistort)); % shift each frequency's phase by -phaseDistort
    y = ifft([0; X; 0; conj(flipud(X))]);
    dataRAW(:,ch) = y(1:size(dataRAW,1));
end

% --- low-pass filter with 0 group delay (filter both ways)  
% This is faster than filtfilt

dataRAW = iosr.dsp.sincFilter(dataRAW, loratio);
% dataRAW = filter(b,a,dataRAW);
% dataRAW = flipud(dataRAW);
% dataRAW = filter(b,a,dataRAW);
% dataRAW = flipud(dataRAW);

if useGpu
    dataRAW = gather(dataRAW);
end
    
% --- convert to volts
dataRAW = dataRAW*info.bitVolts;

% --- downsample
data = dataRAW(downStep:downStep:end,:);
clear dataRAW
% --- remove artifacts

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

% --- fit and remove line noise
iigood = ~isnan(data(:,1));

if plotIt
    f = figure('Position', [100 100 500 800]); clf
    f(2) = figure('Position', [600 100 500 800]); clf
end

lineNoiseChunk = ceil(CHUNK_SIZE/2);
if block==nBlocks
    lineNoiseChunk = 60; %ceil(size(data,1)/info.sampleRate);
end

for ch = 1:Nchan
    if plotIt
        figure(f(1)); clf
    end
    
    if useGpu
        xdata = gpuArray(data(iigood,ch));
    else
        xdata = data(iigood,ch);
    end
    
    newdata = preprocess.removeLineNoiseChunkwise(xdata, NEW_FS, LINE_NOISE_FREQ, 2,lineNoiseChunk , plotIt);
    
    if useGpu
        newdata = gather(newdata);
    end
    
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

% write to disk
dataInt = int16(data'/info.bitVolts);

fwrite(fout, dataInt, 'int16');

fprintf('[%02.2fs]\n', toc)

if block == nBlocks
    break
end

block = block + 1;


end

% --- save output
newInfo = info;
newInfo.sampleRate = NEW_FS;
newInfo.fragments  = info.fragments(:)/downStep;
newInfo.phaseCorrection = OE_HIGHPASS;
newInfo.artifacts = find(bad);
info = newInfo;

fprintf('Saving meta data...\t')
save(flfpInfo, '-v7.3', '-struct', 'newInfo')
fprintf('[%02.2fs]\n', toc)

timestamps = io.convertSamplesToTime((1:nTotSamp)', info.sampleRate, info.timestamps(:), info.fragments(:));
fclose(fid);
fclose(fout);
fprintf('Done\n')