function [PxxMean, fr, totalTime, PxxStd, Pxx] = estimatePSD_welch(x, rois, w, nFFT, Fs, overlap)
% Estimate power spectral density for sub-indexed segments of data (e.g. trials)
%
% [PxxMean, fr, totalTime, PxxStd, Pxx] = estimatePSD_welch(x, rois, w, nFFT, Fs, overlap)
%
% Input
%   x: 1D real-valued signal
%   rois: (nTrial x 2) region of interests in units of integer indices.
%	  Each row consists of (start idx, end idx).
%	  If either is NaN, that segment is ignored.
%	  If empty, the entire signal is considered.
%   w: (nFFT x 1) window function (put [] for no windowing)
%   nFFT: (1) frequency resolution
%   Fs: (1) sampling rate
%   overlap: (1) # of samples to overlap when shifting (Welch's method)
%	     set overlap = [] for no overlap (Barlett's method)
%
% Output
%   Pxx: estimated power spectrum (take 10 * log10(Pxx) for dB)
%   fr: frequency range
%   totalTime: sum of sub-indexed regions
%   PxxStd: standard deviation est.
%   Pxx: Raw power signal over all DFT segments
%
% Usage
%   z = zscore(lfp_data(:, 2));
%   nFFT = 2^13; Fs = lfp_info.adfreq; w = hann(nFFT); overlap = nFFT/2;
%   rois = [timing.motionon, timing.choice];
%   [Pxx, fr, totalTime] = estimatePSD_welch(z, rois, w, nFFT, Fs, overlap);
%   plot(fr, 10 * log10(Pxx), '-r')

% check input
assert(min(size(x)) == 1, 'input signal must be 1D');
if isempty(rois); rois = [1, numel(x)]; end
assert(size(rois, 2) == 2, 'Region of interest indices must be (nTrial x 2)');
assert(numel(nFFT) == 1, 'nFFT is scalar (pref 2^n)');

% Frequency range
fr = ((1:nFFT/2+1)-1)*Fs/nFFT;

fftCount = 1;
totalTime = 0;
nnzTotal = 0;
Pxx = nan(nFFT, size(rois, 1));
for kTrial = 1:size(rois, 1)
    start_idx = rois(kTrial, 1);
    end_idx = rois(kTrial, 2);
    if isnan(start_idx) || isnan(end_idx)
	continue;
    end
    totalTime = totalTime + end_idx - start_idx;

    if ~isempty(overlap)
	zz = buffer(x(start_idx:end_idx), nFFT, overlap);
    else
	zz = buffer(x(start_idx:end_idx), nFFT);
    end
    nnzTotal = nnzTotal + nnz(zz); % buffer adds quite a number of trailing 0's

    if ~isempty(w)
	zz = bsxfun(@times, zz, w);
    end

    zz = abs(fft(zz) / nFFT).^2;
    Pxx(:, fftCount:fftCount+size(zz,2)-1) = zz;
    fftCount = fftCount + size(zz,2);
end
Pxx = 2 * Pxx(1:(nFFT/2+1), 1:fftCount-1); % remove unncessary parts
if ~isempty(w)
    Pxx = Pxx / mean(w.^2); % normalization for proper power 
    % TODO (not quite accurate because of buffer)
end

% normalize to make it a 'density' (integrates to nFFT * power if using Fs = 1)
PxxStd = std(Pxx, [], 2) * nFFT / Fs;
PxxMean = mean(Pxx, 2) * nFFT / Fs;

% compensate for the missing values when using `buffer`
PxxStd = PxxStd * size(Pxx, 2) / nnzTotal * nFFT;
PxxMean = PxxMean * size(Pxx, 2) / nnzTotal * nFFT;
