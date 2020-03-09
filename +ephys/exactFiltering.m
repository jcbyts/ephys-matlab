function x = exactFiltering(x, stopFreqRange)
% Filter 1D signal x using FFT/IFFT
% x = exactFiltering(x, stopFreqRange)
%
% x: 1D signal
% stopFreqRange: (2 x n) frequency ranges to filter out in units of 
%                normalized Nyquist frequency.
%                e.g. If 1000 Hz sampling rate, and want to low pass below 200Hz
%                use [200/1000, 500/1000] == [0.2, 0.5].
%                High pass above 200 Hz would be [0 0.2]

if size(stopFreqRange, 1) == 1
	stopFreqRange = stopFreqRange';
end

x = fft(x);
N = numel(x);
f = ((1:N) - 1) / N;
pidx = false(size(f));
for k = 1:size(stopFreqRange, 2)
	pidx = pidx | (f >= stopFreqRange(1,k) & f < stopFreqRange(2,k)) | (f > (1 - stopFreqRange(2,k)) & f <= (1 - stopFreqRange(1,k)));
end

x(pidx) = 0;
x = ifft(x);
