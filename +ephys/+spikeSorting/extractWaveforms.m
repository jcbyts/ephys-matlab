function [w,index] = extractWaveforms(x, s, win)
% Extract spike waveforms.
%   [w,index] = extractWaveforms(x, s) extracts the waveforms at times s (given in
%   samples) from the filtered signal x using a fixed window around the
%   times of the spikes. The return value w is a 3d array of size
%   length(window) x #spikes x #channels.

if nargin < 3
    win = -8:19;        % window to extract around peak
end
k = size(x, 2);     % number of channels
n = size(s, 1);     % number of spikes
m = length(win);    % length of extracted window
index = bsxfun(@plus, s, win)';

w = reshape(x(index, :), [m n k]);
