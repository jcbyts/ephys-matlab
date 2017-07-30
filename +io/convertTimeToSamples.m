function samples = convertTimeToSamples(times, adfreq, ts, fn)
% Convert timestamps into samples using sampling rate and start-time
% samples = convertTimeToSamples(times, adfreq, ts, fn)
% inputs:
%   times  [any] - vector/matrix of timestamps
%  adfreq  [1 x 1] - sampling rate
%      ts  [m x 1] - vector of recording fragment start timestamps
%      fn  [m x 1] - vector of fragment sample counts
%  ts and fn are necessary to adjust sample times for recordings that were
%  paused or are not continuous

% 20140531  jly     wrote it

if nargin < 4
    fn = round(times(end)*adfreq);
    if nargin < 3
        ts = 0;
        if nargin < 2
            help convertTimeToSamples
            samples = [];
            return
        end
    end
end

nfragments = numel(fn);
fb = ts;
fe = (ts+fn*adfreq);

samples = nan(size(times));
for ff = 1:nfragments
    idx = times >= fb(ff) & times <= fe(ff);
    samples(idx) = floor((times(idx) - fb(ff))*adfreq)+1;
end