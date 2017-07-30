function times = convertSamplesToTime(samples, adfreq, ts, fn)
% Convert samples into timestamps using sampling rate and start-time
% times = convertSamplesToTime(samples, adfreq, ts, fn)
% inputs:
%  samples [n x 1] - vector of sample indices
%  adfreq  [1 x 1] - sampling rate
%      ts  [m x 1] - vector of recording fragment start timestamps
%      fn  [m x 1] - vector of fragment sample counts
%  ts and fn are necessary to adjust sample times for recordings that were
%  paused or are not continuous

% 20140531  jly     wrote it

if nargin < 4
    fn = max(samples);
    if nargin < 3
        ts = 0;
        if nargin < 2
            help convertSamplesToTime
            times = [];
            return
        end
    end
end

nfragments = numel(fn);
fb = ts;
% fe = (ts+fn*adfreq);
sb = [1; cumsum(fn(1:end-1))];
se = cumsum(fn);
times = zeros(numel(samples),1);
for ff = 1:nfragments
    idx = samples >= sb(ff) & samples <= se(ff);
    times(idx) = ((samples(idx) - sb(ff))+1)/adfreq+fb(ff);
end
