function [tstart,tend,idx,delta] = findBlinks(d,varargin)
% find blinks in pupil size (actually area) signal...

% 2016-11-03 - Shaun L. Cloherty <s.cloherty@ieee.org>

args = varargin;
p = inputParser;

p.addParameter('order',8,@(x) validateattributes(x,{'numeric'},{'scalar','even'})); % low-pass filter/differentiator order
p.addParameter('thresh',0.95,@(x) validateattributes(x,{'numeric'},{'scalar','positive'})); % pupil area threshold (as percentage of baseline)
p.addParameter('duration',0.050,@(x) validateattributes(x,{'numeric'},{'scalar','positive'})); % minimum blink duration...?

% p.addParameter('start',[],@(x) validateattributes(x,{'numeric'},{'scalar'})); % start time (absolute?, in seconds?)
% p.addParameter('stop',[],@(x) validateattributes(x,{'numeric'},{'scalar'}));  % stop time

p.addParameter('debug',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

p.parse(args{:});
      
args = p.Results;

% get pupil area
t = d.t;
par = d.parea;

fs = 1./median(diff(t)); % sampling freq. (samples/s)

% compute baseline pupil area...
a = 1;
b = ones(1,args.order)./args.order;
baseline = filter(b,a,par);
baseline(1:args.order) = NaN;
% baseline = circshift(baseline,-round(args.order/2));

% find blinks... i.e., -ve going threshold crossings
% in the pupil area
idx = findZeroCrossings(par-baseline*args.thresh,-1); % area

n = ceil(args.duration*fs); % min. duration (in samples)

% 'debounce'... using a moving maximum
mmx = movingmax(par,n);
idx(mmx(idx) > baseline(idx)*args.thresh) = [];

% ignore blinks too close to the start or end of the recording
% idx(idx < 2*args.order+args.duration*fs) = [];
idx(idx < args.order) = [];
% idx(idx > length(t)-args.duration*fs) = [];

% test for minimum blink duration (?) violations...
%
% FIXME: here we really want a constraint on blink rate or interval...
if numel(idx) > 1
  idx(find(diff(idx) < n)+1) = [];
end

if args.debug
  figure;
  plot(t,par);
  hold on;
  
  plot(t,baseline,'k--');
  plot(t,baseline*args.thresh,'k-');
  
  plot(t,mmx,'m-');
end

if isempty(idx)
  tstart = [];
  tend = [];
  delta = [];
  return
end

delta = NaN(size(idx)); % change in pupil area (normalized)

idx = kron(idx,ones(1,2)); % blink indicies (start, end)
for ii = 1:size(idx,1),
%   t0 = max(idx(ii,1)-n,1);
  t0 = idx(ii,1);
%   t1 = min(idx(ii,2)+n,length(t));
  t1 = length(t);
  
  % window pupil area signal (and apply threshold?)
  tmp = par(t0:t1); % + args.thresh; % area
  
  if args.debug,
    plot(t(t0:t1),par(t0:t1),'k:');

    xx = t([t0,t1,t1,t0]);
    yy = kron(get(gca,'YLim'),[1,1]);
    fh = fill(xx(:),yy(:),zeros(1,3));
    set(fh,'ZData',-1*ones(size(xx)),'FaceAlpha',0.1,'LineStyle','none');
    
    plot(t([t0,t1]),baseline(idx(ii,1))*args.thresh*[1,1],'r-');
  end
  
%   % now search back in time to find blink start
%   k = findZeroCrossings(tmp-baseline(t0:t1)*args.thresh,-1); % area
%   k(k > n+1) = [];
%   idx(ii,1) = max(k) + t0 - 1; % index of blink start
    
  % search forward in time to find blink end
  k = findZeroCrossings(tmp-baseline(idx(ii,1))*args.thresh,1); % area
  if ~isempty(k)
    k(k < n+1) = [];
  end
  if isempty(k)
    % no blink end...!?
    k = (t1-t0)+1;
  end
  idx(ii,2) = min(k) + t0 - 1; % index of blink end
  
%   delta(ii) = range(par(idx(ii,1):idx(ii,2)))./range(par);
  delta(ii) = min(par(idx(ii,1):idx(ii,2)))./baseline(idx(ii,1)); 

  if args.debug && ~isnan(idx(ii,2)),
    xx = t(idx(ii,[1,2,2,1]));
    arrayfun(@(h) set(h,'XData',xx),fh);
  end

end

% remove false positives... e.g., due to noise etc.
if ~isempty(idx)
  ii = diff(idx,1,2) < n; % min. duration
  ii = ii | delta > 0.05; % min. modulation
    
  idx(ii,:) = [];
  delta(ii,:) = [];  
end

% % remove duplicates and/or other false positives...
% ii = 2;
% while ii <= size(idx,1)
%   if idx(ii,1) < idx(ii-1,2),
%     idx(ii,:) = [];
%     continue
%   end
%   ii = ii + 1;
% end

tstart = t(idx(:,1));
tend = t(idx(:,2));

% % nasty hack! - NaN indicated a blink that we cannot determine the
% % end of... e.g., one that ends after the end of the available data
% ii = (idx(:,2) == length(t)) & (baseline(idx(:,2)) < baseline(idx(:,1))*args.thresh);
% idx(ii,2) = NaN;
