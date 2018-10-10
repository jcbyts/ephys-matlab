function [tstart,tend,idx,vel,acc] = findSaccades(d,varargin)
% Find saccades in eye position signals...
%
% Available arguments include:
%
%   order     - low-pass digital differentiating filter order (default: 32)
%   Wn        - low-pass filter corner freq, and transition band as a
%               percentage of the Nyquist frequency (default: [0.1,0.16])
%   accthresh - acceleration threshold (default: 2e4 deg./s^2)
%   velthresh - velocity threshold (default: 10 deg./s)
%   velpeak   - minimum peak velocity (default: 10 deg./s)
%   isi       - minimum inter-saccade interval (default: 0.050s)
%   debug     - show debugging output (default: false)
%
% See also: lpfirdd.

% 2016-11-13 - Shaun L. Cloherty <s.cloherty@ieee.org>

args = varargin;
p = inputParser;

p.addParameter('order',32,@(x) validateattributes(x,{'numeric'},{'scalar','even'})); % low-pass filter/differentiator order
p.addParameter('Wn',[0.1,0.16],@(x) validateattributes(x,{'numeric'},{'vector'})); % filter transition band (percentage of Nyquist frequency)

p.addParameter('accthresh',5e3,@(x) validateattributes(x,{'numeric'},{'scalar','positive'})); % accel. threshold for detection (deg./s^2) (was 2e4)
p.addParameter('velthresh',10,@(x) validateattributes(x,{'numeric'},{'scalar','positive'})); % velocity threshold (deg./s)
p.addParameter('velpeak',10,@(x) validateattributes(x,{'numeric'},{'scalar','positive'})); % min. peak velocity (deg./s)

p.addParameter('isi',0.050,@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
p.addParameter('dt',0.075,@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));

p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

p.parse(args{:});

args = p.Results;

% low-pass FIR digital differentiator coefficients
N = round(args.order/2);
copt = lpfirdd(N, args.Wn(1), args.Wn(2), 1 ,0)';
coeffs = [fliplr(copt), 0, -copt];

% get gaze position data...
t = d.t;
pos = struct('x',d.x,'y',d.y);

fs = d.fs; % sampling freq. (samples/s)

% horiz. (x) and vert. (y) velocities
vel = structfun(@(x) fs*locfilt(coeffs,1,x),pos,'UniformOutput',false);

% scalar eye speed
speed = hypot(vel.x,vel.y);

% estimate baseline (e.g., pursuit) speed using a moving average...
a = 1;
b = ones(1,2*N+1)./(2*N+1);
baseline = locfilt(b,a,speed);
baseline = circshift(baseline,N); % delay by half filter length

% scalar eye acceleration
accel = fs*locfilt(coeffs,1,speed);

% find saccades... i.e., -ve going zero crossings in the scalar eye
% acceleration signal. Note, quantizing the acceleration trace
% (i.e., dividing by 1e4) improves noise immunity
idx = findZeroCrossings(fix(accel)./args.accthresh,-1);

% ignore saccades with peak speed less than args.velpeak
idx(speed(idx) < args.velpeak) = [];

% ignore saccades too close to the start or end of the recording
idx(idx < args.order+args.dt*fs) = [];
idx(idx > length(t)-args.dt*fs) = [];

% sanity check...
idx(speed(idx) < baseline(idx)+args.velthresh) = [];

acc = accel;

if isempty(idx)
    tstart = [];
    tend = [];
    return
end

% sanity check...
k = bsxfun(@minus,idx,[1:args.dt*fs]-1);
idx(~any(sign(fix(accel(k')./args.accthresh)) > 0,1)) = [];

n = ceil(args.isi*fs); % samples

% test for minimum inter-saccade interval (args.isi) violations?
if numel(idx) > 1
    ii = find(diff(idx) < n,1);
    while ~isempty(ii)
        % keep larger/faster of the two...
        [~,jj] = min(speed(idx(ii+[0,1])));
        idx(ii+jj-1) = [];
        ii = find(diff(idx) < n,1);
    end
end

if args.debug
    figure;
    subplot(2,1,1);
    plot(t,speed);
    hold on;
    
    subplot(2,1,2);
    plot(t,accel./args.accthresh);
    hold on;
    
    plot(t,fix(accel./args.accthresh));
end

if isempty(idx)
    tstart = [];
    tend = [];
    return
end

n = ceil(args.dt*d.fs); % samples

idx = kron(idx,ones(1,3)); % saccade indicies (start, mid, end)
for ii = 1:size(idx,1)
    t0 = max(idx(ii,2)-n,1);
    t1 = min(idx(ii,2)+n,length(t));
    
    % window speed signal (and apply threshold)
    tmp = speed(t0:t1) - args.velthresh;
    
    if args.debug
        for jj = 1:2
            subplot(2,1,jj);
            
            xx = t([t0,t1,t1,t0]);
            yy = kron(get(gca,'YLim'),[1,1]);
            fh(jj) = fill(xx(:),yy(:),zeros(1,3));
            set(fh,'ZData',-1*ones(size(xx)),'FaceAlpha',0.1,'LineStyle','none');
        end
        
        subplot(2,1,1);
        plot(t([t0,t1]),args.velthresh*[1,1],'k--');
        plot(t([t0:t1]),baseline(t0:t1)+args.velthresh,'k-');
    end
    
    % now search back in time to find saccade start, accounting
    % for any baseline/pursuit speed
    k = findZeroCrossings(tmp-baseline(t0:t1),1);
    if isempty(k)
        continue
    end
        
    k(k > n+1) = [];
    %   if ~isempty(k)
    idx(ii,1) = max(k) + t0 - 1; % index of saccade start
    %   end
    
    % search forward in time to find saccade end
    k = findZeroCrossings(tmp-baseline(idx(ii,1)),-1);
    k(k < n+1) = [];
    if ~isempty(k)
        idx(ii,3) = min(k) + t0 - 1; % index of saccade end
    end
    
    if args.debug
        xx = t(idx(ii,[1,3,3,1]));
        arrayfun(@(h) set(h,'XData',xx),fh);
    end
    
    % the code below attempts to refine our estimate of saccade start
    % and end by computing eye velocity in the direction of the saccade,
    % i.e., by projecting x and y velocity components onto the saccade
    % vector.
    %
    % if this fails, we should go with our current estimate...
    
    % determine saccade direction vector...
    v = [diff(pos.x(idx(ii,[1,3]))); diff(pos.y(idx(ii,[1,3])))]; % FIXME: should dy be -ve?
    v = v./norm(v); % unit vector
    
    % project velocity onto the saccade vector
    t0 = max(idx(ii,1)-(2*N+1),1);
    t1 = min(idx(ii,3)+(2*N+1),length(t));
    
    sacvel = [vel.x(t0:t1), vel.y(t0:t1)]*v;
    
    % apply threshold
    tmp = sacvel - args.velthresh;
    
    % estimate baseline/pursuit speed at saccade start
    sacvel_ = nanmean(sacvel(1:2*N+1));
    
    % refine estimate of saccade start and end
    k = findZeroCrossings(tmp-sacvel_,1);
    k(k+t0-1 > idx(ii,2)) = [];
    if ~isempty(k)
        idx(ii,1) = max(k) + t0 - 1; % index of saccade start
    else
        % go with our current estimate...
        if args.debug
            warning('Failed to refine estimate of start time for saccade at t = %3f.',d.t(idx(ii,2)));
        end
    end
    
    % estimate baseline/pursuit speed at saccade end
    sacvel_ = nanmean(sacvel(end-(2*N+1):end));
    
    k = findZeroCrossings(tmp-sacvel_,-1);
    k(k+t0-1 < idx(ii,2)) = [];
    if ~isempty(k)
        idx(ii,3) = min(k) + t0 - 1; % index of saccade end
    else
        % go with our current estimate...
        if args.debug
            warning('Failed to refine estimate of end time for saccade at t = %3f.',d.t(idx(ii,2)));
        end
    end
    
    if args.debug
        xx = t(idx(ii,[1,3,3,1]));
        arrayfun(@(h) set(h,'XData',xx),fh);
    end
    
end

tstart = t(idx(:,1));
tend = t(idx(:,3));
idx(:,2) = [];

% vel = speed;
acc = accel;

    function y = locfilt(b,a,x)
        y = filter(b,a,x);
        y(1:2*N) = NaN;
        y = circshift(y,-N);
    end

end