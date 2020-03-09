function [tstart,tend,idx,vel] = findPursuit(d,varargin)
% Find pursuit eye movements in eye position signals...
%
% Available arguments include:
%
%   order    - low-pass digital differentiating filter order (default: 128)
%   Wn       - low-pass filter corner freq, and transition band as a
%              percentage of the Nyquist frequency (default: [0.005,0.01])
%   dt       - padding applied around saccade segments (both before and
%              after the saccade; default: 0.010s)
%   thresh   - velocity threshold (default: 5 deg./s)
%   duration - minimum pursuit duration (default: 0.050s)
%   sargs    - cell array of arguments passed through to @eye.findSaccades()
%   debug    - show debugging output (default: false)
%
% See also: eye.findSaccades, eye.rmSaccades.

% using lpfirdd between 40 and 60Hz works pretty well...

% 1. median filter, n = 32
% 2. resample @ 1kHz

% 3. lpfirdd 50Hz (80Hz) --> saccades
% 4. remove saccades
% 5. lpfirdd 2.5Hz (5.0Hz) --> pursuit

% 2016-10-31 - Shaun L. Cloherty <s.cloherty@ieee.org>

args = varargin;
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('order',128,@(x) validateattributes(x,{'numeric'},{'scalar','even'})); % median filter order
p.addParameter('Wn',[2.5,5.0]/500,@(x) validateattributes(x,{'numeric'},{'vector'})); % filter transition band (percentage of Nyquist frequency)

p.addParameter('dt',0.010,@(x) validateattributes(x,{'numeric'},{'scalar'}));

p.addParameter('thresh',3.0,@(x) validateattributes(x,{'numeric'},{'scalar','positive'})); % pursuit velocity threshold
p.addParameter('duration',0.050,@(x) validateattributes(x,{'numeric'},{'scalar','positive'})); % minimum pursuit duration...?

p.addParameter('debug',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

% argments passed through for saccade detection...
p.addParameter('sargs',{},@iscell);

p.parse(args{:});

args = p.Results;

% % find saccades...
% [~,~,ix] = d.findSaccades(args.sargs);
%
% if 0 & ~isempty(ix),
%   % extend 'window' by args.dt before and after each saccade
%   n = ceil(args.dt./d.dt);
%
%   ix(:,1) = max(ix(:,1)-n,1);
%   ix(:,2) = min(ix(:,2)+n,length(d.t));
% end

% remove saccades
[d,idx] = d.rmSaccades('dt',args.dt,'interp',true,'debug',args.debug,'sargs',args.sargs);

% calculate velocity... low-pass FIR digital differentiator
[~,~,~,vel] = d.findSaccades('order',args.order,'Wn',args.Wn);

% now we're looking for pursuit...
speed = hypot(vel.x,vel.y); % direction and speed

if args.debug
    figure;
    
    ah(1) = subplot(3,1,1);
    plot(d.t,vel.x); hold on;
    
    ylabel('H. Vel. (deg./s)');
    
    ah(2) = subplot(3,1,2);
    plot(d.t,vel.y); hold on;
    
    ylabel('V. Vel. (deg./s)');
    
    ah(3) = subplot(3,1,3);
    plot(d.t,speed); hold on;
    
    xlabel('Time (s)');
    ylabel('Eye speed (deg./s)');
    
    if ~isempty(idx)
        for ii = 1:3.
            set(gcf, 'currentaxes', ah(ii));
            
            % show saccades (shaded)
            for jj = 1:size(idx,1)
                x = kron(d.t(idx(jj,:))',[1,1]);
                y = [ylim, fliplr(ylim)];
                fh = fill(x,y,[1.0,0.90,0.90]); % red
                set(fh,'EdgeColor','none','ZData',-1*ones(size(x)));
            end
            
            xlim([min(d.t),max(d.t)]);
        end
    end
    
    % show threshold
    plot(xlim,args.thresh*[1,1],'r--');
    
    set(ah,'Layer','top');
end

% +ve going zero crossings
idx = findZeroCrossings(speed - args.thresh,1);

n = ceil(args.duration*d.fs); % min. duration (in samples)

% 'debounce'... using a moving minimum
mmn = movingmin(speed,n);
idx(mmn(idx) < args.thresh) = [];

if isempty(idx)
    tstart = [];
    tend = [];
    return
end

idx = kron(idx,ones(1,2)); % pursuit indicies (start, end)
for ii = 1:size(idx,1)
    %   t0 = max(idx(ii,2)-n,1);
    t0 = idx(ii,1);
    %   t1 = min(idx(ii,2)+n,length(t));
    t1 = length(d.t);
    
    % window speed signal (and apply threshold)
    tmp = speed(t0:t1) - args.thresh;
    
    if args.debug
        for jj = 1:3
            set(gcf, 'currentaxes', ah(jj));
            
            xx = d.t([t0,t1,t1,t0]);
            yy = kron(get(gca,'YLim'),[1,1]);
            fh(jj) = fill(xx(:),yy(:),[0.9,1.0,0.9]); % green
            set(fh,'EdgeColor','none','ZData',-2*ones(size(xx)));
        end
    end
    
    % search forward in time to find pursuit end
    k = findZeroCrossings(tmp,-1);
    %   if ~isempty(k)
    %     k(k < n+1) = [];
    %   end
    if isempty(k)
        % no pursuit end...!?
        k = (t1-t0)+1;
    end
    idx(ii,2) = min(k) + t0 - 1; % index of blink end
    
    if args.debug && ~isnan(idx(ii,2))
        xx = d.t(idx(ii,[1,2,2,1]));
        arrayfun(@(h) set(h,'XData',xx),fh);
    end
end

tstart = d.t(idx(:,1));
tend = d.t(idx(:,2));
