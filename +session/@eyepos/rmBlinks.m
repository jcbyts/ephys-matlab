function d = rmBlinks(d,ix,varargin)
% remove blinks from eye position signal(s)...

% 2016-11-04 - Shaun L. Cloherty <s.cloherty@ieee.org>

args = varargin;
p = inputParser;
p.KeepUnmatched = true;

p.addParameter('interp',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
p.addParameter('dt',0.020,@(x) validateattributes(x,{'numeric'},{'scalar'}));

p.addParameter('debug',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

p.parse(args{:});

args = p.Results;

% % find blinks...
% tmp = p.Unmatched;
% tmp.debug = args.debug;
% [~,~,ix] = d.findBlinks(tmp);

if isempty(ix)
    return
end

% % extend 'window' by args.dt before and after each blink
% n = ceil(args.dt./median(diff(d.t)));
%
% ix(:,1) = max(ix(:,1)-n,1);
% ix(:,2) = min(ix(:,2)+n,length(d.t));

% remove blinks for eye position data
for ii = 1:size(ix,1)
    for f = {'x','y'}
        if args.interp
            x = d.(f{1})(ix(ii,:)); % 1x2 vector
            
            % handle blinks that don't "end" within the time window available...
            if ix(ii,2) == length(d.t)
                x(2) = x(1);
            end
            
            d.(f{1})(ix(ii,1):ix(ii,2)) = interp1( ... % x(ts)
                d.t(ix(ii,:)), ... % t
                x, ... % x(t)
                d.t(ix(ii,1):ix(ii,2))); % ts
        else
            d.(f{1})(ix(ii,1):ix(ii,2)) = NaN;
        end
    end
end

if args.debug
    figure;
    subplot(2,1,1);
    plot(d.t,d.x);
    
    subplot(2,1,2);
    plot(d.t,d.y);
end
