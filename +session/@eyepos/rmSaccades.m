function [d,idx] = rmSaccades(d,varargin)
% Remove saccades from eye position signal(s)...
%
% Available arguments include:
%
%   dt     - padding applied around saccade segments (both before and
%            after the saccade; default: 0.010s)
%   interp - fill saccade segments in eye position signals (default: true)
%   sargs  - cell array of arguments passed through to @eye.findSaccades()
%   debug  - show debugging output (default: false)
%
% Arguments controlling saccade detection can be passed through to
% @eye.findSaccades() via the 'sargs' argument, e.g.,
%
%   nosac = d.rmSaccades('dt',0.005,'sargs',{'order',32,'Wn',[0.1,0.16]});
%
% will call d.findSaccades() as follows:
%
%   [~,~,ix] = d.findSaccades('order',32,'Wn',[0.1,0.16]);
%
% saccades with then be extended by 5ms (before and after each saccade)
% and removed from the eye position signals
%
% See also: eye.findSaccades.

% 2016-11-16 - Shaun L. Cloherty <s.cloherty@ieee.org>

args = varargin;
p = inputParser;
p.KeepUnmatched = true;

p.addParameter('dt',0.010,@(x) validateattributes(x,{'numeric'},{'scalar'}));
p.addParameter('interp',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

p.addParameter('debug',true,@(x) validateattributes(x,{'logical'},{'scalar'}));

% argments passed through for saccade detection...
p.addParameter('sargs',{},@iscell);

p.parse(args{:});

args = p.Results;

% find saccades
[~,~,idx] = d.findSaccades(args.sargs{:});

if isempty(idx)
    return
end

% extend 'window' by args.dt before and after each saccade
n = ceil(args.dt./d.dt);

idx(:,1) = max(idx(:,1)-n,1);
idx(:,2) = min(idx(:,2)+n,length(d.t));

% remove saccades from eye position data
for ii = 1:size(idx,1)
    for f = {'x','y'}
        x = d.(f{1})(idx(ii,:)); % 1x2 vector
        
        %     % handle saccades that don't 'end'...
        %     if idx(ii,2) == length(d.t)
        %       x(2) = x(1);
        %     end
        
        d.(f{1})(idx(ii,2):end) = d.(f{1})(idx(ii,2):end) - diff(x);
        
        if args.interp
            d.(f{1})(idx(ii,1):idx(ii,2)) = d.(f{1})(idx(ii,1)); % interpolate...
        else
            d.(f{1})(idx(ii,1):idx(ii,2)) = NaN;
        end
    end
end

if args.debug
    figure;
    ah(1) = subplot(2,1,1);
    plot(d.t,d.x);
    hold on
    
    ah(2) = subplot(2,1,2);
    plot(d.t,d.y);
    hold on
    
    for ii = 1:2.
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
    
    set(ah,'Layer','top');
end
