function d = resample(d,varargin)
% resample eye position signal(s) using cubic interpolation
%
% available options are:
%
%   fs    - desired sample rate (default: 1kHz)
%   order - optional median (pre-)filter order
%           (default: 1, no median filtering)
%   debug - show debugging output (default: false)

% 2016-11-13 - Shaun L. Cloherty <s.cloherty@ieee.org>

args = varargin;
p = inputParser;
p.addParameter('fs',1e3,@(x) validateattributes(x,{'numeric'},{'scalar'}));

p.addParameter('order',1,@(x) validateattributes(x,{'numeric'},{'scalar'})); % median filter order

% p.addParameter('blksz',1e4,@(x) validateattributes(x,{'numeric'},{'scalar','positive'})); % median filter order

p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

p.parse(args{:});
      
args = p.Results;

t = d.tsample; % sample time(s)

pnames = {'x','y','pwdth','phght'}; % properties to resample

% concatenate for ease of filtering etc.
tmp = [d.x, d.y, d.pwdth, d.phght];
  
if args.debug,
  for ii = 1:size(tmp,2),
    subplot(size(tmp,2),1,ii);
    plot(t,tmp(:,ii),'-');
    ylabel(pnames{ii});
    hold on
  end
end

% optional median (pre-)filter
if args.order > 1,
  tmp = medfilt1(tmp,args.order,[],'truncate');
end

if args.debug,
  for ii = 1:size(tmp,2),
    subplot(size(tmp,2),1,ii);
    plot(t,tmp(:,ii),'--');
  end
end

% resample with cubic interpolation
tstart = t(1); 
tend = t(end);

N = 20;
t = [flipud(tstart-(t(1+[1:N])-tstart)); t; flipud(tend-(t(end-N:end-1)-tend))];

tmp = [repmat(tmp(1,:),N,1); tmp; repmat(tmp(end,:),N,1)];

[tmp,t] = resample(tmp,t,args.fs,'spline');

ix = find((t >= tstart) & (t <= tend)); % samples to keep
 
t = t(ix);
tmp = tmp(ix,:);

if args.debug,
  for ii = 1:size(tmp,2),
    subplot(size(tmp,2),1,ii);
    plot(t,tmp(:,ii),'k-');
    hold off
  end
end

d.tsample = t; % < d.t = d.tsample + d.toffset
d.toffset = d.toffset + (t(1) - tstart); % is this appropriate?
for ii = 1:size(tmp,2),
  d.(pnames{ii}) = tmp(:,ii);
end
