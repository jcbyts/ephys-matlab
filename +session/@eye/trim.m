function d = trim(d,tstart,tstop,varargin)
% trim eye position signal(s)

% 2016-11-13 - Shaun L. Cloherty <s.cloherty@ieee.org>

args = varargin;
p = inputParser;

p.addParameter('debug',false,@(x) validateattributes(x,{'logical'},{'scalar'}));

p.parse(args{:});
      
args = p.Results;

t = d.t; % sample time(s)

if isempty(tstart)
  tstart = min(t);
end

if isempty(tstop)
  tstop = max(t);
end

idx = (t >= tstart) & (t <= tstop); % samples to keep...

if ~any(idx),
  return
end

d.tsample = d.tsample(idx);
d.x = d.x(idx);
d.y = d.y(idx);
d.phght = d.phght(idx);
d.pwdth = d.pwdth(idx);
