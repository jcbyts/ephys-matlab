function i = findZeroCrossings(data, mode)
%FINDZEROCROSSINGS Find zero crossing points.
%   I = FINDZEROCROSSINGS(DATA,MODE) returns the indicies into the supplied
%   DATA vector, corresponding to the zero crossings.
%
%   MODE specifies the type of crossing required:
%     MODE < 0 - results in indicies for the -ve going zero crossings,
%     MODE = 0 - results in indicies for ALL zero crossings (default), and
%     MODE > 0 - results in indicies for the +ve going zero crossings.

% $Id: findZeroCrossings.m,v 1.1 2008-07-21 23:31:50 shaunc Exp $

if nargin < 2,
  mode = 0;
end

[i,j,p] = find(data); % ignore zeros in the data vector

switch sign(mode),
  case -1,
    % find -ve going crossings
    ii = find(diff(sign(p))==-2);
  case 0,
    % find all zero crossings
    ii = find(abs(diff(sign(p)))==2);
  case 1,
    % find +ve going crossings
    ii = find(diff(sign(p))==2);
end;

i = round((i(ii)+i(ii+1))/2);