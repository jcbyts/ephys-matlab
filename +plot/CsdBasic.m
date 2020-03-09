function stats = CsdBasic(lfp, eventTimes, lfpInfo, varargin)
% CSDBASIC computes the current source density
% Inputs:
%   LFP     [nTime x nChannel] - raw voltage traces
%   events       [nEvents x 1] - timestamps of the events
%   lfpInfo           [struct] - lfp info struct from io.getLfp
%  
% optional arguments (as argument pairs):
%   'channelDepths'  [nChannels x 1] - array of channel depths
%   'window'         [1 x 2]         - start and stop of analysis window 
%                                      (aligned to event time)
%   'plot'           logical         -  plot if (default: true)
%   'method'         string          - csd method (default: 'spline')
%
% valid csd methods:
%       'standard' - second spatial derivative
%       'step'     - stepwise inverse method (not really sure)
%       'spline'   - interpolated inverse CSD
%                                       
% 2017 jly wrote it

ip = inputParser();
ip.addParameter('window', [-100 200])
ip.addParameter('channelDepths', [])
ip.addParameter('plot', true)
ip.addParameter('method', 'spline')
ip.parse(varargin{:});

eventTimes = eventTimes(:);

% conver times to samples
ev = io.convertTimeToSamples(eventTimes, lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));

% event-triggered LFP
[sta,~, time] = pdsa.eventTriggeredAverage(lfp, ev(:), ip.Results.window);

if isempty(ip.Results.channelDepths)
    ch0 = (-32:-1)*50;
end

switch ip.Results.method
    case 'spline'
        % compute the CSD using the spline inverse method
        CSD = csd.splineCSD(sta', 'el_pos', ch0);
    case 'standard'
        CSD = csd.standardCSD(sta', 'el_pos', ch0);
    case 'step'
        CSD = csd.stepCSD(sta', 'el_pos', ch0);
    otherwise
        error('valid methods are {spline, standard, step}')
end     

% find the sink and reversal point
ix = time > 0 & time < 100; % look over time window after flash

% sink should be the minimum value
[~,id] = min(reshape(CSD(:,ix), [], 1));
% convert to indices
[depthIndex,timeIndex] = ind2sub(size(CSD(:,ix)), id);
% upsample channels to index into them
chUp   = linspace(1, numel(ch0), size(CSD,1));
depthUp= linspace(ch0(1), ch0(end), size(CSD,1));

% find reversal point
CSD_ = CSD(:,ix);
reversalPoints = findZeroCrossings(CSD_(:,timeIndex));

% output structure
stats.STA   = sta';
stats.CSD   = CSD;
stats.time  = time;
stats.depth = depthUp;
stats.chDepths = ch0;
stats.chUp  = chUp;
stats.sinkDepth = depthUp(depthIndex);
stats.sinkChannel = chUp(depthIndex);
stats.reversalPointDepth = depthUp(reversalPoints);

if ip.Results.plot % afterall, it is a plot function
    imagesc(time, ch0, CSD-mean(CSD(:))); axis xy
    colormap jet
    hold on
    plot(time, bsxfun(@plus, sta, ch0), 'Color', repmat(.1, 1, 3))
    xlim(ip.Results.window)
    plot(time([1 end]), stats.sinkDepth*[1 1], 'w--', 'Linewidth', 2)
    plot(time([1 end]), [1; 1]*stats.reversalPointDepth, 'r--', 'Linewidth', 2)
    tmp = abs(stats.reversalPointDepth - stats.sinkDepth);
    tmp = tmp + stats.sinkDepth;
    plot(time([1 end]), [1; 1]*tmp, 'r--', 'Linewidth', 2)
end


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

if nargin < 2
  mode = 0;
end

[i,~,p] = find(data); % ignore zeros in the data vector

switch sign(mode)
  case -1
    % find -ve going crossings
    ii = find(diff(sign(p))==-2);
  case 0
    % find all zero crossings
    ii = find(abs(diff(sign(p)))==2);
  case 1
    % find +ve going crossings
    ii = find(diff(sign(p))==2);
end;

i = round((i(ii)+i(ii+1))/2);
