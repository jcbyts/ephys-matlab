function [data, timestamps, info] = getEyeData(thisSession, PDS, varargin)

ip = inputParser();
ip.addOptional('overwrite', false)
ip.addOptional('eyetracker', [])
ip.parse(varargin{:});

if nargin < 2
    PDS = io.getPds(thisSession);
end

%assert(istable(thisSession), 'getEyeData: first argument must be a meta table entry')

eyetracker = ip.Results.eyetracker;
if isempty(eyetracker)
    
    eyetracker = 'none';
    
    % figure out which eyetracker to use
    hasEyelink   = cellfun(@(x) isfield(x.initialParametersMerged, 'eyelink'), PDS);
    hasArrington = cellfun(@(x) isfield(x.initialParametersMerged, 'arrington'), PDS);
    
	useArrington = cellfun(@(x) x.initialParametersMerged.arrington.useAsEyepos, PDS(hasArrington));
    useEyelink   = cellfun(@(x) x.initialParametersMerged.eyelink.useAsEyepos, PDS(hasEyelink));
    
    if sum(useEyelink) > sum(useArrington)
        eyetracker = 'eyelink';
    end
    
    if sum(useArrington) > sum(useEyelink)
        eyetracker = 'arrington';
    end
    
end

switch eyetracker
    case 'eyelink'
        [data, timestamps, info] = io.getEdf(thisSession,PDS,ip.Results.overwrite);
    case 'arrington'
        [data, timestamps, info] = io.getVpx(thisSession,PDS,ip.Results.overwrite);
    otherwise
        data = [];
        timestamps = [];
        info = [];
        
end
    
    