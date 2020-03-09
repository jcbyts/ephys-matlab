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
        bad = data(3,:)==0;
        data(:,bad) = nan;
        
    case 'arrington'
        [data, timestamps, info] = io.getVpx(thisSession,PDS,ip.Results.overwrite);
    otherwise
        data = [];
        timestamps = [];
        info = [];
        
end

if mean(isnan(timestamps)) > .9 % reading from the eyetracker file failed. data may be corrupt. try to read from PDS
    warning('eyetracker datafile missing data. loading from PDS\n')
    trial = pds.getPdsTrialData(PDS);

    eyeIdx = trial(1).eyelink.eyeIdx;
    useRaw = trial(1).eyelink.useRawData;
    ctr = trial(1).display.ctr(1:2);

    % --- Get eye position from the eyelink buffer
    sampleNames = trial(1).eyelink.sampleIds;
    timestamps = [];
    data = [];

    for iTrial = 1:numel(trial)
        cm = trial(iTrial).eyelink.calibration_matrix;
        
        if useRaw
            ixEyePos = (cellfun(@(x) any(strfind(x, 'EyeRaw')), sampleNames)); % index to find Eye position
            ptbElXY   = trial(iTrial).eyelink.samples(ixEyePos,:);
            XY = (cm(:,:,eyeIdx)*[ptbElXY; ones(1,size(ptbElXY,2))]);
        else 
            ixEyePos=(cellfun(@(x) any(strfind(x, 'EyeX')), sampleNames)) | (cellfun(@(x) any(strfind(x, 'EyeY')), sampleNames)); % index to find Eye position
            XY = trial(iTrial).eyelink.samples(ixEyePos,:);
        end
        ixPupil=cellfun(@(x) any(strfind(x, 'Pupil')), sampleNames);

        ptbElTime = trial(iTrial).eyelink.samples(1,:);
       
        % time in PTB clock
        time = ptbElTime/1e3-trial(iTrial).timing.eyelinkStartTime(2) + trial(iTrial).timing.eyelinkStartTime(1);
        % time in OE clock
        time = PDS{1}.PTB2OE(time);
        
        
        
        pupil = trial(iTrial).eyelink.samples(ixPupil,:);
        
        % if track is lost
        bad = pupil==0;
        
        % convert to degrees
        XY = pds.px2deg(bsxfun(@minus, XY, ctr(:)), trial(1).display.viewdist, trial(1).display.px2w);
        
        XY(:,bad) = nan;
        pupil(bad) = nan;
        
        if diff(time) > (1/500)
            eyePos = session.eyepos(time(:), XY(1,:)', XY(2,:)', pupil', pupil');
            ep = eyePos.resample('fs', 1e3);
            time = ep.tsample(:)';
            XY = [ep.x ep.y]';
            pupil = (ep.pwdth.^2*pi)';
        end
        timestamps = [timestamps time(:)'];
        data = [data [XY; pupil]];
    end
    
%     [timestamps, id] = sort(timestamps);
%     data = data(:,id);
%     
%     bad = timestamps < 0;
%     timestamps(bad) = [];
%     data(:,bad) = [];
end
    
    