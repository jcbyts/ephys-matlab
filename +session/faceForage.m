classdef faceForage < handle
    % Square flash for spatial mapping (Full Field)
    
    properties
        numTrials
        trial
        display
        design
        sessionTrialIdx=[]
    end
    
    methods
        function h = faceForage(PDS, varargin)
            
            ip = inputParser();
            ip.addOptional('eyetrace', [])
            ip.addOptional('saccades', [])
            ip.parse(varargin{:})
            
            stim = 'forage';
            
            [hasStim, numTrialsPerPDS] = io.findPDScontainingStimModule(PDS, stim);
            trialOffset = [0; cumsum(numTrialsPerPDS)];
            
            hasStim = hasStim | io.findPDScontainingStimModule(PDS, 'faceforage');
            
            hasStim = hasStim | session.faceForage.checkForOriginalFaceForage(PDS);
            
            if ~any(hasStim)
                return
            end
            
            h.display = PDS{find(hasStim,1)}.initialParametersMerged.display;
                                    
            for i = find(hasStim(:)')
                
                [trial_, ~, idx] = h.importPDS(PDS{i});
                
                h.sessionTrialIdx = [h.sessionTrialIdx trialOffset(i) + idx(:)'];
                
                if isempty(trial_)
                    continue
                end
                
                h.trial = [h.trial; trial_(:)];
                
            end
            
            
            h.numTrials = numel(h.trial);
            
% %             % --- import eye position
% %             eyepos = ip.Results.eyetrace;
% %             if ~isempty(eyepos)
% %                if size(eyepos,2) < size(eyepos,1)
% %                    eyepos = eyepos';
% %                end
% %                
% %                assert(size(eyepos,1) >= 3, 'the first row must be timestamps')
% %                
% %                
% %                if ~isstruct(ip.Results.saccades)                   
% %                    sampleRate = 1/mode(diff(eyepos(1,:)));
% %                    ix = ~any(isnan(eyepos(2,:)));
% %                    [saccades] = pdsa.detectSaccades(eyepos(1,ix), eyepos(2:3,ix), ...
% %                        'verbose', false, ...
% %                        'filterPosition', 1, ...
% %                        'filterLength', ceil(20/sampleRate*1e3), ... % 40 ms smoothing for velocity computation
% %                        'detectThresh', 200, ...
% %                        'startThresh', 5, ...
% %                        'minIsi', ceil(50/sampleRate*1e3), ...
% %                        'minDur', ceil(4/sampleRate*1e3), ... % 4 ms
% %                        'blinkIsi', ceil(40/sampleRate*1e3));
% %                else
% %                    saccades = ip.Results.saccades;
% %                end
% %                
% %                if size(eyepos,1) == 4
% %                    pupil = eyepos(4,:);
% %                else
% %                    pupil = [];
% %                end
% %                
% %                for kTrial = 1:h.numTrials
% %                    
% %                    % include time before the trial starts
% %                    preTrial = 1.5;
% %                    
% %                    % find valid eye position
% %                    iix = (eyepos(1,:) > (h.trial(kTrial).start - preTrial)) & (eyepos(1,:) < (h.trial(kTrial).start + h.trial(kTrial).duration));
% %                    h.trial(kTrial).eyeSampleTime = eyepos(1,iix);
% %                    h.trial(kTrial).eyeXDeg = eyepos(2,iix);
% %                    h.trial(kTrial).eyeYDeg = eyepos(3,iix);
% %                    if ~isempty(pupil)
% %                         h.trial(kTrial).pupilArea = pupil(iix);
% %                    end
% %                    
% %                    % valid saccade times
% %                    iix = (saccades.start > (h.trial(kTrial).start - preTrial)) & (saccades.end < (h.trial(kTrial).start + h.trial(kTrial).duration));
% %                    fields = fieldnames(saccades);
% %                    for iField = 1:numel(fields)
% %                        newfield = ['sac_' fields{iField}];
% %                        h.trial(kTrial).(newfield) = saccades.(fields{iField})(iix);
% %                    end
% %                end
% % 
% %             end
            
            
        end % constructor
        
        
        function plotTrial(h, kTrial)
            
            if nargin < 2
                kTrial = randi(h.numTrials);
            end
            
        end
        
       
        
        
        
    end
    
    methods (Static)
        
        function [trial, display, trialIdx] = importPDS(PDS)
            % importPDS checks which version of to stimulus code was run
            % and imports to a common format appropriately
            
            pdsDate = PDS.initialParametersMerged.session.initTime;
            if isfield(PDS.initialParametersMerged.git, 'pep')
                
                if any(strfind(PDS.initialParametersMerged.git.pep.status, 'branch cleanup'))
                    
                    % this date checking is just a crude way I've been
                    % tracking versions. We probably should do something
                    % more sophisticated, but I haven't had time. We do
                    % save the exact commit and git diff, but it seems like
                    % too much work to reconstruct from that. I don't know.
                    % -- Jake
                    if pdsDate > datenum(2018, 02, 01)
                        [trial, display, trialIdx] = session.faceForage.importPDS_v2(PDS);
                    else
                        error('unknown version')
                    end
                    
                else
                    warning('faceForage: pds-stimuli git either on wrong branch, or it is not setup to track. running version 2 import')
                    try
                        [trial, display, trialIdx] = session.faceForage.importPDS_v2(PDS);
                    catch 
                        error('version 2 import failed')
                    end
                end
            else
                [trial, display, trialIdx] = session.faceForage.importPDS_v0(PDS);
            end
            
        end
        
        function [trial, display, stimTrials] = importPDS_v2(PDS)
            
            trial   = [];
            display = PDS.initialParametersMerged.display;
            
            stim = 'faceforage';
            
            trialIx = cellfun(@(x) isfield(x, stim), PDS.data);
            trialIx = trialIx | cellfun(@(x) isfield(x, 'forage'), PDS.data);
            stimTrials = find(trialIx);
            
            % --- check for conditions
            % if we were using pldaps modular features, it is likely that
            % some trials do not include the hartley stimulus. Those trials
            % would've been set by the condition field of pldaps. Check for
            % the use of conditions and then check if hartley was used.
            condstim = 'forage';
            if ~isempty(PDS.conditions)
                condIx = cellfun(@(x) isfield(x, condstim), PDS.conditions(stimTrials));
                
                if any(condIx) % conditions were used (ignore trials that hartley wasn't shown)
                    
                    notUsed = cellfun(@(x) isfield(x.(condstim), 'use')  && x.(condstim).use==0, PDS.conditions(stimTrials(condIx)));
                    
                    trialList = stimTrials(condIx);
                    excludeTrials = trialList(notUsed);
                    
                    stimTrials = setdiff(stimTrials, excludeTrials);
                end
            end
            
            if isempty(stimTrials)
                return;
            end
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end)); %#ok<*AGROW>
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).stop       = trial(kTrial).frameTimes(end);
                trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end)) - trial(kTrial).start;
                
                nFrames = numel(trial(kTrial).frameTimes);
                try
                    trial(kTrial).xpos = PDS.data{thisTrial}.(stim).x(1:nFrames,:);
                    trial(kTrial).ypos = PDS.data{thisTrial}.(stim).y(1:nFrames,:);
                    trial(kTrial).hold = PDS.data{thisTrial}.(stim).ctrHold(1:nFrames,:);
                catch
                    stim = 'forage';
                    trial(kTrial).xpos = PDS.data{thisTrial}.(stim).x(1:nFrames,:);
                    trial(kTrial).ypos = PDS.data{thisTrial}.(stim).y(1:nFrames,:);
                    trial(kTrial).hold = PDS.data{thisTrial}.(stim).ctrHold(1:nFrames,:);
                end
% %                 eyepos = io.getEyePosition(PDS, thisTrial);
% %                 
% %                 trial(kTrial).eyeSampleTime = eyepos(:,1);
% %                 trial(kTrial).eyeXPx        = eyepos(:,2);
% %                 trial(kTrial).eyeYPx        = eyepos(:,3);
% %                 trial(kTrial).pupilArea     = eyepos(:,4);
% %                 
% %                 % --- additional eye position analyses
% %                 %                     results = pdsa.detectSaccades(eyepos(:,1)', eyepos(:,2:3)'./ppd, 'verbose', false);
% %                 
% %                 
% %                 
% %                 iix = trial(kTrial).eyeXPx < 200 | trial(kTrial).eyeXPx > 1800;
% %                 iiy = trial(kTrial).eyeYPx < 100 | trial(kTrial).eyeXPx > 1000;
% %                 bad = iix | iiy;
% %                 
% %                 trial(kTrial).eyeXPx(bad) = nan;
% %                 trial(kTrial).eyeYPx(bad) = nan;
% %                 trial(kTrial).pupilArea(bad) = nan;
% %                 
% % %                 % find eye position on each frame
% % %                 t = trial(kTrial).frameTimes - trial(kTrial).start;
% % %                 tdiff = abs(bsxfun(@minus, trial(kTrial).eyeSampleTime, t(:)')) < mean(diff(trial(kTrial).eyeSampleTime));
% % %                 [irow,~] = find(diff(tdiff)==-1);
% %                 
% %                 %                     trial(kTrial).saccades = false(size(trial(kTrial).frameTimes));
% %                 %
% %                 %                     for iSaccade = 1:size(results,2)
% %                 %                         ix = t > results(1, iSaccade) & t < results(2,iSaccade);
% %                 %                         trial(kTrial).saccades(ix) = true;
% %                 %                     end
% %                 
% % %                 trial(kTrial).eyePosAtFrame = PDS.data{thisTrial}.(stim).eyes(1:nFrames,:);
% %                 trial(kTrial).eyePosAtFrame = PDS.data{thisTrial}.behavior.eyeAtFrame';
% % %                 
% %                 iix = trial(kTrial).eyePosAtFrame(:,1) < 200 | trial(kTrial).eyePosAtFrame(:,1) > 1800;
% %                 iiy = trial(kTrial).eyePosAtFrame(:,2) < 100 | trial(kTrial).eyePosAtFrame(:,2) > 1000;
% %                 bad = iix | iiy;
% % %                 
% %                 trial(kTrial).eyePosAtFrame(bad,:) = nan;
% %                 trial(kTrial).eyePosAtFrame = bsxfun(@minus, trial(kTrial).eyePosAtFrame, display.ctr(1:2)) / display.ppd;
% %                 trial(kTrial).eyePosAtFrame(:,2) = -trial(kTrial).eyePosAtFrame(:,2);
            end
            
        end
        
        function [trial, display, stimTrials] = importPDS_v1(PDS)
            
            trial   = [];
            display = PDS.initialParametersMerged.display;
            
            stim = 'SpatialMapping';
            
            trialIx = cellfun(@(x) isfield(x, stim), PDS.data);
            
            stimTrials = find(trialIx);
            
            if isempty(stimTrials)
                return;
            end
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1)); %#ok<*AGROW>
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).stop       = trial(kTrial).frameTimes(end);
                trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - trial(kTrial).start;
                
                trial(kTrial).pos        = PDS.data{thisTrial}.(stim).pos;
                
                eyepos = io.getEyePosition(PDS, thisTrial);
                trial(kTrial).eyeSampleTime = eyepos(:,1);
                trial(kTrial).eyeXPx        = eyepos(:,2);
                trial(kTrial).eyeYPx        = eyepos(:,3);
                trial(kTrial).pupilArea     = eyepos(:,4);
                
                % --- additional eye position analyses
                %                     results = pdsa.detectSaccades(eyepos(:,1)', eyepos(:,2:3)'./ppd, 'verbose', false);
                
                
                
                iix = trial(kTrial).eyeXPx < 200 | trial(kTrial).eyeXPx > 1800;
                iiy = trial(kTrial).eyeYPx < 100 | trial(kTrial).eyeXPx > 1000;
                bad = iix | iiy;
                
                trial(kTrial).eyeXPx(bad) = nan;
                trial(kTrial).eyeYPx(bad) = nan;
                trial(kTrial).pupilArea(bad) = nan;
                
                % find eye position on each frame
                t = trial(kTrial).frameTimes - trial(kTrial).start;
                tdiff = abs(bsxfun(@minus, trial(kTrial).eyeSampleTime, t(:)')) < mean(diff(trial(kTrial).eyeSampleTime));
                [irow,~] = find(diff(tdiff)==-1);
                
                %                     trial(kTrial).saccades = false(size(trial(kTrial).frameTimes));
                %
                %                     for iSaccade = 1:size(results,2)
                %                         ix = t > results(1, iSaccade) & t < results(2,iSaccade);
                %                         trial(kTrial).saccades(ix) = true;
                %                     end
                
                trial(kTrial).eyePosAtFrame = [trial(kTrial).eyeXPx(irow) trial(kTrial).eyeYPx(irow)];
                
            end
            
        end
        
        
        function [trial, display, stimTrials] = importPDS_v0(PDS)
            trial = [];
            d = pds.getPdsTrialData(PDS);
            display = d(1).display;
            stimTrials = find(arrayfun(@(x) isfield(x.stimulus, 'motionType') & strcmp(x.stimulus.type, 'face'), d));
            
            stim = 'stimulus';
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(d(thisTrial).timing.flipTimes(1,1:end))'; %#ok<*AGROW>
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).stop       = trial(kTrial).frameTimes(end);
                trial(kTrial).duration   = PDS.PTB2OE(d(thisTrial).timing.flipTimes(1,end)) - trial(kTrial).start;
                
                nFrames = min(size(d(thisTrial).(stim).x, 1), numel(trial(kTrial).frameTimes));
                trial(kTrial).frameTimes = trial(kTrial).frameTimes(1:nFrames);
                
                trial(kTrial).xpos = d(thisTrial).(stim).x(1:nFrames,:);
                trial(kTrial).ypos = d(thisTrial).(stim).y(1:nFrames,:);
                trial(kTrial).hold = d(thisTrial).(stim).ctrHold(1:nFrames,:);
                trial(kTrial).motionType = d(thisTrial).(stim).motionType;
                trial(kTrial).objectType = 'face';
                
                trial(kTrial).stop       = trial(kTrial).frameTimes(end);
            end
            
        end
        
        function hasStim = checkForOriginalFaceForage(PDS)
            
            nPds = numel(PDS);
            hasStim = false(nPds,1);
            for i = 1:nPds
                d = pds.getPdsTrialData(PDS{i});
                
                hasStim(i) = any(arrayfun(@(x) isfield(x.stimulus, 'motionType'), d));
            end
        end
        
    end
end

