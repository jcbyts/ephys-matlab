classdef psaForage < handle
    % psaForage Pre Saccadic Attention Forage is a class for importing the
    % targetselection forage paradigm
   properties
       numTrials % number of trials run
       trial     % struct-array of trial values
       display   % display parameters from that session
   end
   
   methods
       
       function o = psaForage(PDS, varargin) % constructor
           
           if ~iscell(PDS)
               PDS = {PDS};
           end
           
           % --- find targetselection or dotselection trials
           stim = 'targetselection';
           
           hasStim = io.findPDScontainingStimModule(PDS, stim);
           
           stim = 'dotselection';
           
           hasStim = hasStim | io.findPDScontainingStimModule(PDS, stim);
           
           if ~any(hasStim)
               return
           end
           
           o.display = PDS{find(hasStim, 1, 'first')}.initialParametersMerged.display;
           
           
           for i = find(hasStim(:)')
               
               trial_ = o.importPDS(PDS{i});
               
               if isempty(trial_)
                   continue
               end
               
               o.trial = [o.trial; trial_(:)];
               
           end
           
           
           o.numTrials = numel(o.trial);
           
       end
       
   end
   
   methods (Static)
       function trial = importPDS(PDS)
           
           % importPDS checks which version of to stimulus code was run
            % and imports to a common format appropriately
            
            pdsDate = PDS.initialParametersMerged.session.initTime;
            if isfield(PDS.initialParametersMerged.git, 'pep')
                
                if any(strfind(PDS.initialParametersMerged.git.pep.status, 'branch cleanup'))
                    
                    if pdsDate > datenum(2018, 12, 12)
                        trial = session.psaForage.importPDS_v2(PDS);
                    elseif pdsDate > datenum(2018, 02, 01)
                        trial = session.psaForage.importPDS_v1(PDS);
                    else
                        error('unknown version')
                    end
                    
                else
                    error('unknown version')
                end
            else
                try
                     trial = session.psaForage.importPDS_v1(PDS);
                catch me
                    error('unknown version')
                end
            end
           
       end
       
       function psaTrial = importPDS_v2(PDS)
          error('not implemented yet')
       end
       
       function psaTrial = importPDS_v1(PDS)
           psaTrial = struct();
           
           % --- find CSD flash trials
           stim = 'dotselection';
           
           trialIx = cellfun(@(x) isfield(x, stim), PDS.data);
           
           stimTrials = find(trialIx);
           
           dot1RewardRate = PDS.initialParametersMerged.(stim).rewardDot1Rate;
           dot2RewardRate = PDS.initialParametersMerged.(stim).rewardDot2Rate;
           stimVisible    = PDS.initialParametersMerged.(stim).stimVisible;
           
           for j = 1:numel(stimTrials)
               thisTrial = stimTrials(j);
               
               kTrial = j;
               
               psaTrial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1));
               psaTrial(kTrial).start      = psaTrial(kTrial).frameTimes(1);
               psaTrial(kTrial).duration   = psaTrial(kTrial).frameTimes(end) - psaTrial(kTrial).start;
               
               fixBehavior = PDS.initialParametersMerged.(stim).fixationBehavior;
               
               % --- Use the stimulus object logs to recreate the timing
               
               % fixation point onset
               onsetIndex = PDS.data{thisTrial}.(fixBehavior).hFix.log(1,:)==1;
               psaTrial(kTrial).fixOn      = PDS.PTB2OE(PDS.data{thisTrial}.(fixBehavior).hFix.log(2,onsetIndex));
               
               % fixation point offset
               offsetIndex = PDS.data{thisTrial}.(fixBehavior).hFix.log(1,:)==0;
               psaTrial(kTrial).fixOff     = PDS.PTB2OE(PDS.data{thisTrial}.(fixBehavior).hFix.log(2,offsetIndex));
               
               % fixation entered
               ix = PDS.data{thisTrial}.(fixBehavior).hFix.fixlog(1,:)==1;
               psaTrial(kTrial).fixEntered = PDS.PTB2OE(PDS.data{thisTrial}.(fixBehavior).hFix.fixlog(2,ix));
               
               % final fixation point offset (transition to state 2)
               psaTrial(kTrial).goSignal   = PDS.PTB2OE(PDS.data{thisTrial}.(stim).states.getTxTime(2));
               
               % --- Targets
               psaTrial(kTrial).targets    = PDS.data{thisTrial}.(stim).hTargs;
               psaTrial(kTrial).numTargs   = numel(psaTrial(kTrial).targets);
               
               % --- Loop over targets and track when they came on, off,
               % position, etc.
               
               % preallocate some variables
               nt = max(arrayfun(@(x) size(x.log,2), psaTrial(kTrial).targets));
               psaTrial(kTrial).targsOn  = nan(psaTrial(kTrial).numTargs, nt);
               psaTrial(kTrial).targsOff = nan(psaTrial(kTrial).numTargs, nt);
               
               for iTarg = 1:psaTrial(kTrial).numTargs
                   
                    % --- Target position
                    
                    % in pixels
                    xpx = psaTrial(kTrial).targets(iTarg).position(1);
                    ypx = psaTrial(kTrial).targets(iTarg).position(2);
                    
                    % subtract off center of the screen
                    xpx = xpx - PDS.initialParametersMerged.display.ctr(1);
                    ypx = ypx - PDS.initialParametersMerged.display.ctr(2);
                    
                    % flip Y position so positive is up
                    ypx = -ypx;
                    
                    dxy = pds.px2deg([xpx; ypx], PDS.initialParametersMerged.display.viewdist, PDS.initialParametersMerged.display.px2w);
                    psaTrial(kTrial).targPosX(iTarg) = dxy(1);
                    psaTrial(kTrial).targPosY(iTarg) = dxy(2);
                   
                    % --- Target features
                    psaTrial(kTrial).targDirection(iTarg) = psaTrial(kTrial).targets(iTarg).theta;
                    
                    tmp = pds.px2deg(psaTrial(kTrial).targets(iTarg).radius,PDS.initialParametersMerged.display.viewdist, PDS.initialParametersMerged.display.px2w);
                    psaTrial(kTrial).targRadius(iTarg)   = tmp(1);
                    
                    if isa(psaTrial(kTrial).targets(iTarg), 'stimuli.objects.gaborTarget')
                        psaTrial(kTrial).targSpeed(iTarg) = psaTrial(kTrial).targets(iTarg).tf/psaTrial(kTrial).targets(iTarg).sf;
                    else
                        error('Need to implement this for dots')
                    end
                    
                    % --- Timings
                    
                    % initialize with nan
                    psaTrial(kTrial).targsOn(iTarg,1) = nan;
                    psaTrial(kTrial).targsOff(iTarg,1) = nan;
                    
                    
                   % --- Log Target Onset
                   if isempty(psaTrial(kTrial).targets(iTarg).log) % target never turned on
                       continue
                   end
                   
                   
                   ix = psaTrial(kTrial).targets(iTarg).log(1,:) == 1; 
                   if ~any(ix) % target never turned on
                       continue
                   end
                   
                   targon = PDS.PTB2OE(psaTrial(kTrial).targets(iTarg).log(2,ix));
                   nt = numel(targon);
                   if  nt > 1 % target turned on more than once?
                       % check if it was logged twice within a frame
                       if diff(targon) < PDS.initialParametersMerged.display.ifi
                           nt = 1;
                       end
                   end
                   psaTrial(kTrial).targsOn(iTarg, 1:nt) = targon(1:nt);
                   
                   % --- Log Target Offset
                   ix = psaTrial(kTrial).targets(iTarg).log(1,:) == 0;
                   if ~any(ix)
                       % target turned off at last frame
                        psaTrial(kTrial).targsOff(iTarg, 1) = PDS.PTB2OE(PDS.data{thisTrial}.(stim).states.getTxTime(PDS.data{thisTrial}.(stim).states.stateId));
                   else
                       targoff = PDS.PTB2OE(psaTrial(kTrial).targets(iTarg).log(2,ix));
                       if numel(targoff) > 1 % why did it turn off twice
                           assert(diff(targoff) < PDS.initialParametersMerged.display.ifi, 'Target turned off twice. What is up?')
                       end
                       
                       psaTrial(kTrial).targsOff(iTarg, 1:nt)= targoff(1:nt);
                   end
                   
                   
               end
               
               psaTrial(kTrial).targChosen = PDS.data{thisTrial}.(stim).dotsChosen;
               if isnan(psaTrial(kTrial).targChosen)
                   psaTrial(kTrial).choiceTime = nan;
               else
                   ix = psaTrial(kTrial).targets(psaTrial(kTrial).targChosen).fixlog(1,:)==1;
                   psaTrial(kTrial).choiceTime = PDS.PTB2OE(psaTrial(kTrial).targets(psaTrial(kTrial).targChosen).fixlog(2,ix));
               end
               
               
               
               % check if reward conditions changed
               fnames = fieldnames(PDS.conditions{thisTrial}.(stim));
               if any(strcmp(fnames, 'rewardDot1Rate'))
                    dot1RewardRate = PDS.conditions{thisTrial}.(stim).rewardDot1Rate;
               end
               
               if any(strcmp(fnames, 'rewardDot2Rate'))
                   dot2RewardRate = PDS.conditions{thisTrial}.(stim).rewardDot2Rate;
               end
               
               if any(strcmp(fnames, 'stimVisible'))
                   stimVisible = PDS.conditions{thisTrial}.(stim).stimVisible;
               end
               
               psaTrial(kTrial).rewardRate = [dot1RewardRate dot2RewardRate];
               psaTrial(kTrial).stimVisible = stimVisible;
                
               % check for change in reward rate
               if kTrial > 1
                   if ~all((psaTrial(kTrial).rewardRate - psaTrial(kTrial-1).rewardRate) == 0)
                       psaTrial(kTrial).switchReward = 1;
                   else
                       psaTrial(kTrial).switchReward = 0;
                   end
               else
                   psaTrial(kTrial).switchReward = 1;
               end
               
               psaTrial(kTrial).isRewarded = PDS.data{thisTrial}.(stim).isRewarded;
               
               psaTrial(kTrial).eyePosAtFrame = PDS.data{thisTrial}.behavior.eyeAtFrame';
                   
               % center
               psaTrial(kTrial).eyePosAtFrame = bsxfun(@minus, psaTrial(kTrial).eyePosAtFrame, PDS.initialParametersMerged.display.ctr(1:2));
               
               % flip Y
               psaTrial(kTrial).eyePosAtFrame(:,2) = -psaTrial(kTrial).eyePosAtFrame(:,2);
               
               % pixels 2 degrees
               psaTrial(kTrial).eyePosAtFrame = pds.px2deg(psaTrial(kTrial).eyePosAtFrame', PDS.initialParametersMerged.display.viewdist, PDS.initialParametersMerged.display.px2w)';
               
               % align stimuli that are yoked to the frame to the frame
               % rate
%                psaTrial(kTrial).fixOn = session.psaForage.alignToNextFrame(psaTrial(kTrial).fixOn, psaTrial(kTrial).frameTimes);
%                psaTrial(kTrial).fixOff = session.psaForage.alignToNextFrame(psaTrial(kTrial).fixOff, psaTrial(kTrial).frameTimes);
               
               
           end
           
           
       end
       
       function alignedTimes = alignToNextFrame(cpuTime, frameTimes)
            diffs = bsxfun(@minus, cpuTime(:)', frameTimes(:))<0;
            n = numel(cpuTime);
            alignedTimes = nan(size(cpuTime));
            for i = 1:n
                alignedTimes(i) = frameTimes(find(diffs(:,i), 1, 'first'));
            end
           
       end
   end
end