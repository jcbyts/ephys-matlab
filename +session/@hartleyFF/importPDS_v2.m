function [trial, display, stimTrials] = importPDS_v2(PDS)
trial = [];
display = PDS.initialParametersMerged.display;


stim = 'hartley';

trialIx = cellfun(@(x) isfield(x, stim), PDS.data);

stimTrials = find(trialIx);

% --- check for conditions
% if we were using pldaps modular features, it is likely that
% some trials do not include the hartley stimulus. Those trials
% would've been set by the condition field of pldaps. Check for
% the use of conditions and then check if hartley was used.
if ~isempty(PDS.conditions)
    condIx = cellfun(@(x) isfield(x, stim), PDS.conditions(stimTrials));
    
    if any(condIx) % conditions were used (ignore trials that hartley wasn't shown)
        
        notUsed = cellfun(@(x) x.(stim).use==0, PDS.conditions(stimTrials(condIx)));
        
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
    
    trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1)); %#ok<*AGROW>
    trial(kTrial).start      = trial(kTrial).frameTimes(1);
    trial(kTrial).stop       = trial(kTrial).frameTimes(end);
    trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - trial(kTrial).start;
    
    trial(kTrial).kx         = PDS.data{thisTrial}.(stim).kx;
    trial(kTrial).ky         = PDS.data{thisTrial}.(stim).ky;
    trial(kTrial).on         = PDS.data{thisTrial}.(stim).on;
    trial(kTrial).phi        = PDS.data{thisTrial}.(stim).phi;
    trial(kTrial).tf         = PDS.data{thisTrial}.(stim).tf;
    
    %                 eyepos = io.getEyePosition(PDS, thisTrial);
    %                 trial(kTrial).eyeSampleTime = eyepos(:,1);
    %                 trial(kTrial).eyeXPx        = eyepos(:,2);
    %                 trial(kTrial).eyeYPx        = eyepos(:,3);
    %                 trial(kTrial).pupilArea     = eyepos(:,4);
    
    % --- need to add frozen sequence
    trial(kTrial).frozenSequence = false;
    trial(kTrial).frozenSequenceLength = nan;
    
    if isfield(PDS.conditions{kTrial}, stim) && isfield(PDS.conditions{kTrial}.(stim), 'generativeModel')
        if strcmp(PDS.conditions{kTrial}.(stim).generativeModel, 'frozen')
            trial(kTrial).frozenSequence = true;
            trial(kTrial).frozenSequenceLength = PDS.data{kTrial}.(stim).sequenceFrame/4;
            
            inds = reshape(1:PDS.data{kTrial}.(stim).sequenceFrame, trial(kTrial).frozenSequenceLength, 4);
            tmp_ = PDS.data{kTrial}.(stim).sequence.kx(inds);
            assert(all(all(bsxfun(@eq, tmp_, tmp_(:,1)))), 'something is wrong with the frozen trial')
        end
    end
    
    % --- add face targets
    if isfield(PDS.data{thisTrial}, 'faceforage')
        trial(kTrial).faceforageX = PDS.data{thisTrial}.faceforage.x;
        trial(kTrial).faceforageY = PDS.data{thisTrial}.faceforage.y;
        trial(kTrial).faceforageCtr = PDS.data{thisTrial}.faceforage.ctrHold;
    else
        trial(kTrial).faceforageX = nan;
        trial(kTrial).faceforageY = nan;
        trial(kTrial).faceforageCtr = nan;
    end
    
    
end


end