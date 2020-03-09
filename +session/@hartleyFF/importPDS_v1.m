function [trial, display, stimTrials] = importPDS_v1(PDS)

trial   = [];
display = PDS.initialParametersMerged.display;


stim = 'hartley';

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
    
    if isfield(PDS.conditions{thisTrial}, stim)
        if isfield(PDS.conditions{thisTrial}.(stim), 'setupRNG')
            if strcmp(PDS.conditions{thisTrial}.(stim).setupRNG, 'frozenSequence')
                trial(kTrial).frozenSequence = true;
                trial(kTrial).frozenSequenceLength = PDS.conditions{thisTrial}.(stim).sequenceLength;
            else
                trial(kTrial).frozenSequence = false;
                trial(kTrial).frozenSequenceLength = nan;
            end
            
        else
            trial(kTrial).frozenSequence = false;
            trial(kTrial).frozenSequenceLength = nan;
        end
        
    else
        trial(kTrial).frozenSequence = false;
        trial(kTrial).frozenSequenceLength = nan;
    end
    
    trial(kTrial).kx         = PDS.data{thisTrial}.(stim).kx;
    trial(kTrial).ky         = PDS.data{thisTrial}.(stim).ky;
    trial(kTrial).on         = PDS.data{thisTrial}.(stim).on;
    trial(kTrial).phi        = PDS.data{thisTrial}.(stim).phi;
    trial(kTrial).tf         = PDS.data{thisTrial}.(stim).tf;
    
    eyepos = io.getEyePosition(PDS, thisTrial);
    trial(kTrial).eyeSampleTime = eyepos(:,1);
    trial(kTrial).eyeXPx        = eyepos(:,2);
    trial(kTrial).eyeYPx        = eyepos(:,3);
    trial(kTrial).pupilArea     = eyepos(:,4);
end

end