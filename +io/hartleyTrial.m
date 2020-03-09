classdef hartleyTrial < handle
    
    

function hartleyTrial = hartleyTrial(PDS, varargin)

stim = 'hartley';

hasStim = io.findPDScontainingStimModule(PDS, stim);

hartleyTrial = struct();
trialNum = 0;

for i = find(hasStim(:)')
    
    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 
    
    stimTrials = find(trialIx);
    
    if isempty(stimTrials)
        continue
    end
        
    kxs=PDS{i}.data{stimTrials(1)}.hartley.kxs;
    kys=PDS{i}.data{stimTrials(1)}.hartley.kys;
    
   
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
    
        kTrial = trialNum + j;
        
        hartleyTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
        hartleyTrial(kTrial).start      = hartleyTrial(kTrial).frameTimes(1);
        hartleyTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - hartleyTrial(kTrial).start;
    
        if isfield(PDS{i}.conditions{thisTrial}, stim)
            if isfield(PDS{i}.conditions{thisTrial}.(stim), 'setupRNG')
                if strcmp(PDS{i}.conditions{thisTrial}.(stim).setupRNG, 'frozenSequence')
                    hartleyTrial(kTrial).frozenSequence = true;
                    hartleyTrial(kTrial).frozenSequenceLength = PDS{i}.conditions{thisTrial}.(stim).sequenceLength;
                else
                    hartleyTrial(kTrial).frozenSequence = false;
                    hartleyTrial(kTrial).frozenSequenceLength = nan;
                end
            
            else
                hartleyTrial(kTrial).frozenSequence = false;
                hartleyTrial(kTrial).frozenSequenceLength = nan;
            end
            
        else
                hartleyTrial(kTrial).frozenSequence = false;
                hartleyTrial(kTrial).frozenSequenceLength = nan;
        end
        
        hartleyTrial(kTrial).kx         = PDS{i}.data{thisTrial}.(stim).kx;
        hartleyTrial(kTrial).ky         = PDS{i}.data{thisTrial}.(stim).ky;
        hartleyTrial(kTrial).on         = PDS{i}.data{thisTrial}.(stim).on;
        
        eyepos = io.getEyePosition(PDS{i}, thisTrial);
        hartleyTrial(kTrial).eyeSampleTime = eyepos(:,1);
        hartleyTrial(kTrial).eyeXPx        = eyepos(:,2);
        hartleyTrial(kTrial).eyeYPx        = eyepos(:,3);
        hartleyTrial(kTrial).pupilArea     = eyepos(:,4);
    end
    
    trialNum = kTrial;
    
end