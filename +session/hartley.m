classdef hartleyFF < handle
    
    properties
        numTrials
        trial
    end
    
    methods
        function h = hartleyFF(PDS, varargin)
            
            stim = 'hartley';
            
            
            hasStim = io.findPDScontainingStimModule(PDS, stim);
            
            h.trial = struct();
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
                    
                    h.trial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
                    h.trial(kTrial).start      = h.trial(kTrial).frameTimes(1);
                    h.trial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - h.trial(kTrial).start;
                    
                    if isfield(PDS{i}.conditions{thisTrial}, stim)
                        if isfield(PDS{i}.conditions{thisTrial}.(stim), 'setupRNG')
                            if strcmp(PDS{i}.conditions{thisTrial}.(stim).setupRNG, 'frozenSequence')
                                h.trial(kTrial).frozenSequence = true;
                                h.trial(kTrial).frozenSequenceLength = PDS{i}.conditions{thisTrial}.(stim).sequenceLength;
                            else
                                h.trial(kTrial).frozenSequence = false;
                                h.trial(kTrial).frozenSequenceLength = nan;
                            end
                            
                        else
                            h.trial(kTrial).frozenSequence = false;
                            h.trial(kTrial).frozenSequenceLength = nan;
                        end
                        
                    else
                        h.trial(kTrial).frozenSequence = false;
                        h.trial(kTrial).frozenSequenceLength = nan;
                    end
                    
                    h.trial(kTrial).kx         = PDS{i}.data{thisTrial}.(stim).kx;
                    h.trial(kTrial).ky         = PDS{i}.data{thisTrial}.(stim).ky;
                    h.trial(kTrial).on         = PDS{i}.data{thisTrial}.(stim).on;
                    
                    eyepos = io.getEyePosition(PDS{i}, thisTrial);
                    h.trial(kTrial).eyeSampleTime = eyepos(:,1);
                    h.trial(kTrial).eyeXPx        = eyepos(:,2);
                    h.trial(kTrial).eyeYPx        = eyepos(:,3);
                    h.trial(kTrial).pupilArea     = eyepos(:,4);
                end
                
                trialNum = kTrial;
                
            end
            
            h.numTrials = numel(h.trial);
            
            
        end % constructor
        
        
        function plotTrial(h, kTrial)
            
            if nargin < 2
                kTrial = randi(h.numTrials);
            end
            
        end
        
    end
end

