classdef hartleyFF < handle
    % Hartley FF (Full Field)
    
    properties
        numTrials
        trial
        display
        design
    end
    
    methods
        function h = hartleyFF(PDS, varargin)
            
            stim = 'hartley';
            
            
            hasStim = io.findPDScontainingStimModule(PDS, stim);
            
            if ~any(hasStim)
                return
            end
            
            h.display = PDS{1}.initialParametersMerged.display;
            
            for i = find(hasStim(:)')
                
                trial_ = h.importPDS(PDS{i});
                
                if isempty(trial_)
                    continue
                end
                
                h.trial = [h.trial; trial_(:)];
                
            end
            
            
            h.numTrials = numel(h.trial);
            
            
        end % constructor
        
        function plotTrial(h, kTrial)
            
            if nargin < 2
                kTrial = randi(h.numTrials);
            end
            
            plot(h.trial(kTrial).eyeXPx); hold on
            plot(h.trial(kTrial).eyeYPx);
            
        end
        
        % --- abstract methods
        buildDesignMatrix(h, varargin)
        
        sta = spikeTriggeredAverage(h, spikeTimes)
        
        sta = AsdRf(h, spikeTimes)        
        
        
    end
    
    methods (Static)
        
        [trial, display] = importPDS(PDS)
        
        [trial, display] = importPDS_v2(PDS)
        
        [trial, display] = importPDS_v1(PDS)
        
    end
end

