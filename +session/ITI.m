classdef ITI < handle
    % Square flash for spatial mapping (Full Field)
    
    properties
        numTrials
        trial
        display
        design
    end
    
    methods
        
        function h = ITI(PDS)
           
            if isstruct(PDS)
                PDS = {PDS};
            end
            
            
            % Loop over PDS files, importing the stimulus parameters
            for i = 1:numel(PDS)
                
                % helper function for the import
                trial_ = h.importPDS(PDS{i});
                
                if isempty(trial_)
                    continue
                end
                
                h.trial = [h.trial; trial_(:)];
                
            end
            
            tends = arrayfun(@(x) x.end, h.trial);
            tstarts = arrayfun(@(x) x.start, h.trial);
            
            h.numTrials = numel(h.trial)-1;
            h.trial = [];
            for i = 1:h.numTrials
                h.trial(i).start = tends(i);
                h.trial(i).stop   = tstarts(i+1);
                h.trial(i).duration = h.trial(i).stop - h.trial(i).start;
            end
            % assuming the display parameters didn't change during the
            % session, store them in the object
            h.display = PDS{1}.initialParametersMerged.display;
            
        end
        
    end
    
     methods (Static)
        
        function trial = importPDS(PDS)
            % importPDS checks which version of to stimulus code was run
            % and imports to a common format appropriately
            
            trial = [];
            numTrials = numel(PDS.data);
            
            for i = 1:numTrials
                trial(i).start = PDS.PTB2OE(PDS.data{i}.timing.flipTimes(1,1));
                trial(i).end   = PDS.PTB2OE(PDS.data{i}.timing.flipTimes(1,end));
            end
        end
        
     end
    
    
end