classdef csdFlash < handle
    %CSDFLASH import module for csd flash protocol
    %   Detailed explanation goes here
    
    properties
        numTrials
        trial
        display
    end
    
    methods
        
        function obj = csdFlash(PDS)
            
            % --- find CSD flash trials
            stim = 'csdFlash';
            
            hasStim = io.findPDScontainingStimModule(PDS, stim);
            
            if ~any(hasStim)
                return
            end
            
            obj.display = PDS{find(hasStim, 1, 'first')}.initialParametersMerged.display;
            
            
            for i = find(hasStim(:)')
                
                trial_ = obj.importPDS(PDS{i});
                
                if isempty(trial_)
                    continue
                end
                
                obj.trial = [obj.trial; trial_(:)];
                
            end
            
            
            obj.numTrials = numel(obj.trial);
            
        end
    end
    
    methods (Static)
        function trial = importPDS(PDS)
            
            pdsDate = PDS.initialParametersMerged.session.initTime;
            
            if pdsDate > datenum(2018,02,01)
                trial = [];
            else
                trial = session.csdFlash.importPDS_v1(PDS);
            end
                
        end
           
        
        function csdTrial = importPDS_v1(PDS)
            csdTrial = struct();
            
            % --- find CSD flash trials
            stim = 'csdFlash';
            
            trialIx = cellfun(@(x) isfield(x, stim), PDS.data);
            
            stimTrials = find(trialIx);
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                csdTrial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1));
                csdTrial(kTrial).start      = csdTrial(kTrial).frameTimes(1);
                csdTrial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - csdTrial(kTrial).start;
                
                csdTrial(kTrial).on         = PDS.data{thisTrial}.(stim).on;
                csdTrial(kTrial).onset      = csdTrial(kTrial).frameTimes(diff(csdTrial(kTrial).on)==1);
            end
            
            
        end
    end
end
