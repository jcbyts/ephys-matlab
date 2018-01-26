classdef natImgBackground < handle
    % natImgBackground import module for natural images protocol
    %   Detailed explanation goes here
    
    properties
        numTrials
        trial
        display
    end
    
    methods
        
        function obj = natImgBackground(PDS)
            
            % --- find CSD flash trials
            stim = 'natImgBackground';
            
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
                trial = session.natImgBackground.importPDS_v1(PDS);
            end
                
        end
           
        
        function trial = importPDS_v1(PDS)
            trial = struct();
            
            % --- find CSD flash trials
            stim = 'natImgBackground';
            
            trialIx = cellfun(@(x) isfield(x, stim), PDS.data);
            
            stimTrials = find(trialIx);
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1));
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - trial(kTrial).start;
                trial(kTrial).imgIdx     = PDS.data{thisTrial}.(stim).imgIndex(PDS.data{thisTrial}.(stim).texShown(1:end-1));

            end
            
            % TODO: actually link to images
%             imgsShown = unique(cell2mat(arrayfun(@(x) x.imgIdx(:), trial, 'UniformOutput', false)'));
%             stim = 'natImgBackground';
%             arrayfunPDS.initialParametersMerged.(stim).fileList(imgsShown)
            
            
        end
    end
end
