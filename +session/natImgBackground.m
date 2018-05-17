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
            % importPDS checks which version of to stimulus code was run
            % and imports to a common format appropriately
            
            pdsDate = PDS.initialParametersMerged.session.initTime;
            if isfield(PDS.initialParametersMerged.git, 'pep')
                
                if any(strfind(PDS.initialParametersMerged.git.pep.status, 'branch cleanup'))
                    
                    if pdsDate > datenum(2018, 02, 01)
                        trial = session.natImgBackground.importPDS_v2(PDS);
                    else
                        error('unknown version')
                    end
                    
                else
                    error('unknown version')
                end
            else
                trial = session.natImgBackground.importPDS_v1(PDS);
            end
        end
        
        function trial = importPDS_v2(PDS)
            trial = struct();
            
            % --- find CSD flash trials
            stim = 'natImgBackground';
            
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
                return
            end
            
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
