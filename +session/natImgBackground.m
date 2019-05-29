classdef natImgBackground < handle
    % natImgBackground import module for natural images protocol
    %   Detailed explanation goes here
    
    properties
        numTrials
        trial
        display
        sessionTrialIdx=[]
    end
    
    methods
        
        function obj = natImgBackground(PDS)
            
            % --- find CSD flash trials
            stim = 'natImgBackground';
            
            [hasStim, numTrialsPerPDS] = io.findPDScontainingStimModule(PDS, stim);
            trialOffset = [0; cumsum(numTrialsPerPDS)];
            
            if ~any(hasStim)
                return
            end
            
            obj.display = PDS{find(hasStim, 1, 'first')}.initialParametersMerged.display;
            
            
            for i = find(hasStim(:)')
                
                [trial_, ~, idx] = obj.importPDS(PDS{i});
                obj.sessionTrialIdx = [obj.sessionTrialIdx trialOffset(i) + idx(:)'];
                
                if isempty(trial_) || isempty(fieldnames(trial_))
                    continue
                end
                
                obj.trial = [obj.trial; trial_(:)];
                
            end
            
            
            obj.numTrials = numel(obj.trial);
            
        end
    end
    
    methods (Static)
        function [trial, display, trialIdx] = importPDS(PDS)
            % importPDS checks which version of to stimulus code was run
            % and imports to a common format appropriately
            
            pdsDate = PDS.initialParametersMerged.session.initTime;
            if isfield(PDS.initialParametersMerged.git, 'pep')
                
                if any(strfind(PDS.initialParametersMerged.git.pep.status, 'branch cleanup'))
                    
                    if pdsDate > datenum(2018, 02, 01)
                        [trial, display, trialIdx] = session.natImgBackground.importPDS_v2(PDS);
                    else
                        error('unknown version')
                    end
                    
                else
                   warning('natImgBackground: git tracking failed. using version 2 import')
                   try
                       [trial, display, trialIdx] = session.natImgBackground.importPDS_v2(PDS);
                   catch
                       error('version 2 import failed')
                   end
                end
            else
                [trial, display, trialIdx] = session.natImgBackground.importPDS_v1(PDS);
            end
        end
        
        function [trial, display, stimTrials] = importPDS_v2(PDS)
            trial = struct();
            
            % --- find CSD flash trials
            stim = 'natImgBackground';
            display = [];
            pdstrial = pds.getPdsTrialData(PDS);
            
            trialIx = arrayfun(@(x) x.(stim).use, pdstrial);
            
            stimTrials = find(trialIx);
            
            
            if isempty(stimTrials)
                return
            end
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1));
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).stop       = trial(kTrial).frameTimes(end);
                trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - trial(kTrial).start;
                trial(kTrial).imgIdx     = PDS.data{thisTrial}.(stim).imgIndex(PDS.data{thisTrial}.(stim).texShown(1:end-1));
                
                fl = pdstrial(thisTrial).(stim).fileList;
                trial(kTrial).fileList = fl(unique(trial(kTrial).imgIdx));
                for k = 1:numel(trial(kTrial).fileList)
                    trial(kTrial).fileList(k).folder = pdstrial(thisTrial).(stim).imgDir;
                end
                    
            end
            
            % TODO: actually link to images
            %             imgsShown = unique(cell2mat(arrayfun(@(x) x.imgIdx(:), trial, 'UniformOutput', false)'));
            %             stim = 'natImgBackground';
            %             arrayfunPDS.initialParametersMerged.(stim).fileList(imgsShown)
            
            
        end
           
        
        function [trial, display, stimTrials] = importPDS_v1(PDS)
            trial = struct();
            display = [];
            % --- find CSD flash trials
            stim = 'natImgBackground';
            
            trialIx = cellfun(@(x) isfield(x, stim), PDS.data);
            
            stimTrials = find(trialIx);
            d = pds.getPdsTrialData(PDS);
            
            imgDirs = unique(arrayfun(@(x) x.natImgBackground.imgDir, d(:), 'uni', 0));
            assert(numel(imgDirs)==1, 'I only implemented handling one img directory per PDS file')
            % look for known image directories
            
            
            knownImageDirs = dir('Z:\Data\ImageDatabank');
            knownImageDirs(1:2) = [];
            knownImageDirs = knownImageDirs([knownImageDirs.isdir]);
            
            id = (arrayfun(@(x) any(strfind(imgDirs{1}, x.name)), knownImageDirs));
            imgDir = knownImageDirs(id);
            

            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1));
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).stop      = trial(kTrial).frameTimes(end);
                trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - trial(kTrial).start;
                imgIdx     = PDS.data{thisTrial}.(stim).imgIndex(PDS.data{thisTrial}.(stim).texShown(1:end-1));
                
                uniqueImages = unique(imgIdx);
                imagesShown = d(thisTrial).(stim).fileList(uniqueImages);
                nImagesShown = numel(imagesShown);
                trial(kTrial).images = cell(nImagesShown,1);
                for i = 1:nImagesShown
                    % TODO: this should be optional -- reading this image
                    % in is huge
%                     trial(kTrial).images{i} = imread(fullfile(imgDir.folder, imgDir.name, imagesShown(i).name));
                    trial(kTrial).imgIdx = imgIdx==uniqueImages(i);
                end
                    
            end
            
            % TODO: actually link to images
%             imgsShown = unique(cell2mat(arrayfun(@(x) x.imgIdx(:), trial, 'UniformOutput', false)'));
%             stim = 'natImgBackground';
%             arrayfunPDS.initialParametersMerged.(stim).fileList(imgsShown)
            
            
        end
    end
end
