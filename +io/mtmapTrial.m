function mtmapTrial = mtmapTrial(PDS, varargin)

stim = 'DotMapping';

ppd  = PDS{1}.initialParametersMerged.display.ppd;
ifi  = PDS{1}.initialParametersMerged.display.ifi;

hasStim = io.findPDScontainingStimModule(PDS, stim);
hasFixation = cellfun(@(x) strncmp(x.initialParametersMerged.pldaps.trialFunction, 'stimuli.fixflash', 16), PDS);

hasStim = hasStim(:); % & ~hasFixation(:);

mtmapTrial = struct();
trialNum = 0;

for i = find(hasStim(:)')
    
    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 
    
    stimTrials = find(trialIx);
    
    if isempty(stimTrials)
        continue
    end
    
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
        try
            kTrial = trialNum + j;
            
            % --- Timing
            mtmapTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
            mtmapTrial(kTrial).start      = mtmapTrial(kTrial).frameTimes(1);
            mtmapTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - mtmapTrial(kTrial).start;
            
            % --- Eye position
            eyepos = io.getEyePosition(PDS{i}, thisTrial);
            mtmapTrial(kTrial).eyeSampleTime = eyepos(:,1);
            mtmapTrial(kTrial).eyeXPx        = eyepos(:,2);
            mtmapTrial(kTrial).eyeYPx        = eyepos(:,3);
            mtmapTrial(kTrial).pupilArea     = eyepos(:,4);
            
            results = pdsa.detectSaccades(eyepos(:,1)', eyepos(:,2:3)'./ppd, 'verbose', false);
            
            
            
            iix = mtmapTrial(kTrial).eyeXPx < 200 | mtmapTrial(kTrial).eyeXPx > 1800;
            iiy = mtmapTrial(kTrial).eyeYPx < 100 | mtmapTrial(kTrial).eyeXPx > 1000;
            bad = iix | iiy;
            
            mtmapTrial(kTrial).eyeXPx(bad) = nan;
            mtmapTrial(kTrial).eyeYPx(bad) = nan;
            mtmapTrial(kTrial).pupilArea(bad) = nan;
            
            % find eye position on each frame
            t = mtmapTrial(kTrial).frameTimes - mtmapTrial(kTrial).start;
            tdiff = abs(bsxfun(@minus, mtmapTrial(kTrial).eyeSampleTime, t(:)')) < mean(diff(mtmapTrial(kTrial).eyeSampleTime));
            [irow,icol] = find(diff(tdiff)==-1);
            
            mtmapTrial(kTrial).saccades = false(size(mtmapTrial(kTrial).frameTimes));
            mtmapTrial(kTrial).saccadeTimes = results(1,:)+mtmapTrial(kTrial).start;
            
            sacDx = results(7,:) - results(5,:); % dx
            sacDy = results(8,:) - results(6,:); % dy
            [th, rho] = cart2pol(sacDx, sacDy);
            mtmapTrial(kTrial).saccadeDirection = th/pi*180;
            mtmapTrial(kTrial).saccadeAmp       = rho;
            
            for iSaccade = 1:size(results,2)
                ix = t > (results(1, iSaccade)-0.05) & t < (results(2,iSaccade)+0.05);
                mtmapTrial(kTrial).saccades(ix) = true;
            end
            
            mtmapTrial(kTrial).eyePosAtFrame = [mtmapTrial(kTrial).eyeXPx(irow) mtmapTrial(kTrial).eyeYPx(irow)];
            mtmapTrial(kTrial).pupilAtFrame = mtmapTrial(kTrial).pupilArea(irow);
            
            nApert = numel(mtmapTrial(kTrial).frameTimes);
            mtmapTrial(kTrial).xpos       = PDS{i}.data{thisTrial}.(stim).x(1:nApert,:);
            mtmapTrial(kTrial).ypos       = PDS{i}.data{thisTrial}.(stim).y(1:nApert,:);
            mtmapTrial(kTrial).on         = PDS{i}.data{thisTrial}.(stim).on(1:nApert,:);
            
            mtmapTrial(kTrial).xpos(~mtmapTrial(kTrial).on) = nan;
            mtmapTrial(kTrial).ypos(~mtmapTrial(kTrial).on) = nan;
            
            mtmapTrial(kTrial).direction  = PDS{i}.data{thisTrial}.(stim).direction(1:nApert,:);
            mtmapTrial(kTrial).speed      = PDS{i}.data{thisTrial}.(stim).speed(1:nApert,:) / ifi / ppd;
            
            mtmapTrial(kTrial).dx      = mtmapTrial(kTrial).speed .* cosd(mtmapTrial(kTrial).direction);
            mtmapTrial(kTrial).dy      = mtmapTrial(kTrial).speed .* sind(mtmapTrial(kTrial).direction);
                        
            mtmapTrial(kTrial).xposEye    = bsxfun(@minus, mtmapTrial(kTrial).xpos, mtmapTrial(kTrial).eyePosAtFrame(:,1));
            mtmapTrial(kTrial).yposEye    = -bsxfun(@minus, mtmapTrial(kTrial).ypos, mtmapTrial(kTrial).eyePosAtFrame(:,2));
            
            mtmapTrial(kTrial).apertureSize = [PDS{i}.data{thisTrial}.(stim).hDots.maxRadius]/ppd;
        end
        
    end
    
    trialNum = kTrial;
end