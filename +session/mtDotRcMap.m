classdef mtDotRcMap < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        numTrials
        trial
        display
        design
    end
    
    methods
        function m = mtDotRcMap(PDS, fixationOnly)
            
            if nargin < 2
                fixationOnly = false;
            end
            
            stim = 'DotMapping';
            
            ppd  = PDS{1}.initialParametersMerged.display.ppd;
            ifi  = PDS{1}.initialParametersMerged.display.ifi;
            
            
            hasStim = io.findPDScontainingStimModule(PDS, stim);
            if ~any(hasStim)
                return
            end
            
            m.display = PDS{find(hasStim,1)}.initialParametersMerged.display;
            
            hasFixaton = cellfun(@(x) strncmp(x.initialParametersMerged.pldaps.trialFunction, 'stimuli.fixflash', 16), PDS);
            
            if fixationOnly
                hasStim = hasStim(:) & hasFixaton(:);
            end
            
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
                        [irow,~] = find(diff(tdiff)==-1);
                        
                    
                        mtmapTrial(kTrial).saccades = false(size(mtmapTrial(kTrial).frameTimes));
                        mtmapTrial(kTrial).saccadeTimes = results(1,:)+mtmapTrial(kTrial).start;
                        
                        % 5. StartX
                        % 6. StartY
                        % 7. EndX
                        % 8. EndY
                        sacDx = results(7,:) - results(5,:);
                        sacDy = results(8,:) - results(6,:);
                        [th, rho] = cart2pol(sacDx, sacDy);
                        mtmapTrial(kTrial).saccadeDirection = th/pi*180;
                        mtmapTrial(kTrial).saccadeAmp       = rho;
                        
                        if fixationOnly
                            fixations = diff(results(1,:))>.35;
                            if any(fixations)
                                for iSaccade = find(fixations(:)')
                                    ix = t > results(2,iSaccade) + .05 && t < results(1,iSaccade+1);
                                    
                                    mtmapTrial(kTrial).saccades(ix) = true;
                                end
                            end
                        else
                            
                            for iSaccade = 1:size(results,2)
                                ix = t > (results(1, iSaccade)-0.05) & t < (results(2,iSaccade)+0.05);
                                mtmapTrial(kTrial).saccades(ix) = true;
                            end
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
            
            try
                mtmapTrial(arrayfun(@(x) isempty(x.xpos),mtmapTrial)) = [];
            end
            m.trial = mtmapTrial;
            m.numTrials = numel(mtmapTrial);
            
        end
        
        function plotTrial(m, kTrial)
            
            if nargin < 2
                kTrial = randi(m.numTrials);
            end
            
            ppd = m.display.ppd;
            subplot(1,3,1)
            plot(m.trial(kTrial).xpos, m.trial(kTrial).ypos, '.k')
            hold on
            plot(m.trial(kTrial).eyePosAtFrame(:,1), m.trial(kTrial).eyePosAtFrame(:,2), '.r')
            
            subplot(1,3,2)
            plot(m.trial(kTrial).eyePosAtFrame(:,1), 'r'); hold on
            plot(m.trial(kTrial).xpos, 'k.')
            
            subplot(1,3,3)
            plot(m.trial(kTrial).xposEye/ppd, m.trial(kTrial).yposEye/ppd, '.k')
            
            quiver(m.trial(kTrial).xposEye/ppd, m.trial(kTrial).yposEye/ppd, m.trial(kTrial).dx/5, m.trial(kTrial).dy/5, 'AutoScale', 'off')
            
            
        end
        
        
        
        function buildDesignMatrix(m, varargin)
            
            % --- check for optional arguments
            ip = inputParser();
            ip.addParameter('trialIdx', 1:m.numTrials)
            ip.addParameter('nTimeLags', 30)
            ip.addParameter('xwin', [0 8])
            ip.addParameter('ywin', [-8 2])
            ip.addParameter('binSize', 1)
            ip.addParameter('fullSpatioTemporal', false)
            
            ip.parse(varargin{:})
            
            fullSpatioTemporal = ip.Results.fullSpatioTemporal;
            
            xwin  = ip.Results.xwin;
            ywin  = ip.Results.ywin;
            
            ppd = m.display.ppd;
            
            binSize = ip.Results.binSize;
            xax = xwin(1):binSize:xwin(2);
            yax = ywin(1):binSize:ywin(2);
            
            [xgrid,ygrid]=meshgrid(xax, yax);
            
            nTotalFrames = sum(arrayfun(@(x) numel(x.frameTimes), m.trial));
            
            sz = size(xgrid);
            
            dx     = zeros(nTotalFrames, prod(sz));
            dy     = zeros(nTotalFrames, prod(sz));
            spcnt  = zeros(nTotalFrames, 1);
            
            
            
            iFrame = 0;
            
            nTrials = numel(m.trial);
            for kTrial = 1:nTrials
                fprintf('%d / %d \n', kTrial, nTrials)
                frameIdx = iFrame + (1:numel(m.trial(kTrial).frameTimes));
                
                ixGood = ~m.trial(kTrial).saccades;
                
                if any(ixGood)
                    
                    frameIdx = frameIdx(ixGood);
                    
                    x = m.trial(kTrial).xposEye(ixGood,:)/ppd;
                    y = m.trial(kTrial).yposEye(ixGood,:)/ppd;
                    
                    nGrid  = prod(sz);
                    nFrame = size(x,1);
                    nApert = size(x,2);
                    s = unique(m.trial(kTrial).apertureSize);
                    
                    for iDot = 1:nApert
                        xd = x(:, iDot) - xgrid(:)';
                        yd = y(:, iDot) - ygrid(:)';
                        r = sqrt(xd.^2 + yd.^2);
                        r(isnan(r)) = inf;
                        idx = r <= s | r <= binSize;
                        [i, j] = find(idx);
                        %         ind = find(idx);
                        
                        ind = sub2ind([nFrame nGrid], i, j);
                        dxtmp = zeros(nFrame, nGrid);
                        dytmp = zeros(nFrame, nGrid);
                        dxtmp(ind) = m.trial(kTrial).dx(i,iDot);
                        dytmp(ind) = m.trial(kTrial).dy(i,iDot);
                        
                        dx(frameIdx,:) = dx(frameIdx,:) + dxtmp;
                        dy(frameIdx,:) = dy(frameIdx,:) + dytmp;
                        
                    end
                end
                
%                 spcnt(frameIdx) = histc(sp.st, m.trial(kTrial).frameTimes);
                
                iFrame = iFrame + nFrame;
                
            end
            
            %% set up for regression
            ifi = m.display.ifi;
            nkt = ceil(.15/ifi);
            m.design.X = [dx dy];
            
            flipTimes = [m.trial.frameTimes];
            m.design.rowTimes = flipTimes(:);
            
            if fullSpatioTemporal
                X = double(dx > 0 |  dy > 0);
                Xd = rfmap.makeStimRowsDense(X , nkt);
                
                m.design.biasCol = size(Xd,2)+1;
                m.design.Xd = [Xd ones(nTotalFrames,1)];
                
                m.design.nkt = nkt;
                
                m.design.XX = m.design.Xd'*m.design.Xd;
                m.design.nkTime = nkt;
            end
            m.design.xax = xax;
            m.design.yax = yax;
            m.design.sz = [numel(xax) numel(yax)];
            
            
        end
        
    end
end


%
%
% %% plot sample trial
% kTrial = 1;
% figure(1); clf;
%
%
% %%
%
% % build grid
%
% %%
% saccades  = cell2mat(arrayfun(@(x) x.saccades(:), mtmapTrial, 'UniformOutput', false)');
%
% % C = (Xd(~saccades,:)'*Xd(~saccades,:));
%
% %%
% saccades = false(size(saccades));
% sta = Xd(~saccades,:)'*spcnt(~saccades);
%
% % sta = (C + 10e7*speye(nkt*nGrid*2+1))\sta;
%
% figure(2); clf
% plot(sta(1:end-1))
%
% I = reshape(sta(1:end-1), nkt, []);
% subplot(1,2,1)
% imagesc(I)
%
% subplot(1,2,2)
% imagesc(reshape(mean(Xd(~saccades,1:end-1)), nkt, []));
%
% figure(1); clf
% [u,s,v] = svd(I);
% timeK  = u(:,1);
% spaceK = v(:,1);
% signFlip = sign(sum(timeK));
%
% timeK = timeK * signFlip;
% spaceK = spaceK * signFlip;
% subplot(1,2,1)
% quiver(xgrid(:), ygrid(:), spaceK(1:nGrid), spaceK(nGrid+1:end))
% subplot(1,2,2)
% plot( (1:nkt) * ifi, timeK)