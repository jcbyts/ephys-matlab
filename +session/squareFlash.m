classdef squareFlash < handle
% Square flash for spatial mapping (Full Field)

    properties
        numTrials
        trial
        display
        design
    end
    
    methods
        function h = squareFlash(PDS, varargin)
            
            stim = 'SpatialMapping';
            
            
            hasStim = io.findPDScontainingStimModule(PDS, stim);
            
            h.display = PDS{1}.initialParametersMerged.display;
            
            h.trial = struct();
            trialNum = 0;
            
            ppd = PDS{1}.initialParametersMerged.display.ppd;
            
            for i = find(hasStim(:)')
                
                trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data);
                
                stimTrials = find(trialIx);
                
                if isempty(stimTrials)
                    continue
                end
                
              
                
                for j = 1:numel(stimTrials)
                    thisTrial = stimTrials(j);
                    
                    kTrial = trialNum + j;
                    
                    h.trial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
                    h.trial(kTrial).start      = h.trial(kTrial).frameTimes(1);
                    h.trial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - h.trial(kTrial).start;
                    
                    h.trial(kTrial).pos        = PDS{i}.data{thisTrial}.(stim).pos;
                    
                    eyepos = io.getEyePosition(PDS{i}, thisTrial);
                    h.trial(kTrial).eyeSampleTime = eyepos(:,1);
                    h.trial(kTrial).eyeXPx        = eyepos(:,2);
                    h.trial(kTrial).eyeYPx        = eyepos(:,3);
                    h.trial(kTrial).pupilArea     = eyepos(:,4);
                    
                    % --- additional eye position analyses
%                     results = pdsa.detectSaccades(eyepos(:,1)', eyepos(:,2:3)'./ppd, 'verbose', false);
        
        
        
                    iix = h.trial(kTrial).eyeXPx < 200 | h.trial(kTrial).eyeXPx > 1800;
                    iiy = h.trial(kTrial).eyeYPx < 100 | h.trial(kTrial).eyeXPx > 1000;
                    bad = iix | iiy;
                    
                    h.trial(kTrial).eyeXPx(bad) = nan;
                    h.trial(kTrial).eyeYPx(bad) = nan;
                    h.trial(kTrial).pupilArea(bad) = nan;
                    
                    % find eye position on each frame
                    t = h.trial(kTrial).frameTimes - h.trial(kTrial).start;
                    tdiff = abs(bsxfun(@minus, h.trial(kTrial).eyeSampleTime, t(:)')) < mean(diff(h.trial(kTrial).eyeSampleTime));
                    [irow,~] = find(diff(tdiff)==-1);
                    
%                     h.trial(kTrial).saccades = false(size(h.trial(kTrial).frameTimes));
%                     
%                     for iSaccade = 1:size(results,2)
%                         ix = t > results(1, iSaccade) & t < results(2,iSaccade);
%                         h.trial(kTrial).saccades(ix) = true;
%                     end
                    
                    h.trial(kTrial).eyePosAtFrame = [h.trial(kTrial).eyeXPx(irow) h.trial(kTrial).eyeYPx(irow)];
                    
                end
                
                trialNum = kTrial;
                
            end
            
            h.numTrials = numel(h.trial);
            
            
        end % constructor
        
        
        function plotTrial(h, kTrial)
            
            if nargin < 2
                kTrial = randi(h.numTrials);
            end
            
            x = squeeze(mean(h.trial(kTrial).pos([1 3], :, :), 1));
            y = squeeze(mean(h.trial(kTrial).pos([2 4], :, :), 1));
            fr = repmat((1:size(x,2))', 1, size(x,1))';
            
            plot3(h.trial(kTrial).frameTimes(fr(:))-h.trial(kTrial).start, x(:), y(:), '.'); hold on
           
            plot3(h.trial(kTrial).eyeSampleTime-h.trial(kTrial).eyeSampleTime(1), h.trial(kTrial).eyeXPx, h.trial(kTrial).eyeYPx, 'Linewidth', 1)
            
            plot3(h.trial(kTrial).frameTimes-h.trial(kTrial).start, h.trial(kTrial).eyePosAtFrame(:,1), h.trial(kTrial).eyePosAtFrame(:,2), '.');
            
            
        end
        
        
        function S = binSpace(h, varargin)
            % binSpace bins the stimulus at a specified resolution
            % Inputs (as argument pairs)
            %   'trialIdx'        - range of trials (default: 1:numTrials)
            %   'correctEyePos'   - retinal or screen coords (default: true)
            %   'window'  [1 x 2] - bounds of the window to analyze (degrees)
            %   'binSize' [1 x 1] - size of each bin (degrees)
            % Output
            %   S@struct 
            %       X [time x bins] - stimulus binned in space. Each row is
            %                         the stimulus frame reshaped to be a
            %                         single vector of size diff(window)/binSize
            %       bins [time x 1] - the time of each row of X
            %       xax  [1 x bins] - x axis (in degrees)
            %       yax  [1 x bins] - y axis (in degrees)
            %       size [1 x 2]    - size of spatial stimulus
            %
            % Example Call:
            % S = obj.binSpace('correctEyePos', true, 'window', [-2 2], ...
            %       'binSize', 1);
            %
            % % plot a frame
            %   imagesc(S.xax, S.yax, reshape(X(1,:), S.sz));
            
            
            ip = inputParser();
            ip.addParameter('trialIdx', 1:h.numTrials)
            ip.addParameter('nTimeLags', 30)
            ip.addParameter('correctEyePos', true)
            ip.addParameter('window', [-5 5])
            ip.addParameter('binSize', 1)
            ip.parse(varargin{:})
            

            flipTimes = cell2mat(arrayfun(@(x) x.frameTimes, h.trial(ip.Results.trialIdx), 'UniformOutput', false))';
            
            ppd = h.display.ppd;
            
            win = ip.Results.window*ppd;

            binSize = ip.Results.binSize*ppd;
            binx = @(x) (x==0) + ceil(x/binSize);
            biny = @(x) ceil(x/binSize);
            
            % --- find total stimulus space
            sz = [biny(diff(win)) binx(diff(win))];
            nTotalFrames = sum(arrayfun(@(x) numel(x.frameTimes), h.trial));           
            
            
            h.design.trialIdx = ip.Results.trialIdx;
            h.design.nkTime   = ip.Results.nTimeLags;
            
            
            iFrame = 0;
            
            X = zeros(nTotalFrames, prod(sz));
            

            nTrials = numel(h.trial);
            for kTrial = 1:nTrials
                fprintf('%d / %d \n', kTrial, nTrials)
                frameIdx = iFrame + (1:numel(h.trial(kTrial).frameTimes));
                
                nFrames = numel(frameIdx);
                
                xtmp = zeros(nFrames, prod(sz));
                
                x = squeeze(mean(h.trial(kTrial).pos([1 3],:,:)))';
                y = squeeze(mean(h.trial(kTrial).pos([2 4],:,:)))';
                
                eyeX = h.trial(kTrial).eyePosAtFrame(:,1);
                eyeY = h.trial(kTrial).eyePosAtFrame(:,2);
                
                x_ = bsxfun(@minus, x, eyeX);
                y_ = bsxfun(@minus, y, eyeY);
                
                % flip y, because pixels count from top left to bottom
                % right
                y_ = -y_;
                off = (x_ < win(1) | x_ > win(2)) | (y_ < win(1) | y_ > win(2));
                x_(off) = nan;
                y_(off) = nan;
                
                figure(1); clf
                xbin = biny(x_ - win(1));
                ybin = biny(y_ - win(1));
                               
                for k = 1:size(xbin,2)
                    stimOn = find(~(isnan(xbin(:,k)) | isnan(ybin(:,k))));
                
                    xyind = sub2ind(sz, ybin(stimOn,k), xbin(stimOn,k));
                    
                    ind = sub2ind([nFrames, prod(sz)], stimOn, xyind);
                    
                    xtmp(ind) = xtmp(ind) + 1;
                    
                end
                
                X(frameIdx,:) = xtmp;               
                
                iFrame = iFrame + nFrames;
                
            end
            
            S.X    = X;
            S.bins = flipTimes(:);
            S.size = sz;
            S.xax  = (1:sz(2)) * binSize / ppd + ip.Results.window(1);
            S.yax  = (1:sz(1)) * binSize / ppd + ip.Results.window(1);
            
%             tic
%             fprintf('Unwrapping convolution for STA, regression... \t')
%             Xd  = rfmap.makeStimRowsDense(X, h.design.nkTime);
%             fprintf('[%02.0fs]\n', toc)
%             
%             h.design.biasCol = size(Xd,2)+1;
%             
%             h.design.rowTimes = flipTimes(:);
%             
%             h.design.Xd  = [Xd ones(size(X,1),1)];
%             
%             % sample covariance matrix (for regression)
%             h.design.XX = h.design.Xd'*h.design.Xd;
%             
%             fprintf('Done\n')
%             fprintf('%02.0f total seconds of stimulus\n', size(h.design.Xd,1)*h.display.ifi)

        end
        
        function [Xd, XX] = buildDesignMatrix(h, S, nkTime)
            
            Xd  = rfmap.makeStimRowsDense(S.X, nkTime);           
            
            Xd  = [Xd ones(size(S.X,1),1)];
            
            % sample covariance matrix (for regression)
            XX = Xd'*Xd;
            
        end
        
        function sta = spikeTriggeredAverage(h, S, Xd, spikeTimes, rho)
            
            if nargin < 5
                rho = 10e2;
            end
            
           bins = S.bins;

           % bin spikes at the frame rate
            y = histc(spikeTimes, bins);
            
            % remove spikes from large gaps in frames (hacky way to get rid
            % of trial gaps)
            y(diff(bins)>h.display.ifi) = 0;
            
            % rehsape into column vector
            y = y(:);
            
            sta.w  = Xd'*y;
            sta.xy = sta.w;
            sta.yy = y'*y;
            sta.ny = numel(y);
            
            XX = Xd'*Xd;
            sta.w = (XX + rho * speye(size(XX,2)))\sta.w;
            sta.fullRF  = full(reshape(sta.w(1:end-1), [], prod(S.size)));
            
            % do SVD to get space-time separable RF
            [u,~,v] = svd(sta.fullRF);
            % align temporal RF to 0 (offsets due to bias)
            u(:,1) = u(:,1) - mean(u(1:5,1));
            
            [~, im] = max(abs(u(:,1)));
            sflip = sign(u(im,1));
    
            sta.RF = reshape(sflip*v(:,1), S.size);
            sta.time = (1:numel(u(:,1)))*h.display.ifi;
            sta.RFtime = flipud(sflip*u(:,1));

            
        end
        
        
        function sta = AsdRf(h, spikeTimes)
            
           bins = (1:size(h.design.Xd))*h.display.ifi + h.design.rowTimes(1);
           % bin spikes at the frame rate
            y = histc(spikeTimes, bins);
            
            y = y(:);
            
            sta.w = fastASD_2D_hyper_dual(h.design.Xd(:,setdiff(1:size(h.design.Xd,2), h.design.biasCol)), y-mean(y), [h.design.nkTime h.design.nkx*h.design.nky], 1);
            sta.fullRF  = full(reshape(sta.w, [h.design.nkTime h.design.nkx*h.design.nky]));
            
            [u,~,v] = svd(sta.fullRF);
            u(:,1) = u(:,1) - mean(u(1:5,1));
            [~, im] = max(abs(u(:,1)));
            sflip = sign(u(im,1));
    
            sta.kxs = h.design.kxs;
            sta.kys = h.design.kys;
            sta.RF = reshape(sflip*v(:,1), h.design.nkx, h.design.nky);
            sta.time = (1:h.design.nkTime)*h.display.ifi;
            sta.RFtime = sflip*u(:,1);
            
        end

           
        
    end
end

