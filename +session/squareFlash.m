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
            
            hasStim = hasStim | io.findPDScontainingStimModule(PDS, 'spatialSquares');
            
            if ~any(hasStim)
                return
            end
            
            h.display = PDS{find(hasStim,1)}.initialParametersMerged.display;
                                    
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
            
            
            flipTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), h.trial(ip.Results.trialIdx), 'UniformOutput', false))';
            
            ppd = h.display.ppd;
            
            win = ip.Results.window*ppd;
            if numel(win) == 4
                winx = win([1 3]);
                winy = win([2 4]);
            else
                winx = win;
                winy = win;
            end
            
            binSize = ip.Results.binSize*ppd;
            
            % spatial coordinates (in pixels)
            xax = winx(1):binSize:winx(2);
            yax = winy(1):binSize:winy(2);
            [xx,yy] = meshgrid(xax, yax);
            
%             binx = @(x) (x==0) + ceil(x/binSize);
%             biny = @(x) ceil(x/binSize);
            
            % --- find total stimulus space
%             sz = [biny(diff(winy)) binx(diff(winx))];
            
            sz = size(xx);
            nTotalFrames = sum(arrayfun(@(x) numel(x.frameTimes), h.trial));
            
            
            h.design.trialIdx = ip.Results.trialIdx;
            h.design.nkTime   = ip.Results.nTimeLags;
            
            
            frameCounter = 0;
            
            X = zeros(nTotalFrames, prod(sz));
            
            
            nTrials = numel(h.trial);
            for kTrial = 1:nTrials
                fprintf('%d / %d \n', kTrial, nTrials)
                frameIdx = frameCounter + (1:numel(h.trial(kTrial).frameTimes));
                
                nFrames = numel(frameIdx);
                
                xtmp = zeros(nFrames, prod(sz));
                
%                 x = squeeze(mean(h.trial(kTrial).pos([1 3],:,:)))';
%                 y = squeeze(mean(h.trial(kTrial).pos([2 4],:,:)))';
                
                % get upper left corner of the squares
                xUL = squeeze(h.trial(kTrial).pos(1,:,:))';
                yUL = squeeze(h.trial(kTrial).pos(2,:,:))';
                % lower right corner
                xLR = squeeze(h.trial(kTrial).pos(3,:,:))';
                yLR = squeeze(h.trial(kTrial).pos(4,:,:))';
                
                % eye position at frame flip (in pixels)
                eyeX = h.trial(kTrial).eyePosAtFrame(:,1);
                eyeY = h.trial(kTrial).eyePosAtFrame(:,2);
                
                if ip.Results.correctEyePos
%                     x_ = bsxfun(@minus, x, eyeX);
%                     y_ = bsxfun(@minus, y, eyeY);
                    
                    xUL_ = bsxfun(@minus, xUL, eyeX);
                    yUL_ = bsxfun(@minus, yUL, eyeY);
                    xLR_ = bsxfun(@minus, xLR, eyeX);
                    yLR_ = bsxfun(@minus, yLR, eyeY);
                else
%                     x_ = x - h.display.ctr(1);
%                     y_ = y - h.display.ctr(2);

                    xUL_ = xUL - h.display.ctr(1);
                    yUL_ = yUL - h.display.ctr(2);
                    xLR_ = xLR - h.display.ctr(1);
                    yLR_ = yLR - h.display.ctr(2);
                end
                
                % flip y, because pixels count from top left to bottom
                % right
%                 y_ = -y_;
                yUL_ = -yUL_;
                yLR_ = -yLR_;
                
%                 figure(1); clf
%                 plot(xx(:), yy(:), '.k'); hold on
%                 for iSquare = 1:size(xUL_,2)
%                     for iFrame = 1:size(xUL_,1)
%                         plot([xUL_(iFrame,iSquare) xLR_(iFrame,iSquare)], [yUL_(iFrame,iSquare) yLR_(iFrame,iSquare)], '-'); hold on
%                     end
%                 end
%                 
% %                 figure(2); clf
%                 tmp2 = zeros(size(xx));
%                 for iSquare = 1:size(xUL_,2)
%                     for iFrame = 1:size(xUL_,1)
%                         tmp = (xUL_(iFrame,iSquare) <= xx) & (xLR_(iFrame,iSquare) >= xx) ...
%                             & (yUL_(iFrame,iSquare) >= yy) & (yLR_(iFrame,iSquare) <= yy);
%                             
% %                          tmp =    (yUL_(iFrame,iSquare) >= yy) & (yLR_(iFrame,iSquare) <= yy);
%                         tmp2 = tmp2 + tmp;
%                 
%                     end
%                 end
%                 imagesc(tmp2)
%                 drawnow
                
%                 
                
                for iSquare = 1:size(xUL_,2)
                    for iFrame = 1:size(xUL_,1)
                        tmp = (xUL_(iFrame,iSquare) <= xx) & (xLR_(iFrame,iSquare) >= xx) ...
                            & (yUL_(iFrame,iSquare) >= yy) & (yLR_(iFrame,iSquare) <= yy);
                        xtmp(iFrame,:) = xtmp(iFrame,:) + tmp(:)';
                    end
                end
%                 
                figure(1); clf
                subplot(1,2,1)
                imagesc(reshape(sum(xtmp), sz))
                subplot(1,2,2)
                imagesc(xtmp)
                drawnow
                
%                 % ignore values that are outside the window
%                 off = (x_ < winx(1) | x_ > winx(2)) | (y_ < winy(1) | y_ > winy(2));
%                 x_(off) = nan;
%                 y_(off) = nan;
%                 
%                 figure(1); clf
%                 xbin = binx(x_ - winx(1));
%                 ybin = biny(y_ - winy(1));
%                 
%                 for k = 1:size(xbin,2)
%                     stimOn = find(~(isnan(xbin(:,k)) | isnan(ybin(:,k))));
%                     
%                     xyind = sub2ind(sz, ybin(stimOn,k), xbin(stimOn,k));
%                     
%                     ind = sub2ind([nFrames, prod(sz)], stimOn, xyind);
%                     
%                     xtmp(ind) = xtmp(ind) + 1;
%                     
%                 end
                
                X(frameIdx,:) = xtmp;
                
                frameCounter = frameCounter + nFrames;
                
            end
            
            S.X    = X;
            S.bins = flipTimes(:);
            S.size = sz;
            S.xax  = (1:sz(2)) * binSize / ppd + winx(1)/ppd;
            S.yax  = (1:sz(1)) * binSize / ppd + winy(1)/ppd;
            
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
        
        function D = lssmooth(h, S, Xd, spikeTimes, rho)
            
            if nargin < 5
                rho = 10e2;
            end
            
            sz = size(Xd);
            nt = round((sz(2)-1)/prod(S.size));
            
            inds = reshape(1:(sz(2)-1), nt, [])';
            indsspace = inds(:);
            
            if ~iscell(spikeTimes)
                spikeTimes = {spikeTimes};
            end
            
            for i = 1:sz(1)
                Xd(i,:) = [Xd(i,indsspace) Xd(i,end)];
            end
            
            XX = Xd'*Xd;
            q = qfsmooth(S.size(1), S.size(2));
            
            q = repmat({q}, nt, 1);
            Cinv = blkdiag(q{:}, .01);
            Cinv = Cinv + .001*rho* speye(size(XX,2));
            
            bins = S.bins;
            nUnits = numel(spikeTimes);
            D = cell(nUnits,1);
            
            for kUnit = 1:nUnits
                % bin spikes at the frame rate
                y = histc(spikeTimes{kUnit}, bins);
                
                % remove spikes from large gaps in frames (hacky way to get rid
                % of trial gaps)
                y(diff(bins)>h.display.ifi) = 0;
                
                % rehsape into column vector
                y = y(:);
                
                sta.w  = Xd'*y;
                sta.xy = sta.w;
                sta.yy = y'*y;
                sta.ny = numel(y);
                
                
                
                sta.w = (XX + rho * Cinv )\sta.xy;
                %                 figure(1); clf; plot(sta.w)
                
                sta.fullRF  = full(reshape(sta.w(1:end-1), prod(S.size), nt))';
                
                % do SVD to get space-time separable RF
                [u,~,v] = svd(sta.fullRF);
                % align temporal RF to 0 (offsets due to bias)
                u(:,1) = u(:,1) - mean(u(1:5,1));
                
                [~, im] = max(abs(u(:,1)));
                sflip = sign(u(im,1));
                
                sta.RF = reshape(sflip*v(:,1), S.size);
                sta.time = (1:numel(u(:,1)))*h.display.ifi;
                sta.RFtime = flipud(sflip*u(:,1));
                
                D{kUnit} = sta;
            end
            %             figure(1); clf
            %             subplot(1,2,1)
            %             imagesc(sta.RF)
            %             subplot(1,2,2)
            %             plot(sta.time, sta.RFtime)
            
            
        end
        
        function D = AutoRidgeSpaceTimeSep(h, stim, spikeTimes, maxiter)
            if nargin < 4
                maxiter = 200;  % minimum length scale along each dimension
            end
            
            bins = stim.bins;
            
            for kUnit = 1:numel(spikeTimes)
                y = histc(spikeTimes{kUnit}, bins-h.display.ifi);
                
                % remove spikes from large gaps in frames (hacky way to get rid
                % of trial gaps)
                y(diff(bins)>h.display.ifi) = 0;
                
                % rehsape into column vector
                y = y(:);
                y = y-mean(y);
                ys = smooth(y, 20); % smooth y to start
                
                sta.w  = stim.X'*ys;
                sta.xy = sta.w;
                sta.yy = ys'*ys;
                sta.ny = numel(ys);
                
                
                fprintf('\n\n...Running Automatic ridge regression on space...\n');
                
                nks = [stim.size];
                dd.xx = stim.X'*stim.X;   % stimulus auto-covariance
                dd.xy = (stim.X'*ys); % stimulus-response cross-covariance
                dd.yy = ys'*ys;   % marginal response variance
                dd.nx = numel(dd.xy);     % number of dimensions in stimulus
                dd.ny = numel(ys);  % total number of samples
                
                % Run ridge regression using fixed-point update of hyperparameters
                tic;
                kridge = autoRidgeRegress_fp(dd,maxiter);
                toc;
                
                %                 spaceRF = reshape(kasd, nks);
                
                xtime = stim.X*kridge(:);
                
                nktime = 20;
                Xd = rfmap.makeStimRowsDense(xtime, nktime);
                
                cc.xx = Xd'*Xd;   % stimulus auto-covariance
                cc.xy = (Xd'*y); % stimulus-response cross-covariance
                cc.yy = y'*y;   % marginal response variance
                cc.nx = numel(cc.xy);     % number of dimensions in stimulus
                cc.ny = numel(y);  % total number of samples
                
                ktime = autoRidgeRegress_fp(cc,maxiter);
                
                ktime = flipud(ktime(:));
                
                yu = flipud(y);
                yu = filter(ktime/sum(ktime), 1, yu);
                yu = flipud(yu);
                
                dd.xy = (stim.X'*yu); % stimulus-response cross-covariance
                dd.yy = yu'*yu;   % marginal response variance
                
                % Run ridge regression using fixed-point update of hyperparameters
                tic;
                kridge = autoRidgeRegress_fp(dd,maxiter);
                toc;
                
                sta.RF = reshape(kridge, nks);
                sta.time = (1:nktime)*h.display.ifi;
                sta.RFtime = ktime;
                
                D{kUnit} = sta;
                
            end
        end
        
        %         function sta = glmsmooth(h,S,Xd,spikeTimes, rho)
        %
        %         end
        
        function D = LsSmoothSpaceTimeSep(h, stim, spikeTimes, rho)
            if nargin < 4
                rho = 200;  % minimum length scale along each dimension
            end
            
            bins = stim.bins;
            
            for kUnit = 1:numel(spikeTimes)
                y = histc(spikeTimes{kUnit}, bins);
                
                % remove spikes from large gaps in frames (hacky way to get rid
                % of trial gaps)
                y(diff(bins)>h.display.ifi) = 0;
                
                % rehsape into column vector
                y = y(:);
                %                 y = y-mean(y);
                y = zscore(y);
                ys = smooth(y, 20); % smooth y to start
                ny = numel(y);
                
                sta.w  = stim.X'*ys;
                sta.xy = sta.w;
                sta.yy = ys'*ys;
                sta.ny = numel(ys);
                
                
                %                 fprintf('\n\n...Running Automatic ridge regression on space...\n');
                
                nks = [stim.size];
                
                qf = qfsmooth(nks(1), nks(2));
                
                % %                 folds = getcvfolds(numel(y), 5, 10001);
                % %                 X = [stim.X ones(ny,1)];
                % %                 opts = struct();
                % %                 opts.family = 'poissexp';
                % %
                % %                 kspace = cvglmfitqp(round(ys), X, qf, folds, opts);
                
                dd.xx = stim.X'*stim.X;   % stimulus auto-covariance
                dd.xy = (stim.X'*ys); % stimulus-response cross-covariance
                dd.yy = ys'*ys;   % marginal response variance
                dd.nx = numel(dd.xy);     % number of dimensions in stimulus
                dd.ny = numel(ys);  % total number of samples
                
                kspace = (dd.xx + rho*qf)\dd.xy;
                
                xtime = stim.X*kspace(:);
                
                nktime = 20;
                Xd = rfmap.makeStimRowsDense(xtime, nktime);
                
                cc.xx = Xd'*Xd;   % stimulus auto-covariance
                cc.xy = (Xd'*y); % stimulus-response cross-covariance
                cc.yy = y'*y;   % marginal response variance
                cc.nx = numel(cc.xy);     % number of dimensions in stimulus
                cc.ny = numel(y);  % total number of samples
                
                qftime = .1*qfsmooth1D(nktime) + .9*eye(nktime);
                ktime = (cc.xx + (rho)*qftime)\cc.xy;
                
                ktime = flipud(ktime(:));
                
                yu = flipud(y);
                yu = filter(ktime/sum(ktime), 1, yu);
                yu = flipud(yu);
                
                dd.xy = (stim.X'*yu); % stimulus-response cross-covariance
                dd.yy = yu'*yu;   % marginal response variance
                
                kspace = (dd.xx + rho*qf)\dd.xy;
                
                sta.RF = reshape(kspace, nks);
                sta.time = (1:nktime)*h.display.ifi;
                sta.RFtime = ktime;
                
                D{kUnit} = sta;
                
            end
        end
        
        function D = AsdSpaceTimeSep(h, stim, spikeTimes, minlen)
            if nargin < 4
                minlen = 2.5;  % minimum length scale along each dimension
            end
            
            bins = stim.bins;
            
            for kUnit = 1:numel(spikeTimes)
                y = histc(spikeTimes{kUnit}, bins-h.display.ifi);
                
                % remove spikes from large gaps in frames (hacky way to get rid
                % of trial gaps)
                y(diff(bins)>h.display.ifi) = 0;
                
                % rehsape into column vector
                y = y(:);
                y = y-mean(y);
                ys = smooth(y, 20); % smooth y to start
                
                sta.w  = stim.X'*ys;
                sta.xy = sta.w;
                sta.yy = ys'*ys;
                sta.ny = numel(ys);
                
                
                fprintf('\n\n...Running ASD_2D on space...\n');
                
                nks = [stim.size];
                tic;
                kasd = fastASD(stim.X,ys,nks,minlen);
                toc;
                
                %                 spaceRF = reshape(kasd, nks);
                
                xtime = stim.X*kasd(:);
                
                nktime = 20;
                Xd = rfmap.makeStimRowsDense(xtime, nktime);
                
                [ktime,~] = fastASD(Xd,y,[nktime 1],minlen);
                
                ktime = flipud(ktime(:));
                
                yu = flipud(y);
                yu = filter(ktime/sum(ktime), 1, yu);
                yu = flipud(yu);
                
                [kasd,~] = fastASD(stim.X,yu,nks,minlen);
                
                sta.RF = reshape(kasd, nks);
                sta.time = (1:nktime)*h.display.ifi;
                sta.RFtime = ktime;
                
                D{kUnit} = sta;
                
            end
        end
        
        function D = AsdRf(h, stim, Xd, spikeTimes, minlen)
            
            if nargin < 5
                minlen = 2.5;  % minimum length scale along each dimension
            end
            
            sz = size(Xd);
            nt = round((sz(2)-1)/prod(stim.size));
            bins = stim.bins;
            
            for kUnit = 1:nUnits
                % bin spikes at the frame rate
                y = histc(spikeTimes{kUnit}, bins);
                
                % remove spikes from large gaps in frames (hacky way to get rid
                % of trial gaps)
                y(diff(bins)>h.display.ifi) = 0;
                
                % rehsape into column vector
                y = y(:);
                y = y-mean(y);
                sta.w  = Xd'*y;
                sta.xy = sta.w;
                sta.yy = y'*y;
                sta.ny = numel(y);
                
                
                fprintf('\n\n...Running ASD_2D...\n');
                
                
                
                nks = [nt stim.size];
                tic;
                [kasd,asdstats] = fastASD(Xd(:,1:end-1),y,nks,minlen);
                toc;
                
                kasd_tns = reshape(kasd,nks(1), prod(nks(2:3)));
                
                for j = 1:min(4,nks(3));
                    subplot(3,4,j); imagesc(ktns(:,:,j));
                    title(sprintf('slice %d',j));
                    subplot(3,4,j+4); imagesc(kridge_tns(:,:,j));
                    subplot(3,4,j+8); imagesc(kasd_tns(:,:,j));
                end
                subplot(3,4,1); ylabel('\bf{true k}');
                subplot(3,4,5); ylabel('\bf ridge');
                subplot(3,4,9); ylabel('\bf ASD');
                
                
                sta.w = kasd;
                %                 figure(1); clf; plot(sta.w)
                
                sta.fullRF  = full(reshape(sta.w, nt, prod(stim.size)));
                
                % do SVD to get space-time separable RF
                [u,~,v] = svd(sta.fullRF);
                % align temporal RF to 0 (offsets due to bias)
                u(:,1) = u(:,1) - mean(u(1:5,1));
                
                [~, im] = max(abs(u(:,1)));
                sflip = sign(u(im,1));
                
                sta.RF = reshape(sflip*v(:,1), stim.size);
                sta.time = (1:numel(u(:,1)))*h.display.ifi;
                sta.RFtime = flipud(sflip*u(:,1));
                
                figure(1); clf
                subplot(1,2,1)
                imagesc(sta.RF)
                subplot(1,2,2)
                plot(sta.time, sta.RFtime)
                
                
                D{kUnit} = sta;
            end
            
        end
        
        %
        %             bins = (1:size(h.design.Xd))*h.display.ifi + h.design.rowTimes(1);
        %            % bin spikes at the frame rate
        %             y = histc(spikeTimes, bins);
        %
        %             y = y(:);
        %
        %             fprintf('\n\n...Running ASD_2D...\n');
        %
        %             minlen = 2.5;  % minimum length scale along each dimension
        %
        %             nks = [stim.
        %             tic;
        %             [kasd,asdstats] = fastASD(stim.x,y,nks,minlen);
        %             toc;
        %
        %
        %             sta.w = fastASD_2D_hyper_dual(h.design.Xd(:,setdiff(1:size(h.design.Xd,2), h.design.biasCol)), y-mean(y), [h.design.nkTime h.design.nkx*h.design.nky], 1);
        %             sta.fullRF  = full(reshape(sta.w, [h.design.nkTime h.design.nkx*h.design.nky]));
        %
        %             [u,~,v] = svd(sta.fullRF);
        %             u(:,1) = u(:,1) - mean(u(1:5,1));
        %             [~, im] = max(abs(u(:,1)));
        %             sflip = sign(u(im,1));
        %
        %             sta.kxs = h.design.kxs;
        %             sta.kys = h.design.kys;
        %             sta.RF = reshape(sflip*v(:,1), h.design.nkx, h.design.nky);
        %             sta.time = (1:h.design.nkTime)*h.display.ifi;
        %             sta.RFtime = sflip*u(:,1);
        %
        %         end
        
        
        
    end
    
    methods (Static)
        
        function [trial, display] = importPDS(PDS)
            % importPDS checks which version of to stimulus code was run
            % and imports to a common format appropriately
            
            pdsDate = PDS.initialParametersMerged.session.initTime;
            if isfield(PDS.initialParametersMerged.git, 'pep')
                
                if any(strfind(PDS.initialParametersMerged.git.pep.status, 'branch cleanup'))
                    
                    if pdsDate > datenum(2018, 02, 01)
                        [trial, display] = session.squareFlash.importPDS_v2(PDS);
                    else
                        error('unknown version')
                    end
                    
                else
                    error('unknown version')
                end
            else
                [trial, display] = session.squareFlash.importPDS_v1(PDS);
            end
            
        end
        
        function [trial, display] = importPDS_v2(PDS)
            
            trial   = [];
            display = PDS.initialParametersMerged.display;
            
            stim = 'spatialSquares';
            
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
                return;
            end
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end)); %#ok<*AGROW>
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end)) - trial(kTrial).start;
                
                trial(kTrial).pos        = PDS.data{thisTrial}.(stim).pos;
                if size(trial(kTrial).pos, 3) ~= numel(trial(kTrial).frameTimes)
                    keyboard
                end
                eyepos = io.getEyePosition(PDS, thisTrial);
                
                trial(kTrial).eyeSampleTime = eyepos(:,1);
                trial(kTrial).eyeXPx        = eyepos(:,2);
                trial(kTrial).eyeYPx        = eyepos(:,3);
                trial(kTrial).pupilArea     = eyepos(:,4);
                
                % --- additional eye position analyses
                %                     results = pdsa.detectSaccades(eyepos(:,1)', eyepos(:,2:3)'./ppd, 'verbose', false);
                
                
                
                iix = trial(kTrial).eyeXPx < 200 | trial(kTrial).eyeXPx > 1800;
                iiy = trial(kTrial).eyeYPx < 100 | trial(kTrial).eyeXPx > 1000;
                bad = iix | iiy;
                
                trial(kTrial).eyeXPx(bad) = nan;
                trial(kTrial).eyeYPx(bad) = nan;
                trial(kTrial).pupilArea(bad) = nan;
                
%                 % find eye position on each frame
%                 t = trial(kTrial).frameTimes - trial(kTrial).start;
%                 tdiff = abs(bsxfun(@minus, trial(kTrial).eyeSampleTime, t(:)')) < mean(diff(trial(kTrial).eyeSampleTime));
%                 [irow,~] = find(diff(tdiff)==-1);
                
                %                     trial(kTrial).saccades = false(size(trial(kTrial).frameTimes));
                %
                %                     for iSaccade = 1:size(results,2)
                %                         ix = t > results(1, iSaccade) & t < results(2,iSaccade);
                %                         trial(kTrial).saccades(ix) = true;
                %                     end
                
                trial(kTrial).eyePosAtFrame = PDS.data{thisTrial}.behavior.eyeAtFrame';
                
                iix = trial(kTrial).eyePosAtFrame(:,1) < 200 | trial(kTrial).eyePosAtFrame(:,1) > 1800;
                iiy = trial(kTrial).eyePosAtFrame(:,2) < 100 | trial(kTrial).eyePosAtFrame(:,2) > 1000;
                bad = iix | iiy;
                
                trial(kTrial).eyePosAtFrame(bad,:) = nan;
                
            end
            
        end
        
        function [trial, display] = importPDS_v1(PDS)
            
            trial   = [];
            display = PDS.initialParametersMerged.display;
            
            stim = 'SpatialMapping';
            
            trialIx = cellfun(@(x) isfield(x, stim), PDS.data);
            
            stimTrials = find(trialIx);
            
            if isempty(stimTrials)
                return;
            end
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1)); %#ok<*AGROW>
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - trial(kTrial).start;
                
                trial(kTrial).pos        = PDS.data{thisTrial}.(stim).pos;
                
                eyepos = io.getEyePosition(PDS, thisTrial);
                trial(kTrial).eyeSampleTime = eyepos(:,1);
                trial(kTrial).eyeXPx        = eyepos(:,2);
                trial(kTrial).eyeYPx        = eyepos(:,3);
                trial(kTrial).pupilArea     = eyepos(:,4);
                
                % --- additional eye position analyses
                %                     results = pdsa.detectSaccades(eyepos(:,1)', eyepos(:,2:3)'./ppd, 'verbose', false);
                
                
                
                iix = trial(kTrial).eyeXPx < 200 | trial(kTrial).eyeXPx > 1800;
                iiy = trial(kTrial).eyeYPx < 100 | trial(kTrial).eyeXPx > 1000;
                bad = iix | iiy;
                
                trial(kTrial).eyeXPx(bad) = nan;
                trial(kTrial).eyeYPx(bad) = nan;
                trial(kTrial).pupilArea(bad) = nan;
                
                % find eye position on each frame
                t = trial(kTrial).frameTimes - trial(kTrial).start;
                tdiff = abs(bsxfun(@minus, trial(kTrial).eyeSampleTime, t(:)')) < mean(diff(trial(kTrial).eyeSampleTime));
                [irow,~] = find(diff(tdiff)==-1);
                
                %                     trial(kTrial).saccades = false(size(trial(kTrial).frameTimes));
                %
                %                     for iSaccade = 1:size(results,2)
                %                         ix = t > results(1, iSaccade) & t < results(2,iSaccade);
                %                         trial(kTrial).saccades(ix) = true;
                %                     end
                
                trial(kTrial).eyePosAtFrame = [trial(kTrial).eyeXPx(irow) trial(kTrial).eyeYPx(irow)];
                
            end
            
        end
        
    end
end

