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
            % CONSTRUCTOR
            % h = squareFlash(PDS) creates the squareFlash session object
            % PDS is a PDS struct or cell array of PDS structs
            %
            % Optional Arguments:
            %   'eyetrace' [n x 4]  the raw eye position data
            %                       col 1: timestamp
            %                       col 2: x position (degrees)
            %                       col 3: y position (degrees)
            %                       col 4: pupil
            
            % parse optional arguments
            ip = inputParser();
            ip.addOptional('eyetrace', [])
            ip.addOptional('saccades', [])
            ip.addOptional('eyePosPreTrial', 1.5) % seconds before the trial to store
            ip.parse(varargin{:})
            
            
            stim = 'SpatialMapping'; % this is the name we use during data colletion. It should always be the same!
            
            
            if isstruct(PDS)
                PDS = {PDS};
            end
            
            % check for PDS files that contain this stimulus module
            hasStim = io.findPDScontainingStimModule(PDS, stim);
            hasStim = hasStim | io.findPDScontainingStimModule(PDS, 'spatialSquares');
            
            if ~any(hasStim) % exit program if this stimulus wasn't run
                return
            end
            
            % assuming the display parameters didn't change during the
            % session, store them in the object
            h.display = PDS{find(hasStim,1)}.initialParametersMerged.display;
            
            % Loop over PDS files, importing the stimulus parameters
            for i = find(hasStim(:)')
                
                % helper function for the import
                trial_ = h.importPDS(PDS{i});
                
                if isempty(trial_)
                    continue
                end
                
                h.trial = [h.trial; trial_(:)];
                
            end
            
            
            h.numTrials = numel(h.trial);
            
            % --- import eye position
            eyepos = ip.Results.eyetrace;
            if ~isempty(eyepos)
                if size(eyepos,2) < size(eyepos,1)
                    eyepos = eyepos';
                end
                
                assert(size(eyepos,1) >= 3, 'the first row must be timestamps')
                
                if ~isstruct(ip.Results.saccades)% detect saccades
                    sampleRate = 1/mode(diff(eyepos(1,:)));
                    ix = ~any(isnan(eyepos(2,:)));
                    [saccades] = pdsa.detectSaccades(eyepos(1,ix), eyepos(2:3,ix), ...
                        'verbose', false, ...
                        'filterPosition', 1, ...
                        'filterLength', ceil(20/sampleRate*1e3), ... % 40 ms smoothing for velocity computation
                        'detectThresh', 200, ...
                        'startThresh', 5, ...
                        'minIsi', ceil(30/sampleRate*1e3), ...
                        'minDur', ceil(4/sampleRate*1e3), ... % 4 ms
                        'blinkIsi', ceil(40/sampleRate*1e3));
                else
                    saccades = ip.Results.saccades;
                end
                
                if size(eyepos,1) >= 4
                    pupil = eyepos(4,:);
                else
                    pupil = [];
                end
                
                % Loop over trials and add eye position
                for kTrial = 1:h.numTrials
                    
                    % include time before the trial starts
                    preTrial = ip.Results.eyePosPreTrial;
                    
                    % find valid eye position
                    iix = (eyepos(1,:) > (h.trial(kTrial).start - preTrial)) & (eyepos(1,:) < (h.trial(kTrial).start + h.trial(kTrial).duration));
                    h.trial(kTrial).eyeSampleTime = eyepos(1,iix);
                    h.trial(kTrial).eyeXDeg = eyepos(2,iix);
                    h.trial(kTrial).eyeYDeg = eyepos(3,iix);
                    if ~isempty(pupil)
                        h.trial(kTrial).pupilArea = pupil(iix);
                    end
                    
                    % valid saccade times
                    iix = (saccades.tstart > (h.trial(kTrial).start - preTrial)) & (saccades.tend < (h.trial(kTrial).start + h.trial(kTrial).duration));
                    fields = fieldnames(saccades);
                    for iField = 1:numel(fields)
                        newfield = ['sac_' fields{iField}];
                        h.trial(kTrial).(newfield) = saccades.(fields{iField})(iix);
                    end
                    
                    % update the eye position at frame to use the offline
                    % estimate
                    nFrames = numel(h.trial(kTrial).frameTimes);
                    for iFrame = 1:nFrames
                        iix = h.trial(kTrial).eyeSampleTime > (h.trial(kTrial).frameTimes(iFrame)) & ( h.trial(kTrial).eyeSampleTime < (h.trial(kTrial).frameTimes(iFrame) +  h.display.ifi));
                        if any(iix)
                            h.trial(kTrial).eyePosAtFrame(iFrame,:) = [nanmedian(h.trial(kTrial).eyeXDeg(iix)) nanmedian(h.trial(kTrial).eyeYDeg(iix))];
                        end
                    end
                    
                end
                
            end
            
            
        end % constructor
        
        function spks = binSpikes(h, stim, sp, units)
            % spks = binSpikes(h, stim, sp, units)
            T = size(stim.X,1);
            spks = zeros(T,1);
            if nargin < 4
                units = sp.cids(sp.cgs >= 1);
            end
            
            nUnits = numel(units);
            for i = 1:nUnits
                cnt = histc(sp.st(sp.clu == units(i)), stim.bins);
                cnt((diff(stim.bins) > (h.display.ifi*1.5))) = 0;
                spks(:,i) = cnt;
            end
            
        end % binSpikes
        
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
            
            
        end % plotTrial
        
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
            % S = obj.binSpace('correctEyePos', 'no', 'window', [-2 2], ...
            %       'binSize', 1);
            %
            % % plot a frame
            %   imagesc(S.xax, S.yax, reshape(X(1,:), S.sz));
            
            
            ip = inputParser();
            ip.addParameter('trialIdx', 1:h.numTrials)
            ip.addParameter('correctEyePos', 'no')
            ip.addParameter('window', [-5 5])
            ip.addParameter('binSize', 1)
            ip.addParameter('plot', false)
            ip.addParameter('saccadeBuffer', .050) % 50 ms
            ip.parse(varargin{:})
            
            
            flipTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), h.trial(ip.Results.trialIdx), 'UniformOutput', false))';
            
            win = ip.Results.window; % convert window into pixels
            
            % check if window has the same size for x and y
            if numel(win) == 4
                winx = win([1 3]);
                winy = win([2 4]);
            else
                winx = win;
                winy = win;
            end
            
            binSize = ip.Results.binSize; % convert the bin size to pixels
            
            % spatial coordinates (in pixels)
            xax = winx(1):binSize:winx(2);
            yax = winy(1):binSize:winy(2);
            [xx,yy] = meshgrid(xax, yax);
            
            trialIdx = ip.Results.trialIdx;
            % get the dimensions of the space
            sz = size(xx);
            nTotalFrames = sum(arrayfun(@(x) numel(x.frameTimes), h.trial(trialIdx)));
            
            
            h.design.trialIdx = trialIdx;
            
            frameCounter = 0;
            
            X = zeros(nTotalFrames, prod(sz));
            valid = false(nTotalFrames, 1); % is the frame valid to analyze
            
            nTrials = numel(trialIdx);
            for iiT = 1:nTrials
                fprintf('%d / %d \n', iiT, nTrials)
                
                kTrial = trialIdx(iiT);
                
                % index into X
                frameIdx = frameCounter + (1:numel(h.trial(kTrial).frameTimes));
                
                % number of frames this trial
                nFrames = numel(frameIdx);
                
                % X for this trial
                Xtrial = zeros(nFrames, prod(sz));
                
                % get upper left corner of the squares
                xUL = squeeze(h.trial(kTrial).pos(1,:,:))';
                yUL = squeeze(h.trial(kTrial).pos(2,:,:))';
                
                % lower right corner
                xLR = squeeze(h.trial(kTrial).pos(3,:,:))';
                yLR = squeeze(h.trial(kTrial).pos(4,:,:))';
                
                % eye position at frame flip (in pixels)
                eyeX = h.trial(kTrial).eyePosAtFrame(:,1);
                eyeY = h.trial(kTrial).eyePosAtFrame(:,2);
                
                switch ip.Results.correctEyePos
                    case 'simple'
                        xUL_ = bsxfun(@minus, xUL, eyeX);
                        yUL_ = bsxfun(@minus, yUL, eyeY);
                        xLR_ = bsxfun(@minus, xLR, eyeX);
                        yLR_ = bsxfun(@minus, yLR, eyeY);
                        Vtrial = true(size(xLR_,1),1);
                    case 'no'
                        xUL_ = xUL;
                        yUL_ = yUL;
                        xLR_ = xLR;
                        yLR_ = yLR;
                        Vtrial = true(size(xLR_,1),1);
                    case 'centralNo'
                        xUL_ = xUL;
                        yUL_ = yUL;
                        xLR_ = xLR;
                        yLR_ = yLR;
                        
                        eyeDist = sqrt( (eyeX).^2 + (eyeY).^2);
                        Vtrial = eyeDist(:) < 5; % central 5 degrees
                        
                    case 'centralYes'
                        xUL_ = bsxfun(@minus, xUL, eyeX);
                        yUL_ = bsxfun(@minus, yUL, eyeY);
                        xLR_ = bsxfun(@minus, xLR, eyeX);
                        yLR_ = bsxfun(@minus, yLR, eyeY);
                        
                        eyeDist = sqrt( (eyeX).^2 + (eyeY).^2);
                        Vtrial = eyeDist(:) < 6;
                        
                end
                
                % flip y, because pixels count from top left to bottom
                % right
                %                 yUL_ = -yUL_;
                %                 yLR_ = -yLR_;
                
                % loop over squares and recreate binned stimulus
                for iSquare = 1:size(xUL_,2)
                    for iFrame = 1:size(xUL_,1)
                        tmp = (xUL_(iFrame,iSquare) <= xx) & (xLR_(iFrame,iSquare) >= xx) ...
                            & (yUL_(iFrame,iSquare) >= yy) & (yLR_(iFrame,iSquare) <= yy);
                        Xtrial(iFrame,:) = Xtrial(iFrame,:) + tmp(:)';
                    end
                end
                %
%                                 figure(1); clf
%                                 subplot(1,2,1)
%                                 imagesc(reshape(sum(Xtrial), sz))
%                                 subplot(1,2,2)
%                                 imagesc(Xtrial)
%                                 drawnow
%                 
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
                %                     Xtrial(ind) = Xtrial(ind) + 1;
                %
                %                 end
                
                % if clipping out saccades
                if isfield(h.trial, 'sac_start')
                    sac_buffer = ip.Results.saccadeBuffer;
                    sacStartInds = bsxfun(@gt, h.trial(kTrial).frameTimes(:), h.trial(kTrial).sac_start(:)');
                    sacStopInds  = bsxfun(@lt, h.trial(kTrial).frameTimes(:), h.trial(kTrial).sac_end(:)'+sac_buffer);
                    fixations = ~any(sacStartInds & sacStopInds, 2);
                    iix = fixations(:) & Vtrial(:);
                    X(frameIdx(iix),:) = Xtrial(iix,:);
                    valid(frameIdx(iix)) = true;
                else
                    iix = Vtrial(:);
                    X(frameIdx(iix),:) = Xtrial(iix,:);
                    valid(frameIdx(iix)) = true;
                end
                
                frameCounter = frameCounter + nFrames;
                
            end
            
            S.X    = X;
            S.valid = valid;
            S.bins = flipTimes(:);
            S.size = sz;
            S.xax  = xax + binSize/2; % recenter bins
            S.yax  = yax + binSize/2; %(1:sz(1)) * binSize + winy(1);
            
        end
        
        function [Xd, XX] = buildDesignMatrix(~, S, nkTime)
            
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
                
                cc.xx = Xd'*Xd;     % stimulus auto-covariance
                cc.xy = (Xd'*y);    % stimulus-response cross-covariance
                cc.yy = y'*y;       % marginal response variance
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
                    warning('squareFlash: git tracking failed. Trying import version 2')
                    try
                        [trial, display] = session.squareFlash.importPDS_v2(PDS);
                    catch
                        error('Version 2 import failed')
                    end
                end
            else
                [trial, display] = session.squareFlash.importPDS_v1(PDS);
            end
            
            % --- convert pos to degrees
            for iTrial = 1:numel(trial)
                nSquares = size(trial(iTrial).pos,2);
                for iSquare = 1:nSquares
                    
                    % get position in pixels
                    pospx = squeeze(trial(iTrial).pos(1:4,iSquare,:));
                    
                    % subtract the center of the screen
                    pospx = bsxfun(@minus, pospx, display.ctr');
                    
                    % flip y pos so up is positive
                    pospx(2,:) = -pospx(2,:);
                    pospx(4,:) = -pospx(4,:);
                    
                    % convert to degrees
                    trial(iTrial).pos(1:2,iSquare,:) = pds.px2deg(pospx(1:2,:), display.viewdist, display.px2w);
                    trial(iTrial).pos(3:4,iSquare,:) = pds.px2deg(pospx(3:4,:), display.viewdist, display.px2w);
                end
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
            
            pdsTrial = pds.getPdsTrialData(PDS);
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end)); %#ok<*AGROW>
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end)) - trial(kTrial).start;
                
                if j==1
                    trial = mergeStruct(trial(kTrial), pdsTrial(thisTrial).(stim));
                else
                    trial(kTrial) = mergeStruct(trial(kTrial), pdsTrial(thisTrial).(stim));
                end
                %                 trial(kTrial).pos        = PDS.data{thisTrial}.(stim).pos;
                if size(trial(kTrial).pos, 3) ~= numel(trial(kTrial).frameTimes)
                    keyboard
                end
                
                % Eye position at the frame time
                nFrames = numel(trial(kTrial).frameTimes);
                
                % Eye position in pixels
                trial(kTrial).eyePosAtFrame = pdsTrial(thisTrial).behavior.eyeAtFrame(:,1:nFrames)';
%                 trial(kTrial).eyePosAtFrame = PDS.data{thisTrial}.behavior.eyeAtFrame(:,1:nFrames)';
                
                % Exclude edges of the screen
                iix = trial(kTrial).eyePosAtFrame(:,1) < 200 | trial(kTrial).eyePosAtFrame(:,1) > 1800;
                iiy = trial(kTrial).eyePosAtFrame(:,2) < 100 | trial(kTrial).eyePosAtFrame(:,2) > 1000;
                bad = iix | iiy;
                
                trial(kTrial).eyePosAtFrame(bad,:) = nan;
                
                % subtract the center of the screen
                pospx = bsxfun(@minus, trial(kTrial).eyePosAtFrame, display.ctr(1:2));
                
                % flip y pos so up is positive
                pospx(:,2) = -pospx(:,2);
                
                % convert to degrees
                trial(kTrial).eyePosAtFrame = pds.px2deg(pospx', display.viewdist, display.px2w)';
                
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
            
            pdsTrial = pds.getPdsTrialData(PDS);
            
            for j = 1:numel(stimTrials)
                thisTrial = stimTrials(j);
                
                kTrial = j;
                
                trial(kTrial).frameTimes = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,1:end-1)); %#ok<*AGROW>
                trial(kTrial).start      = trial(kTrial).frameTimes(1);
                trial(kTrial).duration   = PDS.PTB2OE(PDS.data{thisTrial}.timing.flipTimes(1,end-1)) - trial(kTrial).start;
                
                % convert pos to degrees
                %                 trial(kTrial).pos        = PDS.data{thisTrial}.(stim).pos;
                if j==1
                    trial = mergeStruct(trial(kTrial), pdsTrial(kTrial).(stim));
                else
                    trial(kTrial) = mergeStruct(trial(kTrial), pdsTrial(thisTrial).(stim));
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

