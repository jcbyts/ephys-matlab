classdef hartleyFF < handle
% Hartley FF (Full Field)

    properties
        numTrials
        trial
        display
        design
    end
    
    methods
        function h = hartleyFF(PDS, varargin)
            
            stim = 'hartley';
            
            
            hasStim = io.findPDScontainingStimModule(PDS, stim);
            
            if ~any(hasStim)
                return
            end
            
            h.display = PDS{1}.initialParametersMerged.display;
            
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
            
            plot(h.trial(kTrial).eyeXpx); hold on
            plot(h.trial(kTrial).eyeYpx);
            
        end
        
        
        function buildDesignMatrix(h, varargin)
            
            ip = inputParser();
            ip.addParameter('trialIdx', 1:h.numTrials)
            ip.addParameter('nTimeLags', 30)
            ip.parse(varargin{:})
            
            
            flipTimes = cell2mat(arrayfun(@(x) x.frameTimes, h.trial(ip.Results.trialIdx), 'UniformOutput', false))';
            kx        = cell2mat(arrayfun(@(x) x.kx', h.trial(ip.Results.trialIdx), 'UniformOutput', false))';
            ky        = cell2mat(arrayfun(@(x) x.ky', h.trial(ip.Results.trialIdx), 'UniformOutput', false))';
            on        = ~(isnan(kx) | isnan(ky));
            
            % --- find total stimulus space
            kxs = unique(kx(on));
            kys = unique(ky(on));
            
            nTotalFrames = sum(arrayfun(@(x) numel(x.frameTimes), h.trial));
            
            h.design.trialIdx = ip.Results.trialIdx;
            h.design.nkTime   = ip.Results.nTimeLags;
            h.design.nkx      = numel(kxs);
            h.design.nky      = numel(kys);
            h.design.kxs      = kxs;
            h.design.kys      = kys;
            
            
            iFrame = 0;
            
            sz = [h.design.nkx h.design.nky];
            X = zeros(nTotalFrames, prod(sz));
            
            
            
            nTrials = numel(h.trial);
            for kTrial = 1:nTrials
                fprintf('%d / %d \n', kTrial, nTrials)
                frameIdx = iFrame + (1:numel(h.trial(kTrial).frameTimes));
                
                nFrames = numel(frameIdx);
                
                xtmp = zeros(nFrames, prod(sz));
                
                stimOn = find(h.trial(kTrial).on);
                [kxj, ~] = find( bsxfun(@eq, h.trial(kTrial).kx(stimOn), h.design.kxs(:)')');
                tmp = bsxfun(@eq, h.trial(kTrial).ky(stimOn), h.design.kys(:)');
                [kyj, ~] = find( tmp' );
                
                kxy = sub2ind(sz, kxj, kyj);
                
                ind = sub2ind([nFrames prod(sz)], stimOn, kxy);
                
                xtmp(ind) = 1;
                
                X(frameIdx,:) = xtmp;
                
                
                iFrame = iFrame + nFrames;
                
            end
            
            tic
            fprintf('Unwrapping convolution for STA, regression... \t')
            Xd  = rfmap.makeStimRowsDense(X, h.design.nkTime);
            fprintf('[%02.0fs]\n', toc)
            
            h.design.biasCol = size(Xd,2)+1;
            
            h.design.rowTimes = flipTimes(:);
            
            h.design.Xd  = [Xd ones(size(X,1),1)];
            
            % sample covariance matrix (for regression)
            h.design.XX = h.design.Xd'*h.design.Xd;
            
            fprintf('Done\n')
            fprintf('%02.0f total seconds of stimulus\n', size(h.design.Xd,1)*h.display.ifi)
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             
%             x=arrayfun(@(x) find(x==kxs), kx(on));
%             y=arrayfun(@(x) find(x==kys), ky(on));
%             
%             ind=sub2ind([numel(kys) numel(kxs)], y, x);
%             
%             
%             fprintf('Building the Design Matrix for %d trials\n', numel(ip.Results.trialIdx))
%             
%             t0 = flipTimes(1);
%             binsize = h.display.ifi;
%             
%             temporalBinningFunction = @(t) (t==0) + ceil(t/binsize);
%             
%             stimOnset = temporalBinningFunction(flipTimes(on)-t0);
%             
%             X = sparse(stimOnset, ind, ones(size(ind,1),1), max(stimOnset), max(ind));
%             
%             ntk = ip.Results.nTimeLags;
%             Xd  = rfmap.makeStimRowsSparse(X, ntk);
%             h.design.biasCol = size(Xd,2)+1;
%             
%             h.design.rowTimes = flipTimes(:);
%             h.design.binfun   = temporalBinningFunction;
%             
%             h.design.Xd  = [Xd ones(size(X,1),1)];
%             
%             h.design.XX = h.design.Xd'*h.design.Xd;
        end
        
        function sta = spikeTriggeredAverage(h, spikeTimes)
            
%            bins = (1:size(h.design.Xd))*h.display.ifi + h.design.rowTimes(1);
            bins = h.design.rowTimes;
           % bin spikes at the frame rate
            y = histc(spikeTimes, bins);
            
%             figure(1); clf
%             plot(y)
            y(diff(bins)>h.display.ifi) = 0;
            
            y = y(:);
            
            sta.w = (h.design.Xd'*y);
            sta.xy = sta.w;
            sta.yy = y'*y;
            sta.ny = numel(y);
            sta.w = (h.design.XX + 10e2 * speye(size(h.design.XX,2)))\sta.w;
            sta.fullRF  = full(reshape(sta.w(1:end-1), [h.design.nkTime h.design.nkx*h.design.nky]));
            
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
        
        function plotFrozenRaster(h, s, clustId)
            
            if nargin < 3
                clustId = s.cids;
            end
            nUnits = numel(clustId);
            
            hartleyTrial = h.trial;
            frozenTrials = find([h.trial.frozenSequence]);
            if any(frozenTrials)
                tmp = arrayfun(@(x) x.frameTimes(1:x.frozenSequenceLength:end)', hartleyTrial(frozenTrials), 'UniformOutput', false)';
                tmp = cellfun(@(x) x(:), tmp(:), 'uni', false);
                sequenceStarts = cell2mat(tmp);
                
                cmap = lines(nUnits);
                
                for kUnit = 1:nUnits
                    st = s.st(s.clu==clustId(kUnit));
                    
                    seqLength = hartleyTrial(frozenTrials(1)).frozenSequenceLength;
                    
                    [spcnt, bcenters] = pdsa.binSpTimes(st, sequenceStarts, [0 seqLength/120], 1e-3);
                    
                    [i, j] = find(spcnt);
                    
                    uIx = s.cids==clustId(kUnit);
                    
                    if isfield(s, 'cgs') && s.cgs(uIx)==3
                        plot(bcenters(j), i-numel(sequenceStarts)*(kUnit-1), '.', 'Color', cmap(kUnit,:)); hold on
                    else
                        plot(bcenters(j), i-numel(sequenceStarts)*(kUnit-1), '.', 'Color', repmat(.5, 1, 3)); hold on
                    end
                    
                end
            end
            
%             for i = frozenTrials(:)'
%                 
%                plot(h.trial(i).eyeSampleTime,  [h.trial(i).eyeXPx, h.trial(i).eyeYPx], 'k');
%             end
        end
        
           
        
    end
    
    methods (Static)
        
        function [trial, display] = importPDS(PDS)
            
            pdsDate = PDS.initialParametersMerged.session.initTime;
            
            if pdsDate > datenum(2018,02,01)
                trial = [];
                display = [];
            else
                
                [trial, display] = session.hartleyFF.importPDS_v1(PDS);
            end
                
        end
        
        
        function [trial, display] = importPDS_v1(PDS)
            
            trial   = [];
            display = PDS.initialParametersMerged.display;
            
            stim = 'hartley';
            
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
                    
                    if isfield(PDS.conditions{thisTrial}, stim)
                        if isfield(PDS.conditions{thisTrial}.(stim), 'setupRNG')
                            if strcmp(PDS.conditions{thisTrial}.(stim).setupRNG, 'frozenSequence')
                                trial(kTrial).frozenSequence = true;
                                trial(kTrial).frozenSequenceLength = PDS.conditions{thisTrial}.(stim).sequenceLength;
                            else
                                trial(kTrial).frozenSequence = false;
                                trial(kTrial).frozenSequenceLength = nan;
                            end
                            
                        else
                            trial(kTrial).frozenSequence = false;
                            trial(kTrial).frozenSequenceLength = nan;
                        end
                        
                    else
                        trial(kTrial).frozenSequence = false;
                        trial(kTrial).frozenSequenceLength = nan;
                    end
                    
                    trial(kTrial).kx         = PDS.data{thisTrial}.(stim).kx;
                    trial(kTrial).ky         = PDS.data{thisTrial}.(stim).ky;
                    trial(kTrial).on         = PDS.data{thisTrial}.(stim).on;
                    trial(kTrial).phi        = PDS.data{thisTrial}.(stim).phi;
                    trial(kTrial).tf         = PDS.data{thisTrial}.(stim).tf;
                    
                    eyepos = io.getEyePosition(PDS, thisTrial);
                    trial(kTrial).eyeSampleTime = eyepos(:,1);
                    trial(kTrial).eyeXPx        = eyepos(:,2);
                    trial(kTrial).eyeYPx        = eyepos(:,3);
                    trial(kTrial).pupilArea     = eyepos(:,4);
                end
            
        end
        
        
    end
end

