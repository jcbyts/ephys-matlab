function [Coarse, F] = CoarseToFineRF(PDS, spikes)
% checks if spatial RF mapping was run, uses a coarse initialization of the
% RF to measure
%%
spatialMap = session.squareFlash(PDS);

%%
SPATIAL_MAPPING_WINDOW        = [-15 -10 15 10];

% SPATIAL_MAPPING_BIN_SIZE      = 1;
% SPATIAL_MAPPING_WINDOW        = [-1 -1 1 1];

SPATIAL_MAPPING_TEMPORAL_LAGS = 1:8;

squareSizesTrial = arrayfun(@(x) x.size, spatialMap.trial);
squareSizes = unique(squareSizesTrial);

SPATIAL_MAPPING_BIN_SIZE      = max(squareSizes);

trialIdx = find(squareSizesTrial==SPATIAL_MAPPING_BIN_SIZE);

SPATIAL_MAPPING_BIN_SIZE      = SPATIAL_MAPPING_BIN_SIZE  * 2;

% binning is the same for all units on this electrode
stim = spatialMap.binSpace('window', SPATIAL_MAPPING_WINDOW, 'binSize', SPATIAL_MAPPING_BIN_SIZE, 'correctEyePos', 'simple', 'trialIdx', trialIdx);
spks = spatialMap.binSpikes(stim, spikes, spikes.cids);

baseLambda = 150; % initlize the spatial smoothness parameter

% 2D nonparametric RF with spatial smoothing
[RFcorr, dcCorr, rfCorr_params] = rfmap.spatialRfAutoSmooth(stim.X(stim.valid,:), spks(stim.valid,:), SPATIAL_MAPPING_TEMPORAL_LAGS, stim.size, 'crossvalidation', true, 'lambda0', baseLambda);
s = var(stim.X(stim.valid,:));



figure(2); clf
plot(RFcorr)
% pause

% fold r2 is the variance explained on test data during fitting
goodIx = sum(bsxfun(@gt, abs(RFcorr), s(:))) > 5;
% r2_test = mean(rfCorr_params.fold_r2);

% good RFs
% goodIx = r2_test > 0.015;
% fit 2D gabors and gaussians to the corrected RFs
TwoDRF = rfmap.fitParametric2dRf(stim, RFcorr(:,goodIx));

Coarse.ROI = SPATIAL_MAPPING_WINDOW;
Coarse.binSize = SPATIAL_MAPPING_BIN_SIZE;
Coarse.RFcorr = RFcorr;
Coarse.dcCorr = dcCorr;
Coarse.rfparams = rfCorr_params;
Coarse.goodIx = goodIx;
Coarse.TwoDRF = TwoDRF;
Coarse.stim   = rmfield(stim, {'X', 'valid', 'bins'});

if nargout == 1
    return
end

figure(1); clf
unitCounter = 0;
nUnits = size(spks,2);

% --- Store output of analyses to the neuron struct array
for kUnit = 1:nUnits
    thisUnit  = unitCounter + kUnit;
    sn = 'spatialSquares';
    
    %         if r2_test(kUnit) > 0.01
    %             xy = [TwoDRF(kUnit).gaussian.means.x0, TwoDRF(kUnit).gaussian.means.y0]; % mean
    %             C = [TwoDRF(kUnit).gaussian.means.sigmax 0; 0 TwoDRF(kUnit).gaussian.means.sigmay]; % covariance
    %         else
    xy = nan;
    C  = nan;
    %         end
    
    neuron(thisUnit).(sn) = struct('xax', stim.xax, ...
        'yax', stim.yax, ...
        'RF', reshape(RFcorr(:,kUnit),stim.size), ...
        'center', xy, ...
        'covariance', C, ...
        'dc', dcCorr(kUnit));
    
    figure(1);
    subplot(ceil(sqrt(nUnits)), round(sqrt(nUnits)), kUnit)
    imagesc(stim.xax, stim.yax, neuron(thisUnit).(sn).RF);
    hold on
    parRF = (find(goodIx)==thisUnit);
    if any(parRF)
        xy = [TwoDRF(parRF).gaussian.means.x0, TwoDRF(parRF).gaussian.means.y0];
        plot.plotellipse(xy, [TwoDRF(parRF).gaussian.means.sigmax 0; 0 TwoDRF(parRF).gaussian.means.sigmay], 1, 'Linewidth', 2, 'Color', 'r');
    end
    %
    axis xy
    drawnow
    
end % end loop over units to store RF mapping info

%%
figure(2); clf
for i =1:numel(TwoDRF)
    xy = [TwoDRF(i).gaussian.means.x0, TwoDRF(i).gaussian.means.y0];
    plot.plotellipse(xy, [TwoDRF(i).gaussian.means.sigmax 0; 0 TwoDRF(i).gaussian.means.sigmay], 1, 'Linewidth', 2); hold on
end

xd = xlim;
yd = ylim;

%% ---- FINE MAPPING

% trial = pds.getPdsTrialData(PDS);
% ROI = round([xd(1) yd(1) xd(2) yd(2)]*trial(1).display.ppd);
% ROI = [xd(1) yd(1) xd(2) yd(2)];
ROI = [-2 -2 2 2];
SPATIAL_MAPPING_WINDOW        = ROI;

SPATIAL_MAPPING_BIN_SIZE      = min(squareSizes);

trialIdx = 1:numel(squareSizesTrial)-1; %find(squareSizesTrial > 0); %find(squareSizesTrial==SPATIAL_MAPPING_BIN_SIZE);

SPATIAL_MAPPING_BIN_SIZE      = SPATIAL_MAPPING_BIN_SIZE  * 2;

% binning is the same for all units on this electrode
stim = spatialMap.binSpace('window', SPATIAL_MAPPING_WINDOW, 'binSize', SPATIAL_MAPPING_BIN_SIZE, 'correctEyePos', 'simple', 'trialIdx', trialIdx);
spks = spatialMap.binSpikes(stim, spikes, spikes.cids);

%%
baseLambda = 5e3;
% 2D nonparametric RF with spatial smoothing
[RFcorr, dcCorr, rfCorr_params] = rfmap.spatialRfAutoSmooth(stim.X(stim.valid,:), spks(stim.valid,:), SPATIAL_MAPPING_TEMPORAL_LAGS, stim.size, 'crossvalidation', true, 'lambda0', baseLambda);
% XX = stim.X'*stim.X;
%%
% Cinv = 5*eye(size(XX,1)) + qfsmooth2nd(stim.size(1), stim.size(2));
% RFcorr = (XX + Cinv*3e3)\(stim.X'*spks);
s = var(stim.X);
figure(1); clf
plot(RFcorr)
hold on
plot(s, 'k')
%
% fold r2 is the variance explained on test data during fitting
r2_test = mean(rfCorr_params.fold_r2);

% goodIx = bsxfun(@gt, abs(RFcorr), s(:));
% RFcorr(~goodIx) = 0;

% good RFs
% goodIx = r2_test > 0.015;
% fit 2D gabors and gaussians to the corrected RFs
% TwoDRF = rfmap.fitParametric2dRf(stim, RFcorr(:,goodIx));

F.ROI = ROI;
F.binSize = SPATIAL_MAPPING_BIN_SIZE;
F.RFcorr = RFcorr;
F.dcCorr = dcCorr;
F.rfparams = rfCorr_params;
% F.goodIx = goodIx;
% F.TwoDRF = TwoDRF;
F.stim   = rmfield(stim, {'X', 'valid', 'bins'});


figure(3); clf
unitCounter = 0;
nUnits = size(spks,2);

% --- Store output of analyses to the neuron struct array
for kUnit = 1:nUnits
    thisUnit  = unitCounter + kUnit;
    sn = 'spatialSquares';
    
    %         if r2_test(kUnit) > 0.01
    %             xy = [TwoDRF(kUnit).gaussian.means.x0, TwoDRF(kUnit).gaussian.means.y0]; % mean
    %             C = [TwoDRF(kUnit).gaussian.means.sigmax 0; 0 TwoDRF(kUnit).gaussian.means.sigmay]; % covariance
    %         else
    xy = nan;
    C  = nan;
    %         end
    
    neuron(thisUnit).(sn) = struct('xax', stim.xax, ...
        'yax', stim.yax, ...
        'RF', reshape(RFcorr(:,kUnit),stim.size), ...
        'center', xy, ...
        'covariance', C, ...
        'r2', r2_test(kUnit), ...
        'dc', dcCorr(kUnit));
    
    figure(1);
    subplot(ceil(sqrt(nUnits)), round(sqrt(nUnits)), kUnit)
    imagesc(stim.xax, stim.yax, neuron(thisUnit).(sn).RF);
    hold on
%     parRF = (find(goodIx)==thisUnit);
%     if any(parRF)
%         xy = [TwoDRF(parRF).gaussian.means.x0, TwoDRF(parRF).gaussian.means.y0];
%         plot.plotellipse(xy, [TwoDRF(parRF).gaussian.means.sigmax 0; 0 TwoDRF(parRF).gaussian.means.sigmay], 1, 'Linewidth', 2, 'Color', 'r');
%     end
    %
    axis xy
    drawnow
    
end % e
%     natImgIdx = find(~arrayfun(@(x) isempty(x.natImgBackground), trial) | arrayfun(@(x) ~x.natImgBackground.use, trial));
%     natImgIdx = natImgIdx(arrayfun(@(x) x.natImgBackground.use==1, trial(natImgIdx)));
    %%
%     I = pds.replayTrial(trial, natImgIdx, 'ROI', ROI, 'gazeContingent', true, 'showOverlay', false);
    
    %%
    