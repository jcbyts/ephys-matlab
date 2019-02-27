function RF = CoarseRFMap(PDS, sp)
% checks if spatial RF mapping was run, uses a coarse initialization of the
% RF to measure
%%

spatialMap = session.squareFlash(PDS);

if spatialMap.numTrials > 10
    
    SPATIAL_MAPPING_WINDOW        = [-15 -10 15 10];
    SPATIAL_MAPPING_TEMPORAL_LAGS = 1:8;
    
    squareSizesTrial = arrayfun(@(x) x.size, spatialMap.trial);
    squareSizes = unique(squareSizesTrial);
    
    SPATIAL_MAPPING_BIN_SIZE      = max(squareSizes);
    
    trialIdx = find(squareSizesTrial==SPATIAL_MAPPING_BIN_SIZE);
    
    % double the bin size
    SPATIAL_MAPPING_BIN_SIZE      = SPATIAL_MAPPING_BIN_SIZE  * 2;
    
    % 2D nonparametric RF with spatial smoothing
    baseLambda = 10; % initlize the spatial smoothness parameter
    
    % binning is the same for all units on this electrode
    stim = spatialMap.binSpace('window', SPATIAL_MAPPING_WINDOW, 'binSize', SPATIAL_MAPPING_BIN_SIZE, 'correctEyePos', 'simple', 'trialIdx', trialIdx);
    spks = spatialMap.binSpikes(stim, sp, sp.cids);
    
    [RFcorr, dcCorr, rfCorr_params] = rfmap.spatialRfAutoSmooth(stim.X(stim.valid,:), spks(stim.valid,:), SPATIAL_MAPPING_TEMPORAL_LAGS, stim.size, 'crossvalidation', true, 'lambda0', baseLambda);
    
    r2_test = mean(rfCorr_params.fold_r2,1);
    
    % good RFs
    goodIx = r2_test > 0.015;
    
    % fit 2D gabors and gaussians to the corrected RFs
    TwoDRF = rfmap.fitParametric2dRf(stim, RFcorr(:,goodIx));
    
    figure(1); clf;
    nUnits = size(spks,2);
    for kUnit = 1:nUnits
        sn = 'spatialSquares';
        
        parRF = (find(goodIx)==kUnit);
        if any(parRF)
            xy = [TwoDRF(parRF).gaussian.means.x0, TwoDRF(parRF).gaussian.means.y0]; % mean
            C = [TwoDRF(parRF).gaussian.means.sigmax 0; 0 TwoDRF(parRF).gaussian.means.sigmay]; % covariance
        else
            xy = nan;
            C  = nan;
        end
        
        neuron(kUnit).(sn) = struct('xax', stim.xax, ...
            'yax', stim.yax, ...
            'RF', reshape(RFcorr(:,kUnit),stim.size), ...
            'center', xy, ...
            'covariance', C, ...
            'r2', r2_test(kUnit), ...
            'dc', dcCorr(kUnit));
        
        
        subplot(ceil(sqrt(nUnits)), round(sqrt(nUnits)), kUnit)
        imagesc(stim.xax, stim.yax, neuron(kUnit).(sn).RF);
        hold on
        
        if any(parRF)
            xy = [TwoDRF(parRF).gaussian.means.x0, TwoDRF(parRF).gaussian.means.y0];
            plot.plotellipse(xy, [TwoDRF(parRF).gaussian.means.sigmax 0; 0 TwoDRF(parRF).gaussian.means.sigmay], 1, 'Linewidth', 2, 'Color', 'r');
        end
        
        axis xy
        drawnow
        
    end
end

RF = cell2mat(arrayfun(@(x) x.(sn), neuron, 'uni', false));
    
return