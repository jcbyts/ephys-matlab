%% refresh session
% rerun getExperimentsAnd to refresh the meta data
meta = io.getExperimentsAnd(); % get all experiments meta data

%thisSession = meta(148,:); % again grab the session you are working with
thisSession = meta(end-3,:);
disp(thisSession)

PDS = io.getPds(thisSession);

sp = io.getSpikes(thisSession, 'Kilo');

%% plot spike waveforms
ops = io.loadOps(thisSession);

figure(1); clf
nShanks = numel(ops);

for iShank = 1:nShanks
    figure(iShank); clf
    if isfield(sp{iShank}, 'uQ')
        [~, ind] = sort(sp{iShank}.clusterDepths);
        clusterIds = sp{iShank}.cids(ind);
    else
        clusterIds = [];
    end
    fig(iShank) = plot.spikeWaveformsFromOps(ops(iShank), sp{iShank}, 'figure', iShank, 'clusterIds', clusterIds, 'numWaveforms', 100);
end



%% detect saccades
% This is quick and dirty. We should probably replace this with something
% more robust. I doubt we want to score every saccade in the GUI, but we
% might want to score enough of them to know that we trust the algorithm.

% load Eye position
[data, timestamps, elInfo] = io.getEdf(thisSession);
% data(2,:) = -data(2,:);
% --- remove bad samples
ix = any(data(1:2,:) == elInfo.bitDeg(2));
data(1,ix) = nan;
data(2,ix) = nan;

%% --- detect those saccades
ix = ~any(isnan(data));
[saccades] = pdsa.detectSaccades(timestamps(ix), data(1:2,ix), ...
    'verbose', false, ...
    'filterPosition', 1, ...
    'filterLength', ceil(20/elInfo.sampleRate*1e3), ... % 40 ms smoothing for velocity computation
    'detectThresh', 200, ...
    'startThresh', 5, ...
    'minIsi', ceil(50/elInfo.sampleRate*1e3), ...
    'minDur', ceil(4/elInfo.sampleRate*1e3), ... % 4 ms
    'blinkIsi', ceil(40/elInfo.sampleRate*1e3));

%% --- plot output
figure(1); clf
plot(timestamps, data(1,:), 'k'); hold on
plot(timestamps, data(2,:), 'Color', repmat(.5, 1, 3));


tbins = timestamps(1:3e3:end);
cnt = histc(saccades.start', tbins);
[~, id] = max(cnt);

plot([saccades.start'; saccades.start'], ylim, 'r--')
xlim(tbins(id + [0 1]))
xlabel('seconds')
ylabel('degrees')
legend({'x', 'y', 'saccade start'})

figure(2); clf
plot(saccades.duration/elInfo.sampleRate, saccades.size, '.')
xlabel('duration (s)')
ylabel('amplitude (deg)')

figure(3); clf
dx = saccades.endXpos - saccades.startXpos;
dy = saccades.endYpos - saccades.startYpos;

[cnt, bins] = hist3([dx(:) dy(:)], [40 40]);

imagesc(bins{1}, bins{2}, cnt)
axis xy


[cnt, bins] = hist3([data(2,:)', data(1,:)'], [400 400]);

figure(3); clf
imagesc(bins{1}, bins{2}, log(cnt))
axis xy

hold on 
% plot(data(1,:), data(2,:), '.r')
%% build stimulus object
spatialMap = session.squareFlash(PDS, 'eyetrace', [timestamps'; data], 'saccades', saccades);

%% bin space and create the design matrix
%x,y,x,y
window  = [0 -10 20 0]; %[-3, 3]% degrees (window is the same in x and y -- TODO: separate)
binSize = 1; % degrees

stim = spatialMap.binSpace('window', window, 'binSize', binSize, 'correctEyePos', 'simple');
stimBase = spatialMap.binSpace('window', window, 'binSize', binSize, 'correctEyePos', 'no');

%% Bin up spike counts
T = size(stim.X,1);
spks = zeros(T,1);
clusterIds = sp{1}.cids(sp{1}.cgs >=2);
units = clusterIds;
% [~, inds] = sort(sp{1}.clusterDepths);
% units = sp{1}.cids(inds);
nUnits = numel(units);
for i = 1:nUnits
    cnt = histc(sp{1}.st(sp{1}.clu == units(i)), stim.bins);
    cnt((diff(stim.bins) > (spatialMap.display.ifi*1.5))) = 0;
    spks(:,i) = cnt;
end

%% Quick check to make sure the mapping is working

% fit spatial RFs
temporal_lags = 1:5;
baseLambda = 10;
[RFcorr, dcCorr, rfCorr_params] = rfmap.spatialRfAutoSmooth(stim.X(stim.valid,:), spks(stim.valid,:), temporal_lags, stim.size, 'crossvalidation', false, 'lambda0', baseLambda);

[RFbase, dcBase, rfBase_params] = rfmap.spatialRfAutoSmooth(stimBase.X(stimBase.valid,:), spks(stimBase.valid,:), temporal_lags, stimBase.size, 'crossvalidation', false, 'lambda0', baseLambda);

I = RFcorr;
figure(1); clf
sx = ceil(sqrt(nUnits));
sy = round(sqrt(nUnits));
ax = pdsa.tight_subplot(sx, sy, .05, .01);

for i = 1:(sx*sy)
    set(gcf, 'currentaxes', ax(i))
    if i <= nUnits
        imagesc(stim.xax, stim.yax, reshape(I(:,i), stim.size)); colormap gray
        set(ax(i), 'gridcolor', 'y')
        grid on
        axis xy
    else
        set(ax(i), 'Visible', 'off')
    end
end
%% Fit with and without eye correction with cross validation for hyper parameter
baseLambda = 10;
[RFcorr, dcCorr, rfCorr_params] = rfmap.spatialRfAutoSmooth(stim.X(stim.valid,:), spks(stim.valid,:), temporal_lags, stim.size, 'crossvalidation', true, 'lambda0', baseLambda);
[RFbase, dcBase, rfBase_params] = rfmap.spatialRfAutoSmooth(stimBase.X(stimBase.valid,:), spks(stimBase.valid,:), temporal_lags, stimBase.size, 'crossvalidation', true, 'lambda0', baseLambda*1000);

% fit 2D gabors and gaussians to the corrected RFs
results = rfmap.fitParametric2dRf(stim, RFcorr);

%% plot the RFs
fits = {RFcorr, RFbase};
for iFit = 1:numel(fits)
    I = fits{iFit};
    figure(iFit); clf
    sx = ceil(sqrt(nUnits));
    sy = round(sqrt(nUnits));
    ax = pdsa.tight_subplot(sx, sy, .05, .01);
    if isfield(rfCorr_params, 'fold_r2')
        r2_test = mean(rfCorr_params.fold_r2);
    else
        r2_test = ones(1, nUnits);
    end
    
    for i = 1:(sx*sy)
        set(gcf, 'currentaxes', ax(i))
        if i <= nUnits && r2_test(i) > .01
            imagesc(stim.xax, stim.yax, reshape(I(:,i), stim.size)); colormap gray
            set(ax(i), 'gridcolor', 'y')
            grid on
            axis xy
            
            hold on
            xy = [results(i).gaussian.means.x0, results(i).gaussian.means.y0];
            plotellipse(xy, [results(i).gaussian.means.sigmax 0; 0 results(i).gaussian.means.sigmay], 1, 'r')
        else
            set(ax(i), 'Visible', 'off')
        end
    end
end
    
%% Flash / Saccade triggered CSD
ops = io.loadOps(thisSession);
% pick a shank and load LFP
iShank = 1;
[lfp, lfpTime, lfpInfo] = io.getLFP(ops(iShank));
lfpInfo.fragments = round(lfpInfo.fragments);
% --- flash-triggered
csdTrial = io.csdTrial(PDS);
if isfield(csdTrial, 'onset')
    flashTimes = [csdTrial.onset];
else
    flashTimes = 1;
end

if ops(iShank).Nchan > 10
    figure(1); clf
    subplot(1,2,1)
    plot.CsdBasic(lfp, flashTimes, lfpInfo)
    title('Flash Triggered')
    xlabel('ms')
    ylabel('depth')   
end

%% --- saccade-triggered
eventTimes = saccades.start';
% find saccades that happened more than 200 ms after the last saccade
% (uncontaminated)
ixLong = find([0 diff(eventTimes)>.2]);

if ops(iShank).Nchan > 10
    subplot(1,2,2)
    plot.CsdBasic(lfp, eventTimes(ixLong), lfpInfo)
    title('Saccade Triggered')
    xlabel('ms')
    ylabel('depth')
end

% Look for saccade-direction tuning in the LFP
th = cart2pol(saccades.endXpos - saccades.startXpos, saccades.endYpos - saccades.startYpos)'/pi*180;

if ops(iShank).Nchan > 10
    figure(2); clf
end

et = eventTimes(ixLong)';
ev = io.convertTimeToSamples(et, lfpInfo.sampleRate, lfpInfo.timestamps(:), ceil(lfpInfo.fragments(:)));
th = th(ixLong(:));
thBins = 0:45:360;
ch0 = (1:32)*40;
cmap = flipud(hsv(numel(thBins)));
if ops(iShank).Nchan > 10
    subplot(1,2,2)
    for i = 1:numel(thBins)-1
        
        ii = th > thBins(i) & th < thBins(i+1);
        [sta,~,bc] = pdsa.eventTriggeredAverage(lfp, ev(ii), [-300 400]);
        
        plot(bc, bsxfun(@plus, sta, ch0), 'Color', cmap(i,:)); hold on
        %     plot(bc, sta+sd, '--', 'Color', cmap(i,:));
        %     plot(bc, sta-sd, '--', 'Color', cmap(i,:));
        
    end
    
    xlabel('Time from saccade')
    
    axis tight
    xlim([-100 200])
    yd = ylim;
    ax2 = axes('Position', get(gca, 'Position'));
    for i = 1:numel(thBins)-1
        quiver(ax2, 0, 0, cosd(mean(thBins(i + [0 1]))), sind(mean(thBins(i + [0 1]))), 'Color', cmap(i,:), 'AutoScale', 'off'); hold on
        axis off
        
    end
    xlim([-1 10])
    ylim([-10 1])
    title('Saccade-triggered LFP')
    
    subplot(1,2,1)
    % flash triggered LFP overlayed
    flashTimes = [csdTrial.onset]';
    flashSamples = io.convertTimeToSamples(flashTimes, lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));
    [sta,~,bc] = pdsa.eventTriggeredAverage(lfp, flashSamples, [-300 400]);
    plot(bc, bsxfun(@plus, sta, ch0), 'Color', 'k');
    axis tight
    ylim(yd)
    xlim([-100 200])
    xlabel('Time from flash')
    title('Flash-triggered LFP')
    
else
    clf
    for i = 1:numel(thBins)-1
        
        ii = th > thBins(i) & th < thBins(i+1);
        [sta,~,bc] = pdsa.eventTriggeredAverage(lfp, ev(ii), [-300 400]);
        
        plot(bc, sta, 'Color', cmap(i,:)); hold on
        %     plot(bc, sta+sd, '--', 'Color', cmap(i,:));
        %     plot(bc, sta-sd, '--', 'Color', cmap(i,:));
        
    end
    
    xlabel('Time from saccade')
    
    axis tight
    xlim([-100 200])
    yd = ylim;
    ax2 = axes('Position', get(gca, 'Position'));
    for i = 1:numel(thBins)-1
        quiver(ax2, 0, 0, cosd(mean(thBins(i + [0 1]))), sind(mean(thBins(i + [0 1]))), 'Color', cmap(i,:), 'AutoScale', 'off'); hold on
        axis off
        
    end
    xlim([-1 10])
    ylim([-10 1])
    title('Saccade-triggered LFP')
end

%% ********** spiking locked to saccades **********
iShank = 1;
figure(3); clf
s = sp{iShank};
if isfield(s, 'uQ')
    goodUnits = s.cgs > 1;
    udepths = s.clusterDepths(goodUnits);
    uId = s.cids(goodUnits);
    [~, depthId] = sort(udepths);
    clustId = uId(depthId);
else
    clustId = s.cids;
end
nUnits = numel(clustId); % including

sx = ceil(sqrt(nUnits));
sy = round(sqrt(nUnits));

ax = pdsa.tight_subplot(sy, sx, 0.01,  0.01);
for kUnit = 1:nUnits
    st = s.st(s.clu == clustId(kUnit));
    
    set(gcf, 'currentaxes', ax(kUnit))
    for i = 1:numel(thBins)-1
        
        ii = th > thBins(i) & th < thBins(i+1);
        [sta, sd, bc] = pdsa.eventPsth(st, et(ii), [-.3 .3], .01);
        
        plot(bc, sta, 'Color', cmap(i,:)); hold on
        plot(bc, sta+sd, '--', 'Color', cmap(i,:));
        plot(bc, sta-sd, '--', 'Color', cmap(i,:));
        
    end
    plot([0 0], ylim, 'k')
    axis tight
    %axis off
    
end

%% Hartley
hart = session.hartleyFF(PDS);

% build design matrix for STA
hart.buildDesignMatrix();


%% 
iShank = 1;
s = sp{iShank};
if isfield(s, 'cgs')
    clustIds = s.cids(s.cgs>1);
else
    clustIds = s.cids;
end
nUnits = numel(clustIds);

figure(2); clf
ax = pdsa.tight_subplot(nUnits, 2, .01, .01);
for kUnit = 1:nUnits
    spikeTimes = s.st(s.clu==clustIds(kUnit));
%     set(gcf, 'currentaxes', ax((kUnit-1)*2 + 1))
    pdsa.plotRaster(spikeTimes, [hart.trial.start], [-.1 20], .01)
    pause
end

    sta = hart.spikeTriggeredAverage(spikeTimes);
    
%     xdat.xx = hart.design.XX;
%     xdat.xy = sta.xy;
%     xdat.yy = sta.yy;
%     xdat.ny = sta.ny;
% %     autoCorrRidgeRegress(xdat)
% %     sta = hart.AsdRf(spikeTimes);
%     [w_hat, wt, wx, wlin] = bilinearMixRegress_coordAscent(xdat.xx + 10e2*speye(size(xdat.xx,2)), xdat.xy, [hart.design.nkTime hart.design.nkx*hart.design.nky], 1, 1:size(xdat.xx,2)-1);
%     sflip = sign(sum(wt));
%     sta.RF = reshape(sflip*wx, [hart.design.nkx hart.design.nky]);
%     sta.RFtime = sflip*wt;
    set(gcf, 'currentaxes', ax((kUnit-1)*2 + 1))
    imagesc(sta.kxs, sta.kys, sta.RF)
    colormap(gray.^2)
    axis xy
    axis off
    set(gcf, 'currentaxes', ax((kUnit-1)*2 + 2))
    plot(sta.time, sta.RFtime)
    axis tight
    axis off
    drawnow

end
%% plot frozen trials
figure(11); clf
hart.plotFrozenRaster(s);

%% dot revco mapping
try
    mtmap = session.mtDotRcMap(PDS, true);
    mtmap.buildDesignMatrix('xwin', [-8 8], 'ywin', [-8 8], 'binSize', 1);
catch
    error('No MT dot mapping trials')
end
%%
% figure(3); clf
% ax = pdsa.tight_subplot(nUnits, 2, .01, .01);

iShank = 1;
s = sp{iShank};
Xsmooth = filter(ones(5,1)/5, 1, mtmap.design.X); % filter forward
C= (Xsmooth'*Xsmooth);
C = C + 10e6*eye(size(Xsmooth,2));
% plot.spikeWaveformsFromOps(ops(2), s)
for kUnit = 1:nUnits
    spikeTimes = s.st; %(s.clu==clustIds(kUnit));
    
    y = histc(spikeTimes, mtmap.design.rowTimes);
    
    y(diff(mtmap.design.rowTimes) > 1.5*mtmap.display.ifi) = 0;
    y(y>15) = 0;
    figure(1); clf
    subplot(2,2,1:2)
    plot(smooth(y,10))
    
    
    
    subplot(2,2,3)
    sta = C\Xsmooth'*(y-mean(y));
    imagesc(mtmap.design.xax, mtmap.design.yax, reshape(sta, mtmap.design.sz(1), mtmap.design.sz(2)*2))
    
    pause
%     yy = y'*y';
%     ny = numel(y);
%     
%     [w_hat, wt, wx, wlin] = bilinearMixRegress_coordAscent(mtmap.design.XX + 10e2*speye(size(mtmap.design.XX,2)), xy, [mtmap.design.nkTime prod(mtmap.design.sz)], 1, 1:size(mtmap.design.XX,2)-1);
%     sflip = sign(sum(wt));
%     sta.RF = reshape(sflip*wx, mtmap.design.sz);
%     sta.RFtime = sflip*wt;
%     set(gcf, 'currentaxes', ax((kUnit-1)*2 + 1))
%     imagesc(mtmap.design.xax, mtmap.design.yax, sta.RF)
%     colormap(gray.^2)
%     colormap jet
%     axis xy
% %     axis off
%     set(gcf, 'currentaxes', ax((kUnit-1)*2 + 2))
%     plot(sta.RFtime)
%     axis tight
%     axis off
%     drawnow
%     
    
    
end

%% try bilinear RF estimation
figure(3); clf
ax = pdsa.tight_subplot(nUnits, 2, .01, .01);

plot.spikeWaveformsFromOps(ops(2), s)
for kUnit = 1:nUnits
    spikeTimes = s.st; %(s.clu==clustIds(kUnit));
    
    y = histc(spikeTimes, mtmap.design.rowTimes);
    
    y(diff(mtmap.design.rowTimes) > 1.5*mtmap.display.ifi) = 0;
    
    figure(4)
    plot(smooth(y,10))
    pause
    
    figure(3)
    xy = mtmap.design.Xd'*y;
%     yy = y'*y';
    ny = numel(y);
    
    [w_hat, wt, wx, wlin] = bilinearMixRegress_coordAscent(mtmap.design.XX + 10e2*speye(size(mtmap.design.XX,2)), xy, [mtmap.design.nkTime prod(mtmap.design.sz)], 1, 1:size(mtmap.design.XX,2)-1);
    sflip = sign(sum(wt));
    sta.RF = reshape(sflip*wx, mtmap.design.sz);
    sta.RFtime = sflip*wt;
    set(gcf, 'currentaxes', ax((kUnit-1)*2 + 1))
    imagesc(mtmap.design.xax, mtmap.design.yax, sta.RF)
    colormap(gray.^2)
    colormap jet
    axis xy
%     axis off
    set(gcf, 'currentaxes', ax((kUnit-1)*2 + 2))
    plot(sta.RFtime)
    axis tight
    axis off
    drawnow
    
    
    
end



%% gaussian pyramid noise

stim = 'gaussianNoiseBlobs';

ppd  = PDS{1}.initialParametersMerged.display.ppd;
ifi  = PDS{1}.initialParametersMerged.display.ifi;

hasStim = io.findPDScontainingStimModule(PDS, stim);
hasFixation = cellfun(@(x) strncmp(x.initialParametersMerged.pldaps.trialFunction, 'stimuli.fixflash', 16), PDS);

hasStim = hasStim(:); % & ~hasFixation(:);

gpyrTrial = struct();
trialNum = 0;

for i = find(hasStim(:)')
    
    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 
    
    stimTrials = find(trialIx);
    
    if isempty(stimTrials)
        continue
    end
    
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
    
        kTrial = trialNum + j;
        
        % --- Timing
        gpyrTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
        gpyrTrial(kTrial).start      = gpyrTrial(kTrial).frameTimes(1);
        gpyrTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - gpyrTrial(kTrial).start;
        
        % --- Eye position
        eyepos = io.getEyePosition(PDS{i}, thisTrial);
        gpyrTrial(kTrial).eyeSampleTime = eyepos(:,1);
        gpyrTrial(kTrial).eyeXPx        = eyepos(:,2);
        gpyrTrial(kTrial).eyeYPx        = eyepos(:,3);
        gpyrTrial(kTrial).pupilArea     = eyepos(:,4);
        
        results = pdsa.detectSaccades(eyepos(:,1)', eyepos(:,2:3)'./ppd, 'verbose', false);
        
        
        
        iix = gpyrTrial(kTrial).eyeXPx < 200 | gpyrTrial(kTrial).eyeXPx > 1800;
        iiy = gpyrTrial(kTrial).eyeYPx < 100 | gpyrTrial(kTrial).eyeXPx > 1000;
        bad = iix | iiy;
       
        gpyrTrial(kTrial).eyeXPx(bad) = nan;
        gpyrTrial(kTrial).eyeYPx(bad) = nan;
        gpyrTrial(kTrial).pupilArea(bad) = nan;
        
        % find eye position on each frame
        t = gpyrTrial(kTrial).frameTimes - gpyrTrial(kTrial).start;
        tdiff = abs(bsxfun(@minus, gpyrTrial(kTrial).eyeSampleTime, t(:)')) < mean(diff(gpyrTrial(kTrial).eyeSampleTime));
        [irow,icol] = find(diff(tdiff)==-1);
        
        gpyrTrial(kTrial).saccades = false(size(gpyrTrial(kTrial).frameTimes));
        
        for iSaccade = 1:size(results,2)
            ix = t > results(1, iSaccade) & t < results(2,iSaccade);
            gpyrTrial(kTrial).saccades(ix) = true;
        end
        
        gpyrTrial(kTrial).eyePosAtFrame = [gpyrTrial(kTrial).eyeXPx(irow) gpyrTrial(kTrial).eyeYPx(irow)];
        
        gpyrTrial(kTrial).rngs       = PDS{i}.data{thisTrial}.(stim).rngs;
        gpyrTrial(kTrial).n          = PDS{i}.data{thisTrial}.(stim).n;
        gpyrTrial(kTrial).xpos       = PDS{i}.data{thisTrial}.(stim).xpos;
        gpyrTrial(kTrial).ypos       = PDS{i}.data{thisTrial}.(stim).ypos;
        gpyrTrial(kTrial).gridpos    = PDS{i}.data{thisTrial}.(stim).gridpos;
        gpyrTrial(kTrial).contrast   = repmat(PDS{i}.data{thisTrial}.(stim).n.mypars(1,:), numel(gpyrTrial(kTrial).frameTimes), 1);
        
        gpyrTrial(kTrial).xposEye    = bsxfun(@minus, gpyrTrial(kTrial).xpos, gpyrTrial(kTrial).eyePosAtFrame(:,1));
        gpyrTrial(kTrial).yposEye    = -bsxfun(@minus, gpyrTrial(kTrial).ypos, gpyrTrial(kTrial).eyePosAtFrame(:,2));
        
        
    end
    
    trialNum = kTrial;
end


%% concatenate trials

flipTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), gpyrTrial, 'UniformOutput', false)');
xposRel   = cell2mat(arrayfun(@(x) x.xposEye, gpyrTrial, 'UniformOutput', false)');
yposRel   = cell2mat(arrayfun(@(x) x.yposEye, gpyrTrial, 'UniformOutput', false)');
eyeposX   = cell2mat(arrayfun(@(x) x.eyeXPx, gpyrTrial, 'UniformOutput', false)');
eyeposY   = cell2mat(arrayfun(@(x) x.eyeYPx, gpyrTrial, 'UniformOutput', false)');
contrast  = cell2mat(arrayfun(@(x) x.contrast, gpyrTrial, 'UniformOutput', false)');
saccades  = cell2mat(arrayfun(@(x) x.saccades(:), gpyrTrial, 'UniformOutput', false)');
validTimes = ~(any(isnan(xposRel),2) | any(isnan(yposRel),2)) & ~saccades;
t0 = flipTimes(1);
tEnd = flipTimes(end);



flipTimes = flipTimes(validTimes);
xposRel   = xposRel(validTimes, :);
yposRel   = yposRel(validTimes, :);
contrast  = contrast(validTimes,:);

[flipTimes, id]=sort(flipTimes);

xposRel = xposRel(id,:);
yposRel = yposRel(id,:);
contrast = contrast(id,:);


%% Build design matrix

sc = gpyrTrial(end).n.scale*gpyrTrial(end).n.sc;
lvls = sort(unique(sc), 'descend');

xwin  = [-5 5]*ppd;
ywin  = [-5 5]*ppd;

binSize = 1; %lvls(kLevel)*37.4400;
xax = xwin(1):binSize:xwin(2);
yax = ywin(1):binSize:ywin(2);

[xx,yy]=meshgrid(xax, yax);

xax = xax/ppd;
yax = yax/ppd;

sz = ceil([diff(ywin)/binSize diff(xwin)/binSize]);

binfun = @(x) (x==0) + ceil(x/binSize);

Xsmooth = [];

gkern = {};
levels = 1;
for kLevel = 1:numel(levels)
    thisLevel = lvls(levels(kLevel));
    levelIx = sc == thisLevel;
    xpos = xposRel(:,levelIx);
    ypos = yposRel(:,levelIx);
    
    c    = contrast(:,levelIx);

% ixGood = (xpos > xwin(1)) & (xpos < xwin(2)-1) & (ypos > ywin(1)) & (ypos < ywin(2)-1);



gkern{kLevel} = exp(-((xx-mean(xwin)).^2 + (yy-mean(ywin)).^2)/(2*thisLevel.^2));

xgrid = binfun(xpos - xwin(1));
ygrid = binfun(ypos - ywin(1));

xgrid(xgrid<1) = nan;
ygrid(ygrid<1) = nan;

xgrid(xgrid>sz(2)) = nan;
ygrid(ygrid>sz(1)) = nan;

gridpos = sub2ind(sz, ygrid, xgrid);

% gridpos(~ixGood)=nan;

ix = ~isnan(gridpos);

[frameNumber,gaussianNumber] = find(ix);

Xsmooth = [Xsmooth sparse(frameNumber, gridpos(ix), c(ix), numel(flipTimes), prod(sz))];

end

addpath l1_ls_matlab\


%% sanity check
figure(1); clf
kFrame = randi(numel(flipTimes));
tmp = sum(reshape(Xsmooth(kFrame,:), [], numel(levels)),2);
imagesc(xax*ppd, yax*ppd, reshape(tmp, sz)); hold on


colormap gray
plot(xposRel(kFrame,levelIx), yposRel(kFrame,levelIx), 'or')
%% measure STA
addBias = 0;
if addBias
    Xd = [Xsmooth ones(size(Xsmooth,1), 1)]; %X.^2;
else
    Xd = Xsmooth;
end
% C = Xd'*Xd;

figure(1); clf
nf = 1;


inds = unitIds([4 6]);
nn = numel(inds);
ax = pdsa.tight_subplot(nn, nf, .11, .1);

for kUnit = 1:nn
%
st = sp.st(sp.clu==(inds(kUnit)));
%
y = histc(st, flipTimes);
bad = diff(flipTimes)> (2*ifi);
y(bad) = 0;



% C = Xd'*Xd;

%


for ts = 1:nf
    set(gcf, 'currentaxes', ax((kUnit-1)*nf + ts))
%     subplot(nUnits,nf,(kUnit-1)*nf + ts)
    ysamp = y;
    if nf==1
        ysamp = smooth(ysamp, 3);
        ysamp(1:end-(7-1)) = ysamp(7:end);
    else
        ysamp(1:end-(ts-1)) = ysamp(ts:end);
    end
    
   
    ysamp = ysamp - mean(ysamp);
%     ysamp = ysamp ./ std(ysamp);
    w =(Xd'*ysamp);
% %     
% %     % ridge regression
% % %     if binS
%     Cprior = speye(size(Xd,2));
%     Cprior = toeplitz([-1 1 zeros(1,size(Xd,2)-3)];
%     
%     w = (C + 10e2*Cprior)\w;

%     w=l1_ls(Xd, ysamp, 8);

%     w = fastASD(X, ysamp, sz, .1);
    if addBias
        w(end)=[];
    end
    I = zeros(sz);
    for kLevel = 1:numel(gkern)
        nb = prod(sz);
        ix = (kLevel-1)*nb + (1:nb);
        I_ = reshape(w(ix), sz);
        I_ = conv2(I_, gkern{kLevel}, 'same');
        I = I + I_;
    end
    I = I - mean(I(:));
    I = sign(I) .* I.^2;
    
    imagesc(xax, yax, I); %, [-1 1]*90);
    axis xy
    colormap(gray)
    drawnow
    grid on
    set(gca, 'XTick', [-5:1:5], 'YTick', [-5:1:5], 'GridColor', 'y')
%     imagesc(xax, yax, I);
%     title(sprintf('frame -%d', ts))
end

end

%% 


y = histc(spikeTimes, stim.bins);
y(diff(stim.bins) > spatialMap.display.ifi*2) = 0;
y = sqrt(y);
nlags = 20;
xc = zeros(nlags*2+1,size(stim.X,2)); 
for i = 1:size(stim.X,2)
xc(:,i) = xcorr(y,stim.X(:,i), nlags);
end


%%
S = prf.fit_generic_linear_model(stim.X, y, 20, stim.size);


%%

imagesc(reshape(Xd(:,1:end-1)'*y, [], prod(stim.size)))

% clf
% imagesc(reshape(sum(reshape(Xd(:,1:end-1)'*y, [], prod(stim.size))), stim.size))

%%