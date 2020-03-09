%% refresh session
% rerun getExperimentsAnd to refresh the meta data
meta = io.getExperimentsAnd(); % get all experiments meta data

%thisSession = meta(148,:); % again grab the session you are working with
thisSession = meta(end-2,:);
disp(thisSession)

fprintf('Loading Behavioral data from the server...')
tic
PDS = io.getPds(thisSession);
fprintf(' [%02.2f]\n', toc)

fprintf('Loading Spiking data from the server...')
tic
sp = io.getSpikes(thisSession, 'Kilo');
fprintf(' [%02.2f]\n', toc)

%% plot spike waveforms
ops = io.loadOps(thisSession);

nShanks = numel(ops);

for iShank = 1:nShanks
    fig(iShank) = figure(iShank+100); clf
    fig(iShank).Position = [1 1 floor((8.5/11)*1080) 1080];
    
    axes('Position', [.05 .6 .4 .35])
    
    if isfield(sp{iShank}, 'uQ')
        [~, ind] = sort(sp{iShank}.clusterDepths);
        clusterIds = sp{iShank}.cids(ind);
    else
        clusterIds = [];
    end
    fig(iShank) = plot.spikeWaveformsFromOps(ops(iShank), sp{iShank}, 'figure',  fig(iShank), 'clusterIds', clusterIds, 'numWaveforms', 100, 'yscale', 50);
    axis off
end




%% detect saccades
% This is quick and dirty. We should probably replace this with something
% more robust. I doubt we want to score every saccade in the GUI, but we
% might want to score enough of them to know that we trust the algorithm.

% load Eye position
[data, timestamps, elInfo] = io.getEdf(thisSession, PDS, false);

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
    'minIsi', ceil(25/elInfo.sampleRate*1e3), ...
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
subplot(1,2,1)

% saccade vectors
dx = saccades.endXpos - saccades.startXpos;
dy = saccades.endYpos - saccades.startYpos;

[cnt, binEdges] = hist3([dx(:) dy(:)], [40 40]);

binCenters{1} = binEdges{1}+mean(diff(binEdges{1})/2);
binCenters{2} = binEdges{2}+mean(diff(binEdges{2})/2);

imagesc(binCenters{1}, binCenters{2}, log(cnt))
title('Saccade Vectors', 'FontWeight', 'normal')
xlabel('Degrees')
ylabel('Degrees')
axis xy


[cnt, binEdges] = hist3([data(2,:)', data(1,:)'], [400 400]);
binCenters{1} = binEdges{1}+mean(diff(binEdges{1})/2);
binCenters{2} = binEdges{2}+mean(diff(binEdges{2})/2);

subplot(1,2,2)
imagesc(binCenters{1}, binCenters{2}, log(cnt))
title('Eye Position', 'FontWeight', 'normal')
xlabel('Degrees')
ylabel('Degrees')
axis xy

%%
figure(20); clf
imagesc(binCenters{1}, binCenters{2}, (cnt))
axis xy
grid on
set(gca, 'GridColor', 'r', 'GridAlpha', 1)


%% build stimulus object
spatialMap = session.squareFlash(PDS, 'eyetrace', [timestamps'; data], 'saccades', saccades);

%% bin space and create the design matrix
%x,y,x,y
psthWindow  = [0 -20 40 20]; %[-3, 3]% degrees (window is the same in x and y -- TODO: separate)
binSize = 1; % degrees

stim = spatialMap.binSpace('window', psthWindow, 'binSize', binSize, 'correctEyePos', 'simple');

%% Bin up spike counts
T = size(stim.X,1);
spks = zeros(T,1);
% clusterIds = sp{1}.cids(sp{1}.cgs >=2);
units = clusterIds;
% [~, inds] = sort(sp{1}.clusterDepths);
% units = sp{1}.cids(inds);
nUnits = numel(units);
for i = 1:nUnits
    cnt = histc(sp{1}.st(sp{1}.clu == units(i)), stim.bins);
    cnt((diff(stim.bins) > (spatialMap.display.ifi*1.5))) = 0;
    spks(:,i) = cnt;
end

figure(3); clf
subplot(4,2,1:2)
imagesc(stim.X(stim.valid,:)')
title('Stimulus')
subplot(4,2,3:4)
imagesc(filter(ones(100,1)/100, 1, spks(stim.valid,:))')
title('Spikes')
subplot(2,2,3)
imagesc(cov(spks))
subplot(2,2,4)
errorbar(sum(spks), var(spks), 'o')

%% Quick check to make sure the mapping is working
% stim = stimBase;
% fit spatial RFs
temporal_lags = 1:8;
baseLambda = 100;
RFcorr = rfmap.spatialRfAutoSmooth(stim.X(stim.valid,:), spks(stim.valid,:), temporal_lags, stim.size, 'crossvalidation', false, 'lambda0', baseLambda);

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

% uncorrected
stimBase = spatialMap.binSpace('window', psthWindow, 'binSize', binSize, 'correctEyePos', 'no');
[RFbase, dcBase, rfBase_params] = rfmap.spatialRfAutoSmooth(stimBase.X(stimBase.valid,:), spks(stimBase.valid,:), temporal_lags, stimBase.size, 'crossvalidation', true, 'lambda0', baseLambda*1000);

% fit 2D gabors and gaussians to the corrected RFs
results = rfmap.fitParametric2dRf(stim, RFcorr);

%% plot the RFs
fits = {RFcorr, RFbase};
for iFit = 1:numel(fits)
    I = fits{iFit};
    figure(iShank + 100 + (iFit-1));
    sx = ceil(sqrt(nUnits));
    sy = round(sqrt(nUnits));
    ax = pdsa.tight_subplot(sx, sy, .05, [.65 .05], [.5 .01]);
    if isfield(rfCorr_params, 'fold_r2')
        r2_test = mean(rfCorr_params.fold_r2);
    else
        r2_test = ones(1, nUnits);
    end
    
    for i = 1:(sx*sy)
        set(gcf, 'currentaxes', ax(i))
        if i <= nUnits
            imagesc(stim.xax, stim.yax, reshape(I(:,i), stim.size)); colormap gray
            set(ax(i), 'gridcolor', 'y')
            grid on
            axis xy
            hold on
            text(mean(stim.xax), .8*max(stim.yax), sprintf('Unit: %d', i), 'Color', cmap(i,:))
            if  r2_test(i) > .01
                
                xy = [results(i).gaussian.means.x0, results(i).gaussian.means.y0];
                plot.plotellipse(xy, [results(i).gaussian.means.sigmax 0; 0 results(i).gaussian.means.sigmay], 1, 'Color', cmap(i,:), 'Linewidth', 2)
            end
        else
            set(ax(i), 'Visible', 'off')
        end
    end
end

%% check that the eye traces align

psa = session.psaForage(PDS);


%% add spikes to trial

targ1On = arrayfun(@(x) x.targsOn(1), psa.trial);
goodTrial = find(~isnan(targ1On));

% bin at ms resolution
binfun = @(x) (x==0) + ceil(x/1e-3);

clear trial
trial = repmat(struct(), numel(goodTrial),1);
for kTrial = 1:numel(goodTrial)
    thisTrial = goodTrial(kTrial);
    
    trial(kTrial).start = psa.trial(thisTrial).start;
    trial(kTrial).duration = (psa.trial(thisTrial).duration);
    
    trial(kTrial).targ1On = binfun(psa.trial(thisTrial).targsOn(1) - trial(kTrial).start);
    trial(kTrial).targ1Off = binfun(psa.trial(thisTrial).targsOff(1) - trial(kTrial).start);
    trial(kTrial).targ1Dir = psa.trial(thisTrial).targDirection(1);
    
    trial(kTrial).targ2On = binfun(psa.trial(thisTrial).targsOn(2) - trial(kTrial).start);
    trial(kTrial).targ2Off = binfun(psa.trial(thisTrial).targsOff(2) - trial(kTrial).start);
    trial(kTrial).targ2Dir = psa.trial(thisTrial).targDirection(2);
    
    trial(kTrial).choice = psa.trial(thisTrial).targChosen;
    trial(kTrial).choiceTime = binfun(psa.trial(thisTrial).choiceTime - trial(kTrial).start);
    
    for kUnit = 1:nUnits
        st = sp{iShank}.st(sp{iShank}.clu == clusterIds(kUnit));
        st = st - trial(kTrial).start;
        trial(kTrial).(sprintf('sptrain%02.0f', kUnit)) = binfun(st(st > 0 & st < trial(kTrial).duration));
    end
    
    trial(kTrial).start = binfun(trial(kTrial).start);
    trial(kTrial).duration = binfun(psa.trial(thisTrial).duration);
end

t0 = trial(1).start;
for kTrial = 1:numel(trial)
    trial(kTrial).start = trial(kTrial).start - t0;
end

%%    
trial(4)
save('sample_trial.mat', '-v7.3', 'trial')

% trial = 
%%


iix = timestamps > psa.trial(1).start & timestamps < psa.trial(end).start;

[cnt, binCenters] = hist3([data(2,iix)', data(1,iix)'], [400 400]);
% binCenters{1} = binEdges{1}+mean(diff(binEdges{1})/2);
% binCenters{2} = binEdges{2}+mean(diff(binEdges{2})/2);

figure(iShank + 100);
axes('Position', [.05 .35 .25 .2])
imagesc(binCenters{2}, binCenters{1}, log(cnt)); colormap gray
hold on
title('Eye Position', 'FontWeight', 'normal')
xlabel('Degrees')
ylabel('Degrees')
axis xy
grid on
set(gca, 'GridColor', 'y', 'GridAlpha', .25)

% for kTrial = 1:psa.numTrials
%     plot(psa.trial(kTrial).eyePosAtFrame(:,1), psa.trial(kTrial).eyePosAtFrame(:,2), '.'); hold on
% end

numTargs = max(arrayfun(@(x) numel(x.targets), psa.trial));
cmap = lines;
for kUnit = 1:nUnits
	if r2_test(kUnit) < .01
        continue
    end
    
    xy = [results(kUnit).gaussian.means.x0, results(kUnit).gaussian.means.y0];
    plot.plotellipse(xy, [results(kUnit).gaussian.means.sigmax 0; 0 results(kUnit).gaussian.means.sigmay], 1, 'Color', cmap(kUnit,:), 'Linewidth', 2); hold on
end
    
targX = cell2mat(arrayfun(@(x) x.targPosX, psa.trial, 'uni', false));
targY = cell2mat(arrayfun(@(x) x.targPosY, psa.trial, 'uni', false));

cTarg = hsv(2);
for kTarg = 1:numTargs
    plot(targX(:,kTarg), targY(:,kTarg), 'o', 'Color', 'k', 'MarkerFaceColor', cTarg(kTarg,:))
end
    
plot(0, 0, 'or')
grid on
xlim([-15 15])
ylim([-10 10])
xlabel('Degrees')
ylabel('Degrees')
% title('RFs / Task')


targ1On = arrayfun(@(x) x.targsOn(1), psa.trial);
goodTrial = ~isnan(targ1On);
spikeCount = zeros(sum(goodTrial), nUnits);

countingWindow = .2;
psthWindow = [-.2 1]; % centered on the event
binSize = .01;
smoothingKernel = ones(1,3)/3; % boxcar

axes('Position', [.35 .35 .3 .2])
for kUnit = 1:nUnits
    if r2_test(kUnit) < .01
        continue
    end
    clix = sp{1}.clu == clusterIds(kUnit);
    st = sp{1}.st(clix);
    
    [m, s, bc, v, tspcnt] = pdsa.eventPsth(st, targ1On(goodTrial), psthWindow, binSize, smoothingKernel);
    spikeCount(:,kUnit) = pdsa.countSpikes(st, targ1On(goodTrial), countingWindow);
    errorbarFill(bc, m, s, cmap(kUnit,:), 'EdgeColor', 'none', 'FaceAlpha', .5); hold on
    axis tight
end

xlabel('Time (sec)')
ylabel('Spike Rate')
title('PSTH')

set(gcf, 'Color', 'w')


t1cho = arrayfun(@(x) x.targChosen == 1, psa.trial);
t2cho = arrayfun(@(x) x.targChosen == 2, psa.trial);
fprintf('Target 1 chosen %d times\n', sum(t1cho))
fprintf('Target 2 chosen %d times\n', sum(t2cho))

t1direction = arrayfun(@(x) x.targDirection(1), psa.trial(goodTrial));
t2direction = arrayfun(@(x) x.targDirection(2), psa.trial(goodTrial));

directions = unique([t1direction; t2direction]);
nDirections = numel(directions);

Tc1 = nan(nDirections, nUnits);
Tc2 = nan(nDirections, nUnits);
Tc1err = nan(nDirections, nUnits);
Tc2err = nan(nDirections, nUnits);

for iDirection = 1:nDirections
    iix = t1direction==directions(iDirection);
    Tc1(iDirection, :) = mean(spikeCount(iix,:));
    Tc1err(iDirection, :) = std(spikeCount(iix,:))/sqrt(sum(iix));
    
    iix = t2direction==directions(iDirection);
    Tc2(iDirection, :) = mean(spikeCount(iix,:));
    Tc2err(iDirection, :) = std(spikeCount(iix,:))/sqrt(sum(iix));
end
%%
ax1 = axes('Position', [.7 .45 .25 .1]);
ax2 = axes('Position', [.7 .35 .25 .1]);
for kUnit = 1:nUnits
    set(gcf, 'currentaxes', ax1)
    errorbarFill(directions, Tc1(:,kUnit), Tc1err(:,kUnit), cmap(kUnit,:), 'EdgeColor', 'none', 'FaceAlpha', .5); hold on
    plot(directions, Tc1(:,kUnit), 'Color', cmap(kUnit,:), 'Linewidth', 2)
    xlabel('Direction')
    ylabel('Spike Count')
%     title('Target 1')
    axis tight
    
    set(gcf, 'currentaxes', ax2)
    errorbarFill(directions, Tc2(:,kUnit), Tc2err(:,kUnit), cmap(kUnit,:), 'EdgeColor', 'none', 'FaceAlpha', .5); hold on
    plot(directions, Tc2(:,kUnit), 'Color', cmap(kUnit,:), 'Linewidth', 2)
    xlabel('Direction')
    ylabel('Spike Count')
%     title('Target 2')
    axis tight
end

fname = [thisSession.Subject{1} '_' datestr(datenum(thisSession.Date{1}), 'yyyymmdd') '.pdf'];    
set(gcf, 'PaperSize', [8.5 11], 'PaperPosition', [0 0 8.5 11])
saveas(gcf, fname)
%%
























%%


%% Flash / Saccade triggered CSD
ops = io.loadOps(thisSession);
% pick a shank and load LFP
iShank = 1;
[lfp, lfpTime, lfpInfo] = io.getLFP(ops(iShank));
lfpInfo.fragments = round(lfpInfo.fragments);

csd = session.csdFlash(PDS);

if ~isempty(csd.trial)
csd = session.csdFlash(PDS);

ops = io.loadOps(thisSession);
stats = csd.computeCsd(ops);
    
figure(1); clf
subplot(1,2,1)

imagesc(stats.time, stats.depth, stats.CSD); colormap jet
hold on
plot(stats.time, bsxfun(@plus, stats.STA, stats.chDepths'), 'Color', .5*[1 1 1])

subplot(1,2,2)
plot(stats.time, stats.STA)
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

%%
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
     subplot(1,2,2);
%     ax1 = subplot(1,2,2);
%     set(gcf, 'currentaxes', ax1)
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
%     set(gcf, 'currentaxes', ax1)
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
    goodUnits = s.cgs >= 1;
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
    clustIds = s.cids(s.cgs>=1);
else
    clustIds = s.cids;
end
nUnits = numel(clustIds);

figure(2); clf
ax = pdsa.tight_subplot(nUnits, 2, .01, .01);
for kUnit = 1:nUnits
    spikeTimes = s.st(s.clu==clustIds(kUnit));
    set(gcf, 'currentaxes', ax((kUnit-1)*2 + 1))
%     pdsa.plotRaster(spikeTimes, [hart.trial.start], [-.1 20], .01)
%     pause
% end

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
s.cgs(:) = 3;
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