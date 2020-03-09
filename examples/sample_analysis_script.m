
%% Sample Analysis Script
% This script shows an example of how to load data into the workspace
% (assuming it has been properly imported)

% --- add paths
cd C:\Users\Jake\Repos\ephys-matlab\
addEphysMatlab

%% set session directory
oepath = 'C:\Data\Ellie_2017-07-31_15-21-11_Shnkd8';
% oepath = 'C:\Data\Ellie_2017-08-09_13-04-23_ShankD15MT6';
%oepath = 'C:\Data\Ellie_2017-08-08_13-41-34_Shank2D14MT5';
 oepath = uigetdir(); % select using GUI

%% Load relevant data into the workspace
[sess, ops, info] = io.loadSession(oepath);

PDS = io.getPds(sess);

sp = io.getSpikes(sess);

%% plot spike waveforms?
nShanks = numel(ops);

for iShank = 1:nShanks
    figure(iShank); clf
    if isfield(sp{iShank}, 'uQ')
        clusterIds = sp{iShank}.cids(sp{iShank}.uQ>10);
    else
        clusterIds = [];
    end
    fig(iShank) = plot.spikeWaveformsFromOps(ops(iShank), sp{iShank}, 'figure', iShank, 'clusterIds', clusterIds, 'numWaveforms', 100);
end


%% Spatial Mapping: square flash

% build trial structure
spatialMap = session.squareFlash(PDS);

help session.squareFlash/binSpace % tells you how to use the function
%% bin space at the resolution you are interested
window = [-10 10]; %[-3, 3]% degrees (window is the same in x and y -- TODO: separate)
binSize = 2; %.5; % degrees

stim = spatialMap.binSpace('window', window, 'binSize', binSize, 'correctEyePos', true);

nTimeBins = 20; % frames
Xd = spatialMap.buildDesignMatrix(stim, nTimeBins);


%% Plot spatial map for each unit
for iShank = 1
s = sp{iShank};
if isfield(s, 'uQ')
    clustIds = s.cids(s.uQ>10);
else
    clustIds = s.cids;
end
nUnits = numel(clustIds);
figure(iShank*10 + 1); clf
sx = ceil(sqrt(nUnits));
sy = round(sqrt(nUnits));
ax = pdsa.tight_subplot(sy,sx,.001, .001);

figure(iShank*10 + 2); clf
ax2 = pdsa.tight_subplot(sy,sx,.001, .001);
for kUnit = 1:nUnits
spikeTimes = s.st(s.clu==clustIds(kUnit));

sta = spatialMap.spikeTriggeredAverage(stim, Xd, spikeTimes);

figure(iShank*10 + 1);
set(gcf, 'currentaxes', ax(kUnit))
xgrid = stim.xax(1:2:end);


% ax = subplot(1,2,1);
imagesc(stim.xax, stim.yax, sta.RF); axis xy
colormap gray
hold on

% plot([xgrid; xgrid], repmat(xgrid([1 end])', 1, numel(xgrid)), 'y')
% plot(repmat(xgrid([1 end])', 1, numel(xgrid)), [xgrid; xgrid], 'y')
plot(xgrid,zeros(size(xgrid)), '+y');
plot(zeros(size(xgrid)), xgrid, '+y');
grid on
set(gca, 'GridColor','y');
set(gca, 'GridAlpha', .25);
axis off
% title('Spatial Map', 'fontweight', 'normal')
% xlabel('degrees')
% ylabel('degrees')
% subplot(1,2,2)

figure(iShank*10 + 2)
set(gcf, 'currentaxes', ax2( kUnit))
plot(sta.time, sta.RFtime)
% xlabel('Time (seconds)')
% title(sprintf('Unit: %d', kUnit))
axis off
drawnow
end
end
%% detect saccades
% This is quick and dirty. We should probably replace this with something
% more robust. I doubt we want to score every saccade in the GUI, but we
% might want to score enough of them to know that we trust the algorithm.

% load Eye position
[data, timestamps, elInfo] = io.getEdf(sess, PDS);

% --- remove bad samples
ix = any(data(1:2,:)==elInfo.bitDeg(2));
data(1,ix) = nan;
data(2,ix) = nan;

% --- detect those saccades
ix = ~any(isnan(data));
[saccades] = pdsa.detectSaccades(timestamps(ix), data(1:2,ix), 'verbose', false);

% --- plot output
figure(1); clf
plot(timestamps, data(1,:), 'k'); hold on
plot(timestamps, data(2,:), 'Color', repmat(.5, 1, 3));


tbins = timestamps(1:3e3:end);
cnt = histc(saccades(1,:), tbins);
[~, id] = max(cnt);

plot([saccades(1,:); saccades(1,:)], ylim, 'r--')
xlim(tbins(id + [0 1]))
xlabel('seconds')
ylabel('degrees')
legend({'x', 'y', 'saccade start'})

figure(2); clf
plot(saccades(3,:)/elInfo.sampleRate, saccades(4,:), '.')
xlabel('duration (s)')
ylabel('amplitude (deg)')
% lfpInfo.timestamps = lfpInfo.timestamps-lfpInfo.timestamps(1);
%% Flash / Saccade triggered CSD

% pick a shank and load LFP
iShank = 1;
[lfp, lfpTime, lfpInfo] = io.getLFP(ops(iShank));

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
%%
% --- saccade-triggered
eventTimes = [saccades.start]';
% find saccades that happened more than 200 ms after the last saccade
% (uncontaminated)
ixLong = find([0 diff(eventTimes)>.2]);
lfpInfo.fragments = ceil(lfpInfo.fragments);
if ops(iShank).Nchan > 10
    subplot(1,2,2)
    plot.CsdBasic(lfp, eventTimes(ixLong), lfpInfo)
    title('Saccade Triggered')
    xlabel('ms')
    ylabel('depth')
end

% Look for saccade-direction tuning in the LFP
th = cart2pol([saccades.endXpos] - [saccades.startXpos], [saccades.endYpos] - [saccades.startYpos])'/pi*180;

if ops(iShank).Nchan > 10
    figure(2); clf
end

et = eventTimes(ixLong)';
ev = io.convertTimeToSamples(et, lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));
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
    
end

%% ********** spiking locked to saccades **********
iShank = 1;
figure(3); clf
s = sp{iShank};
if isfield(s, 'uQ')
    goodUnits = s.uQ > 1;
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
if isfield(s, 'uQ')
    clustIds = s.cids(s.uQ>1);
else
    clustIds = s.cids;
end
nUnits = numel(clustIds);

figure(1); clf
% plot.spikeWaveformsFromOps(ops(iShank), s, 'clusterIds', clustIds, 'numWaveforms', 10)
axis off
figure(2); clf
ax = pdsa.tight_subplot(nUnits, 2, .01, .01);
for kUnit = 1:nUnits
    spikeTimes = s.st(s.clu==clustIds(kUnit));
    
    sta = hart.spikeTriggeredAverage(spikeTimes);
    
    xdat.xx = hart.design.XX;
    xdat.xy = sta.xy;
    xdat.yy = sta.yy;
    xdat.ny = sta.ny;
%     autoCorrRidgeRegress(xdat)
%     sta = hart.AsdRf(spikeTimes);
    [w_hat, wt, wx, wlin] = bilinearMixRegress_coordAscent(xdat.xx + 10e2*speye(size(xdat.xx,2)), xdat.xy, [hart.design.nkTime hart.design.nkx*hart.design.nky], 1, 1:size(xdat.xx,2)-1);
    sflip = sign(sum(wt));
    sta.RF = reshape(sflip*wx, [hart.design.nkx hart.design.nky]);
    sta.RFtime = sflip*wt;
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

