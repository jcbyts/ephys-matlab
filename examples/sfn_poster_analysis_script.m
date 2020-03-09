
%% Sample Analysis Script
% This script shows an example of how to load data into the workspace
% (assuming it has been properly imported)

% --- add paths
cd C:\Users\Jake\Repos\ephys-matlab\
addpath C:\Users\Jake\Dropbox\MatlabCode\Repos\fitglmqp130\
addEphysMatlab

%% set session directory
%oepath = 'C:\Data\Ellie_2017-07-31_15-21-11_Shnkd8';
% oepath = 'C:\Data\Ellie_2017-08-09_13-04-23_ShankD15MT6';
%oepath = 'C:\Data\Ellie_2017-08-08_13-41-34_Shank2D14MT5';
 oepath = uigetdir(); % select using GUI

%% Load relevant data into the workspace
[sess, ops, info] = io.loadSession(oepath);

PDS = io.getPds(sess);

sp = io.getSpikes(sess);

%% plot spike waveforms?
nShanks = numel(ops);

iShank = 2;
s = sp{iShank};
if isfield(s, 'uQ')
    th =  0;
    clustDepths = s.clusterDepths(s.uQ>th);
    clustIds = s.cids(s.uQ>th);
    [~, id] = sort(clustDepths);
    clustIds = clustIds(id);
else
    clustIds = s.cids;
end
% clustIds = clustIds([2 4 5 6 7 10 11 12 13 14 15]);
    
figure(iShank); clf
% 
fig(iShank) = plot.spikeWaveformsFromOps(ops(iShank), sp{iShank}, 'figure', iShank, 'clusterIds', clustIds, 'numWaveforms',10, 'useMean', false);
% pdsa.fixfigure(gcf, 12, [2 6])
% saveas(gcf, 'waveforms.pdf')

%% Spatial Mapping: square flash

% build trial structure
spatialMap = session.squareFlash(PDS);

help session.squareFlash/binSpace % tells you how to use the function
% bin space at the resolution you are interested
window = [-8 8]; %[-3, 3]% degrees (window is the same in x and y -- TODO: separate)
binSize = 1; %.5; % degrees

stim = spatialMap.binSpace('window', window, 'binSize', binSize, 'correctEyePos', true);

nTimeBins = 20; % frames
% Xd = spatialMap.buildDesignMatrix(stim, nTimeBins);


% Plot spatial map for each unit

% clustIds(1) = [];
nUnits = numel(clustIds);

spikeTimes = cell(nUnits,1);
for kUnit = 1:nUnits
spikeTimes{kUnit} = s.st(s.clu==clustIds(kUnit));
end

% S = spatialMap.AsdSpaceTimeSep(stim, spikeTimes);

% S = spatialMap.AutoRidgeSpaceTimeSep(stim, spikeTimes);
%
S = spatialMap.LsSmoothSpaceTimeSep(stim, spikeTimes, 500);
%  ---- Make Plots ----

figure(iShank); clf
fig(iShank) = plot.spikeWaveformsFromOps(ops(iShank), sp{iShank}, 'figure', iShank, 'clusterIds', clustIds, 'numWaveforms',50, 'useMean', true);

cmap = lines;
figure(33); clf
ax = pdsa.tight_subplot(nUnits, 2, .001, .001);
mx = max(cellfun(@(x) max(x.RF(:)), S));
mn = min(cellfun(@(x) min(x.RF(:)), S));

for kUnit = 1:nUnits
    
    sta = S{kUnit};
    set(gcf, 'currentaxes', ax( (kUnit-1)*2 +1))
    xgrid = stim.xax(1:2:end);
    
    
    % ax = subplot(1,2,1);
    clim = [max(mn, min(sta.RF(:))) min(mx, max(sta.RF(:))*2)];
%     clim = [mn/4 mx/4];
%     imagesc(stim.xax, stim.yax, sta.RF, clim); axis xy
    imagesc(stim.xax, stim.yax, sta.RF); axis xy
    colormap(gray.^2)
    hold on
    
    % plot([xgrid; xgrid], repmat(xgrid([1 end])', 1, numel(xgrid)), 'y')
    % plot(repmat(xgrid([1 end])', 1, numel(xgrid)), [xgrid; xgrid], 'y')
    plot(xgrid,zeros(size(xgrid)), '+y');
    plot(zeros(size(xgrid)), xgrid, '+y');
    grid on
    set(gca, 'GridColor','y');
    set(gca, 'GridAlpha', .25);
    xlim(stim.xax([1 end]))
    ylim(stim.xax([1 end]))
    axis off
    % title('Spatial Map', 'fontweight', 'normal')
    % xlabel('degrees')
    % ylabel('degrees')
    % subplot(1,2,2)
    set(gcf, 'currentaxes', ax( (kUnit-1)*2 +2))
    plot(sta.time, sta.RFtime, 'Color', cmap(kUnit,:))
    % xlabel('Time (seconds)')
    % title(sprintf('Unit: %d', kUnit))
    if kUnit < nUnits
    axis off
    end
    drawnow
end

set(iShank, 'Color', 'w')
set(33, 'Color', 'w')

%%
pdsa.fixfigure(33, 10, [2 nUnits])
saveas(33, 'spatialmaps.pdf')
% detect saccades
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
% Flash / Saccade triggered CSD

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
    plot.CsdBasic((lfp), flashTimes, lfpInfo, 'skipChannels', 2)
    title('Flash Triggered')
    xlabel('ms')
    ylabel('depth')
end

% --- saccade-triggered
eventTimes = saccades(1,:);
% find saccades that happened more than 200 ms after the last saccade
% (uncontaminated)
ixLong = find([0 diff(eventTimes)>.2]);
ch0 = (1:32)*40;
iix =1:2:numel(ch0);
if ops(iShank).Nchan > 10
    subplot(1,2,2)
    plot.CsdBasic((lfp), eventTimes(ixLong), lfpInfo, 'skipChannels', 2)
    title('Saccade Triggered')
    xlabel('ms')
    ylabel('depth')
end

% Look for saccade-direction tuning in the LFP
th = cart2pol(saccades(7,:) - saccades(5,:), saccades(8,:) - saccades(6,:))'/pi*180;

if ops(iShank).Nchan > 10
    figure(2); clf
end

et = eventTimes(ixLong)';
ev = io.convertTimeToSamples(et, lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));
th = th(ixLong(:));
thBins = 0:45:360;

cmap = flipud(hsv(numel(thBins)));
if ops(iShank).Nchan > 10
    subplot(1,2,2)
    for i = 1:numel(thBins)-1
        
        ii = th > thBins(i) & th < thBins(i+1);
        
        [sta,~,bc] = pdsa.eventTriggeredAverage(lfp(:,iix), ev(ii), [-300 400]);
        
        plot(bc, bsxfun(@plus, sta, ch0(iix)), 'Color', cmap(i,:)); hold on
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

pdsa.fixfigure(gcf, 10, [7 4])
saveas(gcf, 'CSD.pdf')
% ********** spiking locked to saccades **********

figure(3); clf
s = sp{iShank};
% if isfield(s, 'uQ')
%     goodUnits = s.uQ > 10;
%     udepths = s.clusterDepths(goodUnits);
%     uId = s.cids(goodUnits);
%     [~, depthId] = sort(udepths);
%     clustId = uId(depthId);
% else
%     clustId = s.cids;
% end

clustId = clustIds;
nUnits = numel(clustId); % including

ax = pdsa.tight_subplot(nUnits, 1, 0.01,  0.01);
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

pdsa.fixfigure(gcf, 10, [2 nUnits])
saveas(gcf, 'sacTuning.pdf')

%

% creat stimulus struct
stim = 'DotMotionMapping';
hasStim = io.findPDScontainingStimModule(PDS, stim);

dotMappingTrial = struct();
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
        dotMappingTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
        dotMappingTrial(kTrial).start      = dotMappingTrial(kTrial).frameTimes(1);
        dotMappingTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - dotMappingTrial(kTrial).start;
        
        
        onset = find(diff(PDS{i}.data{thisTrial}.DotMotionMapping.on)==1) + 1;
        dotMappingTrial(kTrial).motionon  = dotMappingTrial(kTrial).frameTimes(onset);
        offset = min(find(diff(PDS{i}.data{thisTrial}.DotMotionMapping.on)==-1) + 1, numel(dotMappingTrial(kTrial).frameTimes));
        dotMappingTrial(kTrial).motionoff = dotMappingTrial(kTrial).frameTimes(offset);
        
        if PDS{i}.data{thisTrial}.DotMotionMapping.on(find(isnan(PDS{i}.data{thisTrial}.DotMotionMapping.on), 1)-1)
            dotMappingTrial(kTrial).motionoff = [dotMappingTrial(kTrial).motionoff dotMappingTrial(kTrial).frameTimes(end)];
        end
        dotMappingTrial(kTrial).direction = PDS{i}.data{thisTrial}.DotMotionMapping.direction(onset);  % motion directions
        dotMappingTrial(kTrial).speed = PDS{i}.data{thisTrial}.DotMotionMapping.speed(onset);  % motion speed
    end
end

%
%
figure(11); clf
figure(12); clf
figure(13); clf
figure(14); clf
sx = 1; %ceil(sqrt(nUnits));
sy = nUnits; %round(sqrt(nUnits));

clusterIds = clustIds;
for kUnit = 1:nUnits
% kUnit = 0;
% %%
% kUnit = kUnit + 1;
% clusterIds = unique(sp{1}.clu);

spikeTimes = sp{1}.st(sp{1}.clu==clusterIds(kUnit));

motionOnset  = cell2mat(arrayfun(@(x) x.motionon(:), dotMappingTrial(:), 'UniformOutput', false));
motionOffset = cell2mat(arrayfun(@(x) x.motionoff(:), dotMappingTrial(:), 'UniformOutput', false));
direction   = cell2mat(arrayfun(@(x) x.direction(:), dotMappingTrial(:), 'UniformOutput', false));
speed = cell2mat(arrayfun(@(x) x.speed(:), dotMappingTrial(:), 'UniformOutput', false));
speeds = unique(speed);
directions = unique(direction);
nDirections = numel(directions);
nTrials = numel(motionOnset);

cmap = [zeros(nDirections, 1) cosd(directions(:)-225) sind(directions(:)-225)];
cmap = (cmap + 1.2)/2.5;
cmap(:,1) = 0;
[~, trialId] = sort(direction);
binSize = 1e-3;
[spcnt, bcenters, valid] = pdsa.binSpTimes(spikeTimes, motionOnset, [-.2 .6], binSize);

win = [.01 .3];
[spcntTransient] = pdsa.binSpTimes(spikeTimes, motionOnset, win, binSize);
spcntTransient = sum(spcntTransient,2)/diff(win);

win = [.1 .5];
[spcntSustained] = pdsa.binSpTimes(spikeTimes, motionOnset, win, binSize);
spcntSustained = sum(spcntSustained,2)/diff(win);

spsort = spcnt(trialId,:);
% cmap = parula(nDirections);
% cmap = jet(nDirections);

figure(11);
% set(gcf, 'Color', 'w')
subplot(sy, sx, kUnit)
for i = 1:nTrials
    sptimes = find(spsort(i,:));
    plot([1; 1]*bcenters(sptimes), [i; i+3]*ones(1,numel(sptimes)), 'Color', cmap(direction(trialId(i))==directions,:)); hold on
end
drawnow
ds = direction(trialId);
yt = round(linspace(1,nTrials, 10));
set(gca, 'YTick', yt, 'YTickLabel', ds(yt))
% axis tight
ylim([0 nTrials])
ylabel('Direction')
xlabel('Time (s)')
xlim(bcenters([1 end]))

mu = nan(nDirections, 2);
se = nan(nDirections, 2);
if kUnit < nUnits
    set(gca, 'XTickLabel', '')
    xlabel('')
end

figure(12);
% set(gcf, 'Color', 'w')
subplot(sy, sx, kUnit)

for iTh = 1:nDirections
    trialix = direction==directions(iTh);
    skern = ones(50,1)/50;
    sm = filter(skern, 1, mean(spcnt(trialix,:)))/binSize;
    plot(bcenters, sm, 'Color', cmap(iTh,:)); hold on
    axis tight
    
    mu(iTh,1) = mean(spcntTransient(trialix));
    mu(iTh,2) = mean(spcntSustained(trialix));
    se(iTh,1) = std(spcntTransient(trialix))/sqrt(sum(trialix));
    se(iTh,2) = std(spcntSustained(trialix))/sqrt(sum(trialix));
end
xlim(bcenters([1 end]))
ylabel('Firing Rate')
xlabel('Time (s)')
if kUnit < nUnits
    set(gca, 'XTickLabel', '')
    xlabel('')
end


figure(13);
% set(gcf, 'Color', 'w')
subplot(sy, sx, kUnit)
plot(directions, mu(:,1), 'k'); hold on
plot(directions, mu(:,1)-se(:,1), 'k');
plot(directions, mu(:,1)+se(:,1), 'k');
% errorbar(directions, mu(:,1), se(:,1)); hold on
% errorbar(directions, mu(:,2), se(:,2));
xlim([0 360])
set(gca, 'XTick', 0:90:360)
ylabel('Firing Rate')
xlabel('Direction')

if kUnit < nUnits
    set(gca, 'XTickLabel', '')
    xlabel('')
end

figure(14);
% set(gcf, 'Color', 'w')
subplot(sy, sx, kUnit)

hh = polarplot(directions([1:end 1])/180*pi, mu(([1:end 1]),1), 'k'); hold on
% polarplot(directions([1:end 1])/180*pi, mu(([1:end 1]),2));
if kUnit < nUnits
    set(gca, 'ThetaTickLabel', '', 'RTickLabel', '')
end

end


pdsa.fixfigure(11, 18, 2*[1 nUnits], 'FontName', 'Arial')
pdsa.fixfigure(12, 18, 2*[1 nUnits], 'FontName', 'Arial')
pdsa.fixfigure(13, 18, 2*[1 nUnits],'FontName', 'Arial')
set(14, 'paperSize', 2*[1 nUnits], 'PaperPosition', [0 0 2*[1 nUnits]])
set(11, 'renderer', 'painters')
set(12, 'renderer', 'painters')
set(13, 'renderer', 'painters')
set(14, 'renderer', 'painters')
saveas(11, 'MT1.pdf')
saveas(12, 'MT2.pdf')
saveas(13, 'MT3.pdf')
saveas(14, 'MT4.pdf')

