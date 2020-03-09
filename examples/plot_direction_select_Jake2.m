oepath = uigetdir;

[sess, ops, info] = io.loadSession(oepath);

PDS = io.getPds(sess);
sp = io.getSpikes(sess);


% %% resort single channels
% close all
% % list single electrode channels (this can be a vector if > 1 electrode used)
% chanMap = 4;
% shank{1} = hardware.electrode.customChannelMap(chanMap);
% shank{1}.name = 'MtBurrHoleMapping';
% 
% singleTrodes = find(cellfun(@(x) numel(x.channelMap) < 80, shank));
% if any(singleTrodes)
%     clear sp
%     for i = 1:numel(singleTrodes)
%         iShank = singleTrodes(i);
%         sp{i} = preprocess.runSingleChannelSpikeSortThreshold(ops(iShank));
%     end
% end
%% creat stimulus struct
stim = 'MotionMapping';
hasStim = io.findPDScontainingStimModule(PDS, stim);

dotMappingTrial = struct();
trialNum = 0;

limitDirection = false;

for i = find(hasStim(:)')
    
    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 
    
    stimTrials = find(trialIx);
    
    if isempty(stimTrials)
        continue
    end
    
    initialParams = PDS{i}.initialParametersMerged.MotionMapping;
    if limitDirection  
        if isfield(initialParams, 'priorKappa') && initialParams.priorKappa == 0
            continue
        end
    else
        if isfield(initialParams, 'priorKappa') && initialParams.priorKappa > 0
            continue
        end
    end
    
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
        
        if ~PDS{i}.conditions{thisTrial}.MotionMapping.use
            continue
        end
    
        kTrial = trialNum + j;
        
        % --- Timing
        dotMappingTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
        dotMappingTrial(kTrial).start      = dotMappingTrial(kTrial).frameTimes(1);
        dotMappingTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - dotMappingTrial(kTrial).start;
        
        
        onset = find(diff(PDS{i}.data{thisTrial}.MotionMapping.on)==1) + 1;
        dotMappingTrial(kTrial).motionon  = dotMappingTrial(kTrial).frameTimes(onset);
        offset = min(find(diff(PDS{i}.data{thisTrial}.MotionMapping.on)==-1) + 1, numel(dotMappingTrial(kTrial).frameTimes));
        dotMappingTrial(kTrial).motionoff = dotMappingTrial(kTrial).frameTimes(offset);
        
        if PDS{i}.data{thisTrial}.MotionMapping.on(find(isnan(PDS{i}.data{thisTrial}.MotionMapping.on), 1)-1)
            dotMappingTrial(kTrial).motionoff = [dotMappingTrial(kTrial).motionoff dotMappingTrial(kTrial).frameTimes(end)];
        end
        dotMappingTrial(kTrial).direction = PDS{i}.data{thisTrial}.MotionMapping.direction(onset);  % motion directions
        dotMappingTrial(kTrial).speed = PDS{i}.data{thisTrial}.MotionMapping.speed(onset);  % motion speed
    end
end

%%
iShank = 1;
clustIds = sp{iShank}.cids;
%%
clusterIds = clustIds;
nUnits = numel(clustIds);
for kUnit = 1:nUnits
% kUnit = 0;
% %%
% kUnit = kUnit + 1;
% clusterIds = unique(sp{1}.clu);

spikeTimes = sp{iShank}.st(sp{iShank}.clu==clusterIds(kUnit));

motionOnset  = cell2mat(arrayfun(@(x) x.motionon(:), dotMappingTrial(:), 'UniformOutput', false));
motionOffset = cell2mat(arrayfun(@(x) x.motionoff(:), dotMappingTrial(:), 'UniformOutput', false));
direction   = cell2mat(arrayfun(@(x) x.direction(:), dotMappingTrial(:), 'UniformOutput', false));
speed = cell2mat(arrayfun(@(x) x.speed(:), dotMappingTrial(:), 'UniformOutput', false));
speeds = unique(speed);
directions = unique(direction);
nDirections = numel(directions);
nTrials = numel(motionOnset);

[~, trialId] = sort(direction);
binSize = 1e-3;
[spcnt, bcenters, valid] = pdsa.binSpTimes(spikeTimes, motionOnset, [-.2 .6], binSize);

win = [.01 .1];
[spcntTransient] = pdsa.binSpTimes(spikeTimes, motionOnset, win, binSize);
spcntTransient = sum(spcntTransient,2)/diff(win);

win = [.1 .5];
[spcntSustained] = pdsa.binSpTimes(spikeTimes, motionOnset, win, binSize);
spcntSustained = sum(spcntSustained,2)/diff(win);

spsort = spcnt(trialId,:);
cmap = parula(nDirections);

figure(kUnit); clf
set(gcf, 'Color', 'w')
subplot(4,2,1:4)
for i = 1:nTrials
    sptimes = find(spsort(i,:));
    plot(bcenters(sptimes), i*ones(numel(sptimes),1), '.', 'Color', cmap(direction(trialId(i))==directions,:)); hold on
%     plot([1; 1]*bcenters(sptimes), [i; i+3]*ones(1,numel(sptimes)), 'Color', cmap(direction(trialId(i))==directions,:)); hold on
end

ds = direction(trialId);
yt = round(linspace(1,nTrials, 10));
set(gca, 'YTick', yt, 'YTickLabel', ds(yt))
ylabel('Direction')
xlabel('Time (s)')
xlim(bcenters([1 end]))

mu = nan(nDirections, 2);
se = nan(nDirections, 2);
subplot(4,2,5:6)
for iTh = 1:nDirections
    trialix = direction==directions(iTh);
    skern = ones(10,1)/10;
    sm = filter(skern, 1, mean(spcnt(trialix,:)))/binSize;
    plot(bcenters, sm, 'Color', cmap(iTh,:)); hold on
    
    mu(iTh,1) = mean(spcntTransient(trialix));
    mu(iTh,2) = mean(spcntSustained(trialix));
    se(iTh,1) = std(spcntTransient(trialix))/sqrt(sum(trialix));
    se(iTh,2) = std(spcntSustained(trialix))/sqrt(sum(trialix));
end
xlim(bcenters([1 end]))
ylabel('Firing Rate')
xlabel('Time (s)')


subplot(4,2,7)
errorbar(directions, mu(:,1), se(:,1)); hold on
errorbar(directions, mu(:,2), se(:,2));
xlim([0 360])
set(gca, 'XTick', 0:90:360)
ylabel('Firing Rate')
xlabel('Direction')

subplot(4,2,8)
polarplot(directions([1:end 1])/180*pi, mu(([1:end 1]),1)); hold on
polarplot(directions([1:end 1])/180*pi, mu(([1:end 1]),2));

drawnow
end

%% tuning only
iShank = 2;

figure(iShank); clf
clustIds = sp{iShank}.cids;
%%
clusterIds = clustIds;
nUnits = numel(clustIds);
sx = ceil(sqrt(nUnits));
sy = round(sqrt(nUnits));

for kUnit = 1:nUnits

spikeTimes = sp{iShank}.st(sp{iShank}.clu==clusterIds(kUnit));

motionOnset  = cell2mat(arrayfun(@(x) x.motionon(:), dotMappingTrial(:), 'UniformOutput', false));
motionOffset = cell2mat(arrayfun(@(x) x.motionoff(:), dotMappingTrial(:), 'UniformOutput', false));
direction   = cell2mat(arrayfun(@(x) x.direction(:), dotMappingTrial(:), 'UniformOutput', false));
speed = cell2mat(arrayfun(@(x) x.speed(:), dotMappingTrial(:), 'UniformOutput', false));
speeds = unique(speed);
directions = unique(direction);
nDirections = numel(directions);
nTrials = numel(motionOnset);

[~, trialId] = sort(direction);
binSize = 1e-3;
[spcnt, bcenters, valid] = pdsa.binSpTimes(spikeTimes, motionOnset, [-.2 .6], binSize);

win = [.01 .1];
[spcntTransient] = pdsa.binSpTimes(spikeTimes, motionOnset, win, binSize);
spcntTransient = sum(spcntTransient,2)/diff(win);

win = [.1 .5];
[spcntSustained] = pdsa.binSpTimes(spikeTimes, motionOnset, win, binSize);
spcntSustained = sum(spcntSustained,2)/diff(win);

mu = nan(nDirections, 2);
se = nan(nDirections, 2);

for iTh = 1:nDirections
    trialix = direction==directions(iTh);
    skern = ones(10,1)/10;
    sm = filter(skern, 1, mean(spcnt(trialix,:)))/binSize;
    
    mu(iTh,1) = mean(spcntTransient(trialix));
    mu(iTh,2) = mean(spcntSustained(trialix));
    se(iTh,1) = std(spcntTransient(trialix))/sqrt(sum(trialix));
    se(iTh,2) = std(spcntSustained(trialix))/sqrt(sum(trialix));
end

subplot(sy,sx,kUnit)
errorbar(directions, mu(:,1), se(:,1)); hold on
errorbar(directions, mu(:,2), se(:,2));
xlim([0 360])
set(gca, 'XTick', 0:90:360)
ylabel('Firing Rate')
xlabel('Direction')
axis off
drawnow
end


%% Spatial Map?

%% Spatial Mapping: square flash

% build trial structure
spatialMap = session.squareFlash(PDS);

help session.squareFlash/binSpace % tells you how to use the function
%% bin space at the resolution you are interested
window = [-5 5]; %[-3, 3]% degrees (window is the same in x and y -- TODO: separate)
binSize = 1; %.5; % degrees

stim = spatialMap.binSpace('window', window, 'binSize', binSize, 'correctEyePos', true);

nTimeBins = 20; % frames
Xd = spatialMap.buildDesignMatrix(stim, nTimeBins);


%% Plot spatial map for each unit
iShank = 1;
s = sp{iShank};
if isfield(s, 'uQ')
    clustIds = s.cids(s.uQ>1);
else
    clustIds = s.cids;
end
nUnits = numel(clustIds);

kUnit = 1;
spikeTimes = s.st(s.clu==clustIds(kUnit));

sta = spatialMap.spikeTriggeredAverage(stim, Xd, spikeTimes, 10e3);

figure; clf
ax = subplot(1,2,1);
imagesc(stim.xax, stim.yax, sta.RF); axis xy
colormap jet
grid on
ax.GridColor = 'k';
ax.GridAlpha = 1;
title('Spatial Map', 'fontweight', 'normal')
xlabel('degrees')
ylabel('degrees')
subplot(1,2,2)
plot(sta.time, sta.RFtime)
xlabel('Time (seconds)')
title('Temporal RF')
% title(sprintf('Unit: %d', kUnit))
