

meta = io.getExperimentsAnd();

% an example session that has both coarse spatial mapping and fine
thisSession = meta(183,:);

%%
PDS = io.getPds(thisSession);
spikes = io.getSpikes(thisSession,'Kilo');

%%
ops = io.loadOps(thisSession);
wf = io.getSpikeWaveformsFromOps(ops, spikes{1});

%%

[C, F] = rfmap.CoarseToFineRF(PDS, spikes{1});
C.stim.xax = -C.stim.xax;
F.stim.xax = -F.stim.xax;
%% make stimulus movie
frameIndex = [200 300];
Ibig = pds.replayTrial(PDS{2}, 1, 'frameIndex', frameIndex);
Ismall = pds.replayTrial(PDS{3}, 1, 'frameIndex', [200 300]);
%%
V = double(squeeze(mean(Ibig{1},3)));
ctr = PDS{1}.initialParametersMerged.display.ctr;
ppd = PDS{1}.initialParametersMerged.display.ppd;
nFrames = size(V,3);
px = ((1:1920) - ctr(1)) / ppd;
py = ((1:1080) - ctr(2)) / -ppd;
fr = (1:nFrames);
[xx,yy,zz] = meshgrid(px, py, fr);

%%

figure(1); clf


iFrame = 80;
x = (PDS{2}.data{1}.behavior.eyeAtFrame(1,1:frameIndex(1)+iFrame)-ctr(1))/ppd;
y = (PDS{2}.data{1}.behavior.eyeAtFrame(2,1:frameIndex(1)+iFrame)-ctr(2))/-ppd;

ix = x < -20 | y < -10;
x(ix) = nan;
y(ix) = nan;
imagesc(px, py(20:end), Ibig{1}(20:end,:,:,iFrame)); axis xy
hold on
t = 70;
plot(x(end-t:end), y(end-t:end), 'c')
plot(x(end), y(end), '+r')
plot(C.ROI([1 1 3 3 1])+x(end), C.ROI([2 4 4 2 2])+y(end), 'r--') 
ylim([-15 15])
text(-10, 10, 'Gaze Contingent ROI', 'Color', 'r', 'FontSize', 12, 'FontName', 'Arial')
text(-10, 8.5, 'Eye Position', 'Color', 'c', 'FontSize', 12, 'FontName', 'Arial')
set(gca, 'Color', [.5 .5 .5])
axis off
plot(-19 + [-1 0], [-9 -9], 'k', 'Linewidth', 2)
text(-19 + -1, -8, '1°')
set(gcf, 'PaperSize', [5 4], 'PaperPosition', [0 0 5 4])
saveas(gcf, 'fig_spatial_stim_big.pdf')


figure(1); clf


iFrame = 80;
x = (PDS{3}.data{1}.behavior.eyeAtFrame(1,1:frameIndex(1)+iFrame)-ctr(1))/ppd;
y = (PDS{3}.data{1}.behavior.eyeAtFrame(2,1:frameIndex(1)+iFrame)-ctr(2))/-ppd;

ix = x < -20 | y < -10;
x(ix) = nan;
y(ix) = nan;
imagesc(px, py(20:end), Ismall{1}(20:end,:,:,iFrame)); axis xy
hold on
t = 70;
plot(x(end-t:end), y(end-t:end), 'c')
plot(x(end), y(end), '+r')
plot(F.ROI([1 1 3 3 1])+x(end), F.ROI([2 4 4 2 2])+y(end), 'r--') 
ylim([-15 15])
text(-10, 10, 'Gaze Contingent ROI', 'Color', 'r', 'FontSize', 12, 'FontName', 'Arial')
text(-10, 8.5, 'Eye Position', 'Color', 'c', 'FontSize', 12, 'FontName', 'Arial')
set(gca, 'Color', [.5 .5 .5])
axis off
plot(-19 + [-1 0], [-9 -9], 'k', 'Linewidth', 2)
text(-19 + -1, -8, '1°')
set(gcf, 'PaperSize', [5 4], 'PaperPosition', [0 0 5 4])
saveas(gcf, 'fig_spatial_stim_small.pdf')

%% 
figure(2); clf
unitIdx = [5 8 14]; %4:8 12:14];
nUnits = numel(unitIdx);
% ax1 = subplot(nUnits, 2, (find(mod(1:nUnits,2))-1)*2 + 1);

ax = pdsa.tight_subplot(nUnits, 2, 0.01, 0.1, 0.1);
cmap = lines(nUnits);
for i = 1:nUnits
    W=reshape(wf.wf(:,unitIdx(i)), numel(wf.wftax)-1, numel(wf.chy));
    sd = reshape(wf.sd(:,unitIdx(i)), numel(wf.wftax)-1, numel(wf.chy));
%     set(gcf, 'currentaxes', ax1)
%     subplot(nUnits, 2, (find(mod(1:nUnits,2))-1)*2 + 1)
%     plot(wf.wftax(1:end-1)+i*.002, bsxfun(@plus, W, wf.chy'), 'Color', cmap(i,:)); hold on
    
    
   set(gcf, 'currentaxes', ax((i-1)*2 + 1)); 
    plot(wf.wftax(1:end-1), bsxfun(@plus, W/3, wf.chy'), 'Color', cmap(i,:)); hold on
    plot(wf.wftax(1:end-1), bsxfun(@plus, W/3 + sd/3, wf.chy'), 'Color', cmap(i,:)); hold on
    ylim(spikes{1}.clusterDepths(unitIdx(i)) + [-150 110])
    
    
%     axis tight
    axis off
%     subplot(nUnits, 2, (i-1)*2 + 2)
    set(gcf, 'currentaxes', ax((i-1)*2 + 2));
    imagesc(C.stim.xax, C.stim.yax,reshape(C.RFcorr(:,unitIdx(i)), C.stim.size));
    hold on
    plot(0,0, '+y')
    axis xy
    axis off
    colormap gray
end
hold on
plot(C.stim.xax(end) + [1 2], C.stim.yax(1)  + [1 1], 'r')
text(C.stim.xax(end) + [1], C.stim.yax(1)  + 1*4, '1°', 'Color', 'r')
plot(F.ROI([1 1 3 3 1]), F.ROI([2 4 4 2 2]), 'r-') 
set(gcf, 'PaperSize', [3 nUnits*2], 'PaperPosition', [0 0 3 nUnits*2])
saveas(gcf, 'fig_spatial_rf_big.pdf')

%%

figure(3); clf
unitIdx = [5 8 14]; %4:8 12:14];
nUnits = numel(unitIdx);
% ax1 = subplot(nUnits, 2, (find(mod(1:nUnits,2))-1)*2 + 1);

ax = pdsa.tight_subplot(nUnits, 2, 0.01, 0.1, 0.1);
cmap = lines(nUnits);
for i = 1:nUnits
    W=reshape(wf.wf(:,unitIdx(i)), numel(wf.wftax)-1, numel(wf.chy));
    sd = reshape(wf.sd(:,unitIdx(i)), numel(wf.wftax)-1, numel(wf.chy));
%     set(gcf, 'currentaxes', ax1)
%     subplot(nUnits, 2, (find(mod(1:nUnits,2))-1)*2 + 1)
%     plot(wf.wftax(1:end-1)+i*.002, bsxfun(@plus, W, wf.chy'), 'Color', cmap(i,:)); hold on
    
    
    set(gcf, 'currentaxes', ax((i-1)*2 + 1)); 
    plot(wf.wftax(1:end-1), bsxfun(@plus, W/3, wf.chy'), 'Color', cmap(i,:)); hold on
    plot(wf.wftax(1:end-1), bsxfun(@plus, W/3 + sd/3, wf.chy'), 'Color', cmap(i,:)); hold on
    ylim(spikes{1}.clusterDepths(unitIdx(i)) + [-150 110])
    
%     axis tight
    axis off
%     subplot(nUnits, 2, (i-1)*2 + 2)
    set(gcf, 'currentaxes', ax((i-1)*2 + 2));
    imagesc(F.stim.xax, F.stim.yax,reshape(F.RFcorr(:,unitIdx(i)), F.stim.size)); axis xy
    hold on
%     axis off
    axis tight
    set(gca, 'XTick', -2:1:2, 'YTick', -2:1:2)
    plot(0,0, '+y')
    axis off
%     grid('on')
% plot.offsetAxes(gca)
    colormap gray
end
hold on
plot(F.stim.xax(end) + [0 1], F.stim.yax(1)  + [1 1], 'r')
text(F.stim.xax(end) + [0], F.stim.yax(1)  + 1.3, '1°', 'Color', 'r')
set(gcf, 'PaperSize', [3 nUnits*2], 'PaperPosition', [0 0 3 nUnits*2])
saveas(gcf, 'fig_spatial_rf_small.pdf')
% plot(F.ROI([1 1 3 3 1]), F.ROI([2 4 4 2 2]), 'r-') 

