function [ax] = session_epochs(stim, spikes)
% [ax, deets] = plot_session_epochs(stim, spikes)

set(gcf, 'Color', 'w')
bins = 0:.5:max(spikes{1}.st);
[depths, iii] = sort(spikes{1}.clusterDepths);
cids = spikes{1}.cids(iii);
cmap = jet(numel(spikes{1}.cids));
for ii = 1:numel(spikes{1}.cids)
    c = cids(ii);
    cnt = histc(spikes{1}.st(spikes{1}.clu==c), bins);
    plot(bins, smooth(100*cnt/max(cnt), 10) + depths(ii), 'Color', ([1 1 1] + cmap(ii,:))/2); hold on
end
%     plot(spikes{1}.st, spikes{1}.spikeDepths, '.k'); hold on

ts = arrayfun(@(x) x.start, stim);
te = arrayfun(@(x) x.stop, stim);
xx = [ts(:) te(:) te(:) ts(:)];

labels = {stim.label};
uniqueStims = unique(labels);

yy = kron(get(gca,'YLim'),[1,1]);

nL = numel(uniqueStims);
cmap = lines;
lh = [];
for iL = 1:nL
    inds = find(strcmp(uniqueStims{iL}, labels));
    for i = inds(:)'
        fhs = fill(xx(i,:)',yy(:),cmap(iL,:));
        set(fhs,'FaceAlpha',0.1,'LineStyle','none');
    end
    lh = [lh fhs];
end
xlabel('Time (seconds)')
ylabel('Cortical Depth (mm)')
legend(lh, uniqueStims, 'Location', 'BestOutside', 'Box', 'off')
set(gca, 'Box', 'off', 'Tickdir', 'out')
offsetAxes(gca)
ax = gca;