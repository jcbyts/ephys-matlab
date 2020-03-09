function X = keepMaxClusters(results, T, threshold)
% X = keepMaxClusters(results, T, threshold) keeps only those clusters that
%   have their largest waveform on the center channel. Spikes are assumed
%   to be sorted in groups of K channels, with K-1 channels overlap between
%   adjacent groups. At the end, any duplicate spikes are removed. During
%   this process we always keep the spike from the cluster with the larger
%   waveform.

[D, ~, K] = size(results(1).w);
center = (K + 1) / 2;

R = numel(results);
spikes = {}; %#ok<*AGROW>
clusters = {};
mag = [];
w = zeros(D, K, 0);
for i = 1 : R
    r = results(i);
    a = cluster(r.model, r.b);
    for j = 1 : numel(r.model.pi)
        m = permute(sum(mean(r.w(:, a == j, :), 2) .^ 2, 1), [3 2 1]);
        [mm, ndx] = max(m);
        if ndx == center || (i == 1 && ndx < center) || (i == R && ndx > center)
            spikes{end + 1} = r.s(a == j); 
            clusters{end + 1} = numel(spikes) * ones(size(spikes{end}));
            mag(end + 1) = mm;
            w(:, :, end + 1) = permute(mean(r.w(:, a == j, :), 2), [1 3 2]);
        end
    end
end
M = numel(spikes);
total = cellfun(@numel, spikes);

spikes = cat(1, spikes{:});
clusters = cat(1, clusters{:});


% order spikes in time
[spikes, order] = sort(spikes);
clusters = clusters(order);

% remove spikes of smaller size
N = numel(spikes);
keep = true(N, 1);
prev = 1;
refrac = 4; % 1/3 ms
for i = 2 : N
    if spikes(i) - spikes(prev) < refrac
        if mag(clusters(i)) < mag(clusters(prev))
            keep(i) = false;
        else
            keep(prev) = false;
            prev = i;
        end
    else
        prev = i;
    end
end
spikes = spikes(keep);
clusters = clusters(keep);

% remove clusters that lost too many spikes to other clusters
frac = hist(clusters, 1 : M) ./ total;
keep = true(numel(spikes), 1);
for i = 1 : M
    if frac(i) < threshold
        keep(clusters == i) = false;
    end
end
spikes = spikes(keep);
clusters = clusters(keep);
[~, ~, clusters] = unique(clusters);
spikes = round(spikes);
% create spike matrix
X = sparse(spikes, clusters, 1, T, max(clusters));
