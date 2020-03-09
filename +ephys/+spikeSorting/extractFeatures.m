function b = extractFeatures(w)
% Extract features for spike sorting.
%   b = extractFeatures(w) extracts features for spike sorting from the
%   waveforms in w, which is a 3d array of size length(window) x #spikes x
%   #channels. The output b is a matrix of size #spikes x #features.
%
%   This implementation does PCA on the waveforms of each channel
%   separately and uses the first three principal components. Thus, we get
%   a total of 12 features.

[~, n, k] = size(w);
q = 3;                              % number of components per channel
w = bsxfun(@minus, w, mean(w, 2));  % center data
b = zeros(n, k * q);
for i = 1:k
    C = w(:, :, i) * w(:, :, i)';   % covariance matrix
    [V, ~] = eigs(C, q);            % first q eigenvectors
    b(:, (1:q) + q * (i - 1)) = w(:, :, i)' * V;
end
