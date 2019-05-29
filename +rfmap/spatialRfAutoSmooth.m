function [w,dc,params] = spatialRfAutoSmooth(X, y, lags, dims, varargin)
% Fit spatial receptive field with smoothness penalty
% 
% Inputs:
%   X [nFrames x nDims]    - matrix of the stimulus
%   y [nFrames x nNeurons] - spike counts
%   lags [1 x m]           - temporal lags to average over (e.g., 1:3)
%   dims [1 x n]           - [y, x] dimensions of the stimulus
%                            [1 x 1] scalar if 1D
%                            [1 x 2] twoplet if 2D
%                            [1 x 3] triplet if 3D (y, x, t)
%
% 
% Outputs:
%   B  = [nDims x nNeurons] - weights for each neuron
%   dc = [1 x nNeurons]     - dc term
%
% 08/10/2018 jly    wrote it

nobs   = size(X,1);
npred  = size(X,2);
nneuron = size(y,2);


% average over temporal lags
indices = bsxfun(@plus, (1:nobs)', lags);
indices((any(indices > nobs,2)),:) = [];
X = X(1:size(indices,1),:);

nobs = size(X,1); % new stimulus length
Ytmp = zeros(nobs, nneuron);
for i = 1:nneuron
    tmp = y(:,i);
    tmp = tmp(indices); % index into y to get lags
    Ytmp(:,i) = mean(tmp,2); % average over lags
end
y = Ytmp;


ip = inputParser();
ip.addOptional('crossvalidation', false)
ip.addOptional('lambda0', .1)
ip.addOptional('tol', 1e-4)
ip.addOptional('display', 'off')
ip.addOptional('addDC', true)
ip.addOptional('nFolds', 3)
ip.addOptional('spacetimesmoothratio', 1)
ip.parse(varargin{:})

nd = numel(dims);
switch nd
    case 1
        Cinv = qfsmooth1D(dims(1));
    case 2
        Cinv = qfsmooth2D(dims(2), dims(1));
    case 3
        Cinv = qfsmooth3D(permute(dims, [3, 2, 1]), [ip.Results.spacetimesmoothratio 1]);
    otherwise
        error('spatialRfAutoSmooth: number of dimensions must be 1 , 2, or 3')
end
Cinv = Cinv + 5*eye(size(Cinv,1)); % add ridge parameter

constCols = var(X)==0;
Cinv(constCols,:) = [];
Cinv(:,constCols) = [];
X(:,constCols) = [];

% if including a bias term, augment with a column of ones
if ip.Results.addDC
    X = [X ones(nobs,1)];
%     npred = npred + 1;
    Cinv = blkdiag(Cinv, .1); 
end


% use cross validation to estimate the hyper parameters
if ip.Results.crossvalidation
    
    nFolds = ip.Results.nFolds;
    folds = xvalidationIdx(nobs, nFolds, true);
    nFolds = 1;
    best_lambda = zeros(nFolds, nneuron);
    r2_fold     = zeros(nFolds, nneuron);
    
    
    for kFold = 1:nFolds
        
        test_ix  = folds{kFold,2};
        train_ix = folds{kFold,1};
        
        % sample covariance of the stimulus        
        XX = X(train_ix,:)'*X(train_ix,:);        
        
        % loop over neurons
        for i = 1:nneuron
            fprintf('Neuron %d / %d\n', i, nneuron)
            
            XY = X(train_ix,:)'*y(train_ix,i);
            w0 = XX \ XY; % linear regression initial
            
            testfun= @(w) rsquared(y(test_ix, i), X(test_ix,:)*w);
            
            r20 = testfun(w0); % rsquared with initial weights
            
            lambda0 = ip.Results.lambda0; % initialize the smoothness penalty
            
            % maximum a posteriori weights with smoothness
            w1  = (XX + lambda0 * Cinv) \ XY;
            r21 = testfun(w1);
            
            fprintf('Initial difference: %02.2f\n', r21-r20)
            iter = 0;
            maxIter = 100;
            
            while ((r21 - r20) > ip.Results.tol) && (iter < maxIter)
                iter = iter + 1;
                
                r20 = r21;
                
                lambda0 = lambda0 * 2; % double lambda
                w1  = (XX + lambda0 * Cinv) \ XY;
                
                r21 = testfun(w1);
                fprintf('Iter %d difference: %02.2f\n', iter, 1e3*(r21-r20))
            end
            
            best_lambda(kFold, i) = lambda0/2; % previous lambda
            r2_fold(kFold, i) = r20;
        end
        
    end
    
    if nFolds > 1
    lambdaMax = median(best_lambda);
    else
        lambdaMax = best_lambda;
    end
else
    lambdaMax = ip.Results.lambda0*ones(1,nneuron);
end

% --- Do final fitting with selected hyperparameters, estimate noise
%     variance
XX = X'*X;
XY = X'*y;
ny = size(y,1);


nsevar = zeros(1,nneuron);
k = size(X,2);
B = zeros(k, nneuron);
se = zeros(k, nneuron);
for i = 1:nneuron
    yy = y(:,i)'*y(:,i);
    % compute posterior weights
    B(:,i) = (XX + lambdaMax(i)*Cinv)\ XY(:,i); 
    % compute estimate of the noise variance
    nsevar(i) = (yy - 2*B(:,i)'*XY(:,i) + B(:,i)'*XX*B(:,i))/ny;
    % compute variance-covariance matrix for p(Bhat|B)
    BB = inv(XX + lambdaMax(i)*Cinv) * nsevar(i); %#ok<MINV>
    se(:,i) = diag(BB);
end

if any(constCols)
    w = nan(npred+ip.Results.addDC, nneuron);
    se2 = nan(size(se));
    ix = [~constCols(:); ip.Results.addDC];
    for i = 1:nneuron
        w(ix,i) = B(:,i);
        se2(ix,i) = se(:,i);
    end
else
    w = B;
    se2 = se;
end
        
% save out some details
params = ip.Results;
params.nsevar = nsevar;
params.se = se2;
params.w  = w;
params.lambdaMax = lambdaMax;
if ip.Results.crossvalidation
    params.fold_lambda = best_lambda;
    params.fold_r2     = r2_fold;
end

if ip.Results.addDC
    dc = B(end,:);
    B(end,:) = [];
else
    dc = [];
end

function r2 = rsquared(rtrue, rhat)
% r2 = rsquared(rtrue, rhat)
r2 = 1-(sum((rtrue(:)-rhat(:)).^2))/(sum((rtrue(:)-mean(rtrue(:))).^2));

function xidxs = xvalidationIdx(nTotalSamples, kxValidation, isRandomized, isPartition)
% xidxs = xvalidationIdx(nTotalSamples, kxValidation, isRandomized, isPartition)
% Get k-fold cross-validation indices.
%
% Input
%   nTotalSamples: integer number of samples or indices
%   kxValidation: k for the k-cross validation (if 1, training = test)
%   isRandomized: (opt) sequencial or randomized (default) indices?
%   isPartition: (opt) drop extra points at the end to make equal size partitions
% Output
%   xidxs: {kxValidation x 2} training, test indices
%
% Caution: test set is not guarantteed 
%          to be a partition, UNLESS mod(nTotalSamples, kxValidation) == 0
%          In such case, the test set is a partition of nTotalSamples.
%	   Use isPartition to enforce partitioning.
%

if nargin < 3
    isRandomized = true;
end

if nargin < 4
    isPartition = false;
end

if isPartition
    nTotalSamples = nTotalSamples - rem(nTotalSamples, kxValidation);
end

if rem(nTotalSamples,1) ~= 0 || rem(kxValidation,1) ~= 0
    error('xvalidationIdx:arguments should be both integers');
end

if ~numel(nTotalSamples) == 1
    ridx = nTotalSamples;
    nTotalSamples = numel(ridx);
else
    ridx = 1:nTotalSamples;
end

% size of the test set
m = ceil(nTotalSamples / kxValidation);
if m < 1
    error('xvalidationIdx:insufficient_samples', 'Not enough samples (%d) to create %d-fold cross-validation', nTotalSamples, kxValidation);
end

if isRandomized
    sidx = randperm(nTotalSamples);
    ridx = ridx(sidx);
end

if kxValidation == 1
    % No cross-validation. All samples are training, and all samples are test.
    xidxs{1,1} = ridx;
    xidxs{1,2} = ridx;
    return
end

startIdx = ceil(linspace(1, nTotalSamples, kxValidation+1));
for k = 1:kxValidation
    testSet = startIdx(k) + (1:m) - 1;
    trainSet = setdiff(1:nTotalSamples, testSet);
    xidxs{k,1} = ridx(trainSet); %#ok<AGROW>
    xidxs{k,2} = ridx(testSet); %#ok<AGROW>
end

% --- These functions create the smoothness penalties
function [qf] = qfsmooth1D(numx)
    %[qf] = qfsmooth1D(numx)
    %Create a quadratic form for smoothness regularization based on
    %second-order derivative operator, for a one-dimensional signal
    D = zeros(numx+1,numx);
    for ii = 1:numx-1
        D(ii,ii:ii+1) = [-1 1];
    end
    
    qf = D'*D;

function [qf,D,qfx,qfy] = qfsmooth2D(numx,numy)
    %[qf] = qfsmooth(numx,numy)
    %Create a quadratic form for smoothness regularization, with D a
    %first-order partial derivative operator and qf = D'*D
    D = zeros((numx-1)*numy + (numy-1)*numx,numx*numy);
    for jj = 1:numy
        for ii = 1:numx-1
            [xi, yi] = meshgrid(1:numx,1:numy);
            dd = (xi == ii & yi == jj) - (xi == ii + 1 & yi == jj);
            D(ii + (jj-1)*(numx-1),:) = dd(:);
        end
    end
    
    for jj = 1:numy-1
        for ii = 1:numx
            [xi, yi] = meshgrid(1:numx,1:numy);
            dd = (xi == ii & yi == jj) - (xi == ii & yi == jj + 1);
            D((numx -1)*numy + ii + (jj-1)*numx,:) = dd(:);
        end
    end
    
    qf  = D'*D;
    D1  = D(1:end/2,:);
    qfx = D1'*D1;
    D1  = D(end/2+1:end,:);
    qfy = D1'*D1;
    
function qf = qfsmooth3D(dims, lambda)
% qf = qfsmooth3D(dims, lambda)
    
    if nargin < 2
        lambda = 1;
    end
    
    numt = dims(1);
    numx = dims(2);
    numy = dims(3);
    
    if numel(lambda)==2
        lambda_space = lambda(2);
        lambda_time = lambda(1);
    else
        lambda_space = lambda;
        lambda_time = lambda;
    end
        
    Dt = qfsmooth1D(numt);
    Dxy = qfsmooth(numx, numy);
    qf = kron(lambda_space*Dxy,lambda_time*Dt);

   