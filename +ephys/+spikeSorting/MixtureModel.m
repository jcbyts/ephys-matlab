% Mixture model (Gaussian or Student's t).
%
% Selection criteria for split & merge canditates are based on Ueda et al.
% (2000), Neural Computation. We perform splits and merges separately and
% accept them if they lead to an improvement in BIC. This approach is
% inspired by Klustakwik (http://klustakwik.sourceforge.net).
%
% AE 2012-04-26

classdef MixtureModel
    properties
        D           % number of dimensions
        K           % number of clusters
        mu          % cluster means
        Sigma       % cluster covariances
        pi          % cluster priors
        df = Inf    % degrees of freedom (Inf  = Gaussian)
    end

    methods
        
        function self = MixtureModel(mu, Sigma, pi, df)
            % MixtureModel constructor
            %   model = MixtureModel(mu, Sigma, pi) creates a Gaussian
            %   mixture model with means mu, covariance matrices Sigma, and
            %   component priors pi.
            %
            %   model = MixtureModel(mu, Sigma, pi, df) creates a mixture
            %   of Student's t distributions with means mu, covariance
            %   matrices Sigma, component priors pi and df degrees of
            %   freedom.
            %
            %   Let D = #dimensions, k = #components:
            %       mu      K-by-D
            %       Sigma   D-by-D-by-K
            %       pi      1-byK
            %       df      scalar (Inf or ommitted if Gaussian model)
            [self.K, self.D] = size(mu);
            self.mu = mu;
            self.Sigma = Sigma;
            self.pi = pi;
            if nargin > 3
                self.df = df;
            end
        end
        
        
        function p = pdf(self, x)
            % Evaluate probability density function.
            p = sum(self.likelihood(x), 2);
        end
        
        
        function like = likelihood(self, x)
            % Likelihood of the data for each cluster.
            like = zeros(size(x, 1), self.K);
            for k = 1 : self.K
                if isinf(self.df)
                    like(:, k) = self.pi(k) * mvn(x, self.mu(k, :), self.Sigma(:, :, k));
                else
                    like(:, k) = self.pi(k) * mvt(x, self.mu(k, :), self.Sigma(:, :, k), self.df);
                end
            end
        end
        
        
        function post = posterior(self, x)
            % Posterior for class membership of each datapoint.
            likelihood = self.likelihood(x);
            p = sum(likelihood, 2);
            post = bsxfun(@rdivide, likelihood, p);
            post(p == 0, :) = 0;
        end
        
        
        function logLike = logLikelihood(self, x)
            % Log likelihood of the data under the model.
            logLike = sum(mylog(self.pdf(x)));
        end
        
        
        function bic = BIC(self, x)
            % Bayesian information criterion (BIC) for the model.
            N = size(x, 1);
            params = self.K * (0.5 * self.D * (self.D + 1) + self.D) + self.K;
            bic = -2 * self.logLikelihood(x) + params * log(N);
        end
        
        
        function assignment = cluster(self, x)
            % Cluster data using maximum a-posteriori
            [~, assignment] = max(self.posterior(x), [], 2);
        end
        
        
        function plotData(self, x, d1, d2)
            % Plot scatterplot of two selected data dimensions
            if nargin < 3, d1 = 1; d2 = 2; end
            assignment = self.cluster(x);
            c = jet(self.K);
%             figure(22), clf, 
            cla
            hold all
            hdl = zeros(1, self.K);
            for k = 1 : self.K
                plot(x(assignment == k, d1), x(assignment == k, d2), '.', 'markersize', 1, 'color', c(k, :))
            end
            for k = 1 : self.K
                [V, Lambda] = eig(self.Sigma([d1 d2], [d1 d2], k));
                for i = 1 : 2
                    a = sqrt(Lambda(i, i)) * V(:, i) * [-2 2];
                    plot(self.mu(k, d1) + a(1, :), self.mu(k, d2) + a(2, :), 'color', c(k, :))
                end
                plot(self.mu(k, d1), self.mu(k, d2), '.w', 'markersize', 40);
                hdl(k) = plot(self.mu(k, d1), self.mu(k, d2), '.', 'color', c(k, :), 'markersize', 25);
            end
            legend(hdl, arrayfun(@(x) sprintf('Cluster %d', x), 1 : self.K, 'UniformOutput', false))
            axis equal
            drawnow
        end
        
    end
    
    
    methods(Static)
        
        function [model, assignment] = fit(x, varargin)
            % Fit mixture model using EM and split & merge.
            %   model = MixtureModel.fit(x) fits a Gaussian mixture model.
            %
            %   model = MixtureModel.fit(x, df) fits a mixture of Student's
            %   t model with df degrees of freedom.
            %
            %   model = MixtureModel.fit(..., 'name', value) sets optional
            %   parameters (as name/value pairs) for fitting:
            %       verbose     verbose output
            %       nmax        maximum number of data points to use.
            %       tol         mininum relative change in likelihood to
            %                   determine convergence of the EM
            %       seed        initial seed for random number generator.
            %       ridge       level of ridge on covariance matrices (used
            %                   for regularization) in units of standard
            %                   deviations in muV.
            
            % catch degrees of freedom as second input
            if nargin > 1 && isnumeric(varargin{1})
                dof = varargin{1};
                varargin(1) = [];
            else
                dof = Inf;
            end
            
            % parse optional inputs
            p = inputParser;
            p.addOptional('verbose', false, @(x) isscalar(x) && islogical(x));
            p.addOptional('nmax', 10000, @(x) isscalar(x) && isnumeric(x));
            p.addOptional('tol', 0.005, @(x) isscalar(x) && isnumeric(x));
            p.addOptional('seed', 1, @(x) isscalar(x) && isnumeric(x));
            p.addOptional('ridge', 5, @(x) isscalar(x) && isnumeric(x) && x >= 0);
            p.parse(varargin{:});
            par = p.Results;
            
            % fix random number generator seed to get deterministic results
            if ~isnan(par.seed)
                rng(par.seed)
            end
            
            % use a subset of the data for training
            N = size(x, 1);
            Ntrain = min(par.nmax, N);
            rnd = randperm(N);
            xtrain = x(rnd(1:Ntrain), :);

            % initialize model
            model = ephys.spikeSorting.MixtureModel(mean(xtrain, 1), cov(xtrain), 1, dof);
            logLike = model.logLikelihood(xtrain);
            fprintf('Fitting mixture of Gaussians using split & merge\n')
            fprintf('Starting with one component')
            [model, logLike] = model.EM(xtrain, logLike, par);
            fprintf(' done.\n')
            if par.verbose, model.plotData(xtrain), end
            
            % Run split & merge
            % We alternate between trying to split and trying to merge.
            % Both are done until no candidate leads to success. If both
            % splitting and merging didn't lead to success we terminate
            sm = 1;
            success = true(1, 2);
            while any(success)
                switch sm
                    case 1
                        [model, logLike, success(sm)] = model.trySplit(xtrain, logLike, par);
                    case 2
                        [model, logLike, success(sm)] = model.tryMerge(xtrain, logLike, par);
                end
                if success(sm)
                    success(:) = true;
                else
                    sm = 3 - sm;
                end
            end
            
            % Finalize fit on the entire dataset
            fprintf('Done with split & merge. Final model fit on full datatset')
            model = model.EM(x, logLike, par);
            fprintf(' done.\n')
            
            % determine cluster assignments
            if nargout > 1
                assignment = model.cluster(x);
            end
        end
        
    end
    
    methods(Access = public)
%     methods(Access = private)
        
        function [self, logLike, success] = tryMerge(self, x, logLike, par)
            success = false;
            cands = self.getMergeCandidates(x);
            bicBase = self.BIC(x);
            for ij = cands'
                try
                    fprintf('Trying to merge clusters %d and %d', ij(1), ij(2))
                    newModel = self.mergeClusters(ij(1), ij(2));
                    [newModel, logLike] = newModel.EM(x, logLike, par);
                    bicNew = newModel.BIC(x);
                    if bicNew < bicBase
                        fprintf(' success (BIC decreased by %.5g)\n', bicBase - bicNew)
                        self = newModel;
                        success = true;
                        if par.verbose, self.plotData(x), end
                        break
                    else
                        fprintf(' aborted\n')
                    end
                catch err
                    if strcmp(err.identifier, 'MoGsm:starvation')
                        fprintf(' aborted due to component starvation\n')
                    else
                        rethrow(err)
                    end
                end
            end
        end
        
        
        function [self, logLike, success] = trySplit(self, x, logLike, par)
            success = false;
            cands = self.getSplitCandidates(x);
            bicBase = self.BIC(x);
            for i = cands'
                try
                    fprintf('Trying to split cluster %d', i)
                    newModel = self.splitCluster(i, x, par);
                    [newModel, logLike] = newModel.EM(x, logLike, par);
                    bicNew = newModel.BIC(x);
                    if bicNew < bicBase
                        fprintf(' success (BIC decreased by %.5g)\n', bicBase - bicNew)
                        self = newModel;
                        success = true;
                        if par.verbose, self.plotData(x), end
                        break
                    else
                        fprintf(' aborted\n')
                    end
                catch err
                    if strcmp(err.identifier, 'MoGsm:starvation')
                        fprintf(' aborted due to component starvation\n')
                    else
                        rethrow(err)
                    end
                end
            end
        end
        
        
        function [self, logLike] = EM(self, x, logLike, par)
            % Perform EM steps until convergence
            
            N = size(x, 1);
            iter = 1;
            ridge = eye(self.D) * par.ridge ^ 2;
            while iter < 2 || diff(logLike(end - 1 : end)) / delta > par.tol %#ok
                
                % Do a couple of iterations
                for i = 1 : 10
                    % E step
                    posterior = self.posterior(x);
                    
                    % M step
                    Nk = sum(posterior, 1);
                    for k = 1 : self.K
                        self.mu(k, :) = sum(bsxfun(@times, x, posterior(:, k)), 1) / Nk(k);
                        xmu = bsxfun(@minus, x, self.mu(k, :));
                        self.Sigma(:, :, k) = ridge + (xmu' * bsxfun(@times, xmu, posterior(:, k))) / Nk(k);
                    end
                    self.pi = Nk / sum(Nk);
%                     assert(all(self.pi * N > 2 * self.D) && all(self.pi > 0.005), 'MoGsm:starvation', ...
%                         'Component starvation: cluster %d', find(self.pi * N < 2 * self.D, 1))
                end
                
                % evaluate log-likelihood
                logLike(end + 1) = self.logLikelihood(x); %#ok
                if par.verbose, self.plotData(x), end
                
                if iter == 1
                    delta = diff(logLike(end - 1 : end));
                end
                iter = iter + 1;
                fprintf('.')
            end
        end
        
        
        function self = splitCluster(self, k, x, par)
            % Split cluster k into two by running a partial EM.
            deltaMu  = randn(1, self.D) * chol(self.Sigma(:, :, k));
            pmu = [self.mu(k, :) + deltaMu; self.mu(k, :) - deltaMu];
            pSigma = repmat(det(self.Sigma(:, :, k)) ^ (1 / self.D) * eye(self.D), [1 1 2]);
            ppi = [0.5 0.5];
            x = x(self.cluster(x) == k, :);
            partial = ephys.spikeSorting.MixtureModel(pmu, pSigma, ppi, self.df);
            logLike = partial.logLikelihood(x);
            par.verbose = false;
            par.tol = min(par.tol, 0.01);
            partial = partial.EM(x, logLike, par);
            self.mu(k, :) = partial.mu(1, :);
            self.mu(self.K + 1, :) = partial.mu(2, :);
            self.Sigma(:, :, k) = partial.Sigma(:, :, 1);
            self.Sigma(:, :, self.K + 1) = partial.Sigma(:, :, 2);
            self.pi(k) = partial.pi(1) * self.pi(k);
            self.pi(self.K + 1) = partial.pi(2) * self.pi(k);
            self.K = self.K + 1;
        end
        
        function self = mergeClusters(self, i, j)
            % Merge clusters i and j.
            self.mu(i, :) = (self.pi(i) * self.mu(i, :) + self.pi(j) * self.mu(j, :)) / (self.pi(i) + self.pi(j));
            self.Sigma(:, :, i) = (self.pi(i) * self.Sigma(:, :, i) + self.pi(j) * self.Sigma(:, :, j)) / (self.pi(i) + self.pi(j));
            self.pi(i) = self.pi(i) + self.pi(j);
            self.mu(j, :) = [];
            self.Sigma(:, :, j) = [];
            self.pi(j) = [];
            self.K = self.K - 1;
        end
        
        function cand = getSplitCandidates(self, x)
            % Get list of split candidates.
            N = size(x, 1);
            pk = zeros(N, self.K);
            for k = 1 : self.K
                pk(:, k) = mvn(x, self.mu(k, :), self.Sigma(:, :, k));
            end
            posterior = self.posterior(x);
            fk = bsxfun(@rdivide, posterior, sum(posterior, 1));
            Jsplit = sum(fk .* (mylog(fk) - mylog(pk)), 1);
            [~, cand] = sort(Jsplit, 'descend');
            cand = cand(self.pi(cand) * N > 4 * self.D); % don't split small clusters
            cand = cand(:);
        end
        
        function cand = getMergeCandidates(self, x)
            % Get list of merge candidates.
            posterior = self.posterior(x);
            maxCandidates = ceil(self.K * sqrt(self.K) / 2);
            np = sqrt(sum(posterior .* posterior, 1));
            Jmerge = zeros(self.K * (self.K - 1) / 2, 1);
            cand = zeros(self.K * (self.K - 1) / 2, 2);
            k = 0;
            for i = 1 : self.K
                for j = i + 1 : self.K
                    k = k + 1;
                    Jmerge(k) = posterior(i, :) * posterior(j, :)' / (np(i) * np(j));
                    cand(k, :) = [i j];
                end
            end
            [~, order] = sort(Jmerge, 'descend');
            cand = cand(order(1:min(end, maxCandidates)), :);
        end
    end
end


function y = mylog(x)
% Natural logarithm excluding zeros.
y = reallog(x);
y(x == 0) = 0;
end


function p = mvn(x, mu, Sigma)
% Multivariate normal probability density
%   p = mvn(x, mu, Sigma) calculates the density of the multivariate normal
%   distribution with mean mu and covariance matrix Sigma at x. x is
%   assumed to be a row vector or a matrix of multiple row vectors, in
%   which case the result, p, is a column vector.
D = size(Sigma, 1);
[Ch, ~] = chol(Sigma);
xmu = bsxfun(@minus, x, mu);
const = (2 * pi) ^ (-D / 2);
p = const / prod(diag(Ch)) * exp(-1/2 * sum((xmu / Ch) .^ 2, 2));
end


function p = mvt(x, mu, Sigma, df)
% Multivariate Student's t probability density.
%   p = mvt(x, mu, Sigma) calculates the density of the multivariate t
%   distribution mean mu, covariance matrix Sigma, and df degrees of
%   freedom at x. x is assumed to be a row vector or a matrix of multiple
%   row vectors, in which case the result, p, is a column vector.
D = size(Sigma, 1);
[Ch, ~] = chol(Sigma);
xmu = bsxfun(@minus, x, mu);
delta = sum((xmu / Ch) .^ 2, 2);
p = exp(gammaln((df + D) / 2) - gammaln(df / 2) ...
    - ((df + D) / 2) .* log(1 + delta / df) ...
    - sum(log(diag(Ch))) - (D / 2) * log(df * pi));
end

