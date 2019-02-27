function results = fitParametric2dRf(stim, RF)
% Auto Gaussian & Gabor Surface fit 
% --- 
% Functions to fit a 1D Gaussian to a curve and a 2D Gaussian or Gabor to a surface. The routines are automatic in the sense that they do not require the specification of starting guesses for the model parameters. This is done by evaluating the quality of fit for many different choices of parameters then refining the most promising set of params through least-squares (exhaustive search followed by refinement).
% 
% All functions support 2 methods for computing error bars on the parameters: bootstrapping and MCMC.
% 
% autoGaussianSurf(xi,yi,zi) fits a 2D Gaussian to a surface, defined as:
% 
% zi = a*exp(-((xi-x0).^2/2/sigmax^2 + (yi-y0).^2/2/sigmay^2)) + b
% 
% It can also fit a tilted 2d Gaussian or isotropic 2d Gaussian.
% 
% autoGaborSurf(xi,yi,zi) fits a Gabor, defined as:
% 
% zi = a*exp(-(xip,.^2+yip.^2)/2/sigma^2)*cos(2*pi*xip/lambda + phase) + b
% 
% Where: 
% xip = (xi-x0)*cos(theta) + (yi-i0)*sin(theta); 
% yip =-(xi-x0)*sin(theta) + (yi-i0)*cos(theta);
% 
% The Gabor fit calls autoGaussianSurf internally, using the fact that the absolute value of a Gabor in the Fourier domain is a Gaussian.
% 
% autoGaussianCurve(xi,zi) fits a 1D Gaussian to a curve.
% from https://www.mathworks.com/matlabcentral/fileexchange/31485-auto-gaussian-gabor-fits?s_tid=prof_contriblnk

[xi, yi] = meshgrid(stim.xax, stim.yax);
opts.errorbars = 'bootstrap';

nNeurons = size(RF,2);
results = repmat(struct('gabor', [], 'gaussian', []), nNeurons,1);
for i = 1:nNeurons
    zi = RF(:,i) / max(RF(:,i));
    zi = reshape(zi, stim.size);
    results(i).gaussian = autoGaussianSurf(xi,yi,zi,opts);
%     results(i).gabor    = autoGaborSurf(xi, yi,zi, opts);
end


function [results] = autoGaborSurf(xi,yi,zi,opts)
    %function [results] = autoGaborSurf(xi,yi,zi,opts)
    %
    %Fit a Gabor (a sinusoid windowed by a Gaussian) to a surface. 
    %The Gabor is defined as follows:
    %zi = a*exp(-(xip.^2+yip.^2)/2/sigma^2)*cos(2*pi*xip/lambda + phase) + b
    %
    % Where:
    % xip = (xi-x0)*cos(theta) + (yi-i0)*sin(theta);
    % yip =-(xi-x0)*sin(theta) + (yi-i0)*cos(theta);
    %
    %On the assumption of Gaussian noise through maximum likelihood
    %(least-squares).
    %
    %The procedure is automatic in the sense that it is unlikely to get 
    %stuck in local minima which makes it appropriate in fully automated 
    %scenarios. This is accomplished by a trick: the absolute value of a
    %Gabor in the Fourier domain is a Gaussian. Thus the procedure finds
    %a good starting set of values for the frequency and orientation of the
    %Gabor by calling autoGaussianSurf on the Fourier transforms of the
    %data, uses that to get a good start guess for x0 and y0, and refines
    %the params using lsqcurvefit.
    %
    %Currently only regular grids (generated through meshgrid) are accepted
    %for xi and yi. 
    %
    %opts is a struct containing options:
    % opts.errorbars (none | bootstrap | mcmc) - Type of errors bars
    % requested (default:none, see TryAutoGaussianSurf for use)
    %
    %Example use:
    %
    % [xi,yi] = meshgrid(-10:10,-20:20);
    % xip = (xi-4)*cos(pi/6) + yi*sin(pi/6);
    % yip =-(xi-4)*sin(pi/6) + yi*cos(pi/6);
    % zi = exp(-(xip.^2+yip.^2)/2/2^2).*cos(xip*2*pi/5+pi/3) + .2*randn(size(xi));
    % results = autoGaborSurf(xi,yi,zi);
    %
    %See also autoGaussianSurf
    if nargin < 4
        opts = struct();
    end
    
    p = inputParser;
    p.KeepUnmatched = true;
    p.addOptional('errorbars','none');
    p.addOptional('bsIters',500);
    p.addOptional('mcmcIters',5000);
    p.addOptional('precision',.001*ones(8,1));
    p.parse(opts);
    opts = p.Results;
    
    %Start with fitting a Gaussian in the Fourier domain
    zif = fftshift(abs(fft2(zi)));
    %smooth a little, this helps
    g = [.2,1,.2];
    g = g'*g;
    zif = conv2(zif-mean(zif(:)),g,'same');
    
    ss = @(x) x(1:end-1);
    
    %Because of wraparound, do it in two shots
    [xip,yip] = meshgrid(linspace(-.5,0,size(xi,2)/2),ss(linspace(-.5,.5,size(xi,1)+1)));
    r1 = autoGaussianSurf(xip,yip,zif(:,1:size(xi,2)/2));
    
    %
    [xip,yip] = meshgrid(ss(linspace(-.5,.5,size(xi,2)+1)),linspace(-.5,0,size(xi,1)/2));
    r2 = autoGaussianSurf(xip,yip,zif(1:size(xi,1)/2,:));
    
    %Look at the quality of each fit to decide which to use
    
    if r1.sse < r2.sse
        %Go with r1
        coords = [r1.x0,r1.y0,max(r1.sigmax,r1.sigmay)];
    else
        coords = [r2.x0,r2.y0,max(r2.sigmax,r2.sigmay)];
    end
    
    %0 degrees is vertical
    theta = atan2(coords(2),coords(1));
    sf = min(sqrt(coords(1).^2+coords(2).^2),.5);
    dx = xi(1,2)-xi(1,1);
    dy = yi(2,1)-yi(1,1);
    da = .5*(dx+dy);
    lambda = 1/sf*da;
    sigma  = da/coords(3)/pi/2/1.3;
    phase = 0;
    
    %Get Gabor eval'd with current grid
    g = evalGabor([0,0,theta,lambda,sigma,phase],[xi(:)-mean(xi(:)),yi(:)-mean(yi(:))]);
    
    s = conv2(zi,reshape(g,size(xi)),'same');
    [y0,x0] = find(abs(s)==max(abs(s(:))));
    x0 = xi(1,x0);
    y0 = yi(y0,1);
    
    minsigma = da*.2;
    maxsigma = (max(xi(:))-min(xi(:))+max(yi(:))-min(yi(:)))/2*.4;
    
    %Now use lsqcurvefit to finish
    lb = [-Inf,-Inf,min(xi(:)),min(yi(:)),-Inf,-Inf,minsigma /1.01,-Inf]';
    ub = [ Inf, Inf,max(xi(:)),max(yi(:)), Inf, Inf,maxsigma + .01, Inf]';
    
    p0 = [x0,y0,theta,lambda,sigma,phase];
    g = evalGabor(p0,[xi(:),yi(:)]);
    
    p = [g,ones(size(g))]\zi(:);
    
    params0 = [p(1),p(2),p0]';
    varnames = {'a','b','x0','y0','theta','lambda','sigma','phase'};
    iscircular = [0,0,0,0,1,0,0,1]'==1;
    
    results = doFinalOptimization(@evalGabor,[xi(:),yi(:)],zi(:),params0,lb,ub,true(length(lb),1),varnames,iscircular,opts);
    
    %Collect the results
    results.G = reshape(results.G,size(zi));



function g = evalGabor(ps,X)
    xi = X(:,1);
    yi = X(:,2);
    x0 = ps(1);
    y0 = ps(2);
    theta = ps(3);
    lambda = ps(4);
    sigma = ps(5);
    phase = ps(6);
    
    xip =  (xi-x0)*cos(theta) + (yi-y0)*sin(theta);
    yip = -(xi-x0)*sin(theta) + (yi-y0)*cos(theta);
    
    g = exp(-(xip.^2+yip.^2)/2/sigma^2).*cos(2*pi*xip/lambda+phase);



function [results] = autoGaussianSurf(xi,yi,zi,opts)
    %function [results] = autoGaussianSurf(xi,yi,zi,opts)
    %
    %Fit a surface zi = a*exp(-((xi-x0).^2/2/sigmax^2 + ...
    %                           (yi-y0).^2/2/sigmay^2)) + b
    %
    %On the assumption of Gaussian noise through maximum likelihood
    %(least-squares).
    %
    %The procedure is "robust" in the sense that it is unlikely to get 
    %stuck in local minima which makes it appropriate in fully automated 
    %scenarios. This is accomplished through an initial exhaustive search 
    %for the parameters, followed by refinement with lsqcurvefit
    %
    %Currently only regular grids (generated through meshgrid) are accepted
    %for xi and yi. 
    %Example use:
    %
    %[xi,yi] = meshgrid(-10:10,-20:20);
    %zi = exp(-(xi-3).^2-(yi+4).^2) + randn(size(xi));
    %results = autoGaussianSurf(xi,yi,zi)
    %
    %Returns results, a struct with elements a,b,x0,y0,sigmax,sigmay for
    %the estimated model parameters, G which is the Gaussian evaluated with
    %the estimated model params, and sse, sse0 and r2, which are the sum of
    %squared error for the fitted Gaussian, the sum of squared error for a
    %model with only an offset, and the R^2 proportion of variance
    %accounted for.
    %
    %opts is a struct containing options:
    % opts.errorbars (none | bootstrap | mcmc) - Type of errors bars
    % requested (default:none)
    % opts.iso (true|false) - If true, Gaussian is isotropic (sigmax ==
    % sigmay, default false)
    % opts.tilted (true|false) - If true, Gaussian is tilted by an angle,
    % default false). In that case the definition of the Gaussian becomes:
    % 
    % xip = (xi-x0)*cos(theta) + (yi-i0)*sin(theta);
    % yip =-(xi-x0)*sin(theta) + (yi-i0)*cos(theta);
    % zi = a*exp(-(xip.^2/2/sigmax^2 + ...
    %              yip.^2/2/sigmay^2)) + b
    %
    % opts.positive: if true, gain (a) is required to be > 0

    if nargin < 4
        opts = struct();
    end
    
    %parse inputs
    p = inputParser;
    p.KeepUnmatched = true;
    p.addOptional('positive',false);
    p.addOptional('errorbars','none');
    p.addOptional('bsIters',500);
    p.addOptional('mcmcIters',5000);
    p.addOptional('precision',[.001,.001,.001,.001,.001,.001,.001]');
    p.addOptional('iso',false);
    p.addOptional('tilted',false);
    p.parse(opts);
    opts = p.Results;

    sz = size(zi);
    
    %Verify that the grid is regular
    if any(any(abs(diff(xi,2,2)) >=1e3*eps)) || any(any(diff(yi,2,1) >= 1e3*eps))
        error('xi or yi is not a regular grid');
    end
    
    if any(size(zi)~=size(xi)) || any(size(zi)~=size(yi))
        error('xi, yi and zi are not the same size');
    end
    
    if opts.tilted && opts.iso
        error('A Gaussian cannot be both isotropic and tilted');
    end
    
    xi = xi(:);
    yi = yi(:);
    boundx = [min(xi),max(xi)];
    boundy = [min(yi),max(yi)];
    
    %Find a minimum sigma based on number of elements, range of x and y
    rgx = diff(boundx);
    minsigmax = rgx/sz(2)/5;
    maxsigmax = rgx/2;

    rgy = diff(boundy);
    minsigmay = rgy/sz(1)/5;
    maxsigmay = rgy/2;
    
    minsigma = min(minsigmax,minsigmay);
    maxsigma = max(maxsigmax,maxsigmay);
    sigmas = exp(log(minsigma):.3:log(maxsigma));
    
    rgx = [0:sz(2)/2,-ceil(sz(2)/2)+1:-1]';
    rgy = [0:sz(1)/2,-ceil(sz(1)/2)+1:-1]';
    
    res = zeros(length(sigmas),7);
    
    %Run through all the different values for sigma
    for ii = 1:length(sigmas)
        thefiltx = exp(-rgx.^2/2/sigmas(ii));
        thefilty = exp(-rgy.^2/2/sigmas(ii));
        %convolve zi with gaussian filters and find the maximum response
        %(matched filtering)
        zi2 = reflectconv(reflectconv(zi,thefilty)',thefiltx)';
        if opts.positive
            [~,pos] = max(zi2(:));
        else
            [~,pos] = max(abs(zi2(:)));
        end
        x0 = xi(pos);
        y0 = yi(pos);
        %[y0,x0] = ind2sub(sz,pos);
        
        %Determine the residual error for the optimal x, y for this sigma
        G = exp(-((xi-x0).^2+(yi-y0).^2)/2/sigmas(ii)^2);
        X = [G,ones(length(G),1)];
        ps = X\zi(:);
        res(ii,:) = [sum((zi(:) - X*ps).^2),ps(:)',x0,y0,sigmas(ii),sigmas(ii)];
    end
    
    %Find sigma with corresponding least error
    [~,optsigma] = min(res(:,1));
    
    %Fit the parameters again through lsqcurvefit
    if opts.iso
        lb = [-Inf,-Inf,boundx(1),boundy(1),minsigmax /1.01]';
        ub = [ Inf, Inf,boundx(2),boundy(2),maxsigmax + .01]';
        params0 = res(optsigma,2:end-1)';
        varnames = {'a','b','x0','y0','sigma'};
        iscircular = false(5,1);
        thefun = @pointisogaussian;
        opts.precision = opts.precision(1:5);
    elseif ~opts.tilted
        lb = [-Inf,-Inf,boundx(1),boundy(1),minsigmax /1.01,minsigmay /1.01]';
        ub = [ Inf, Inf,boundx(2),boundy(2),maxsigmax + .01,maxsigmay + .01]';
        params0 = res(optsigma,2:end)';
        varnames = {'a','b','x0','y0','sigmax','sigmay'};
        iscircular = false(6,1);
        thefun = @pointgaussian;
        opts.precision = opts.precision(1:6);
    else
        %Fit a Gaussian to the power spectrum of the thing
        theta = getInitialAngle(reshape(xi,sz),reshape(yi,sz),reshape(zi,sz),opts);
        
        %Because of wraparound, do it in two shots
        lb = [-Inf,-Inf,boundx(1),boundy(1),minsigmax /1.01,minsigmay /1.01,-Inf]';
        ub = [ Inf, Inf,boundx(2),boundy(2),maxsigmax + .01,maxsigmay + .01,Inf]';
        params0 = [res(optsigma,2:end-2)';res(optsigma,end-1:end)'.*[1.3,0.7]';theta];
        varnames = {'a','b','x0','y0','sigmax','sigmay','angle'};
        iscircular = [false(6,1);true];
        thefun = @pointtiltedgaussian;
    end
    
    if opts.positive
        lb(1) = 0;
    end
    zi(isnan(zi)) = 0;
    results = doFinalOptimization(thefun,[xi(:),yi(:)],zi(:),params0,lb,ub,true(length(lb),1),varnames,iscircular,opts);
    
    %Collect the results
    results.G = reshape(results.G,size(zi));


function theta = getInitialAngle(xi,yi,zi,opts)
    zif = fftshift(abs(fft2(zi)));
    %smooth a little, this helps
    g = [.2,1,.2];
    g = g'*g;
    zif = conv2(zif-mean(zif(:)),g,'same');

    ss = @(x) x(1:end-1);
    
    opts.tilted = false;
    opts.errorbars = 'none';
    [xip,yip] = meshgrid(linspace(-.5,0,size(xi,2)/2),ss(linspace(-.5,.5,size(xi,1)+1)));
    r1 = autoGaussianSurf(xip,yip,zif(:,1:size(xi,2)/2),opts);

    %
    [xip,yip] = meshgrid(ss(linspace(-.5,.5,size(xi,2)+1)),linspace(-.5,0,size(xi,1)/2));
    r2 = autoGaussianSurf(xip,yip,zif(1:size(xi,1)/2,:),opts);

    %Look at the quality of each fit to decide which to use

    if r1.sse < r2.sse
        %Go with r1
        coords = [r1.x0,r1.y0,max(r1.sigmax,r1.sigmay)];
    else
        coords = [r2.x0,r2.y0,max(r2.sigmax,r2.sigmay)];
    end

    %0 degrees is vertical
    theta = atan2(coords(2),coords(1));



function [thef] = pointisogaussian(x,xdata)
    thef = exp(-.5*sum(bsxfun(@minus,xdata,[x(1),x(2)]).^2,2)*(1./x(3).^2));


function [thef] = pointtiltedgaussian(x,xdata)
    xdat = bsxfun(@minus,xdata,[x(1),x(2)]);
    xdat = xdat*[cos(x(5)),-sin(x(5)); sin(x(5)),cos(x(5))];
    thef = exp(-.5*(xdat.^2)*(1./[x(3);x(4)].^2));


function [thef] = pointgaussian(x,xdata)
    thef = exp(-.5*(bsxfun(@minus,xdata,[x(1),x(2)]).^2)*(1./[x(3);x(4)].^2));


%Convolution with reflected edge handling
function A = reflectconv(A,f)
    A = bsxfun(@times,fft([A(end:-1:1,:);A;A(end:-1:1,:)]),fft([f(1:floor(end/2));zeros(length(f)*2,1);f(floor(end/2)+1:end)]));
    A = ifft(A);
    A = A(length(f)+1:2*length(f),:);



function results = doFinalOptimization(curvefun,X,z,p0,lb,ub,mask,varnames,iscircular,opts)
    %Function internally called by auto* to perform the final optimization
    %or simulation based on initial parameters. This is what calls
    %lsqcurvefit, does bootstrapping or MCMC
    switch opts.errorbars
        case 'none'
            %Maximum-likelihood
            params = getMLestimate(curvefun,'Iter',p0,X,z,lb,ub,mask);
        case 'bootstrap'
%             if parpool('local') == 0
%                 warning('doFinalOptimization:noparallel',...
%                         'Parallel computing disabled. Try calling matlabpool open before using bootstrap estimates, much faster');
%             end
            I = ceil(rand(size(X,1),opts.bsIters)*size(X,1));
            parfor ii = 1:opts.bsIters
                params(:,ii) = getMLestimate(curvefun,'Off',p0,X(I(:,ii),:),z(I(:,ii)),lb,ub,mask);
            end
            
        case 'mcmc'
            if ~exist('dramrun.m','file')
                thepath = fileparts(which('doFinalOptimization'));
                addpath(genpath([thepath '/dram']));
                addpath(genpath([thepath '/econo']));
                if ~exist('dramrun.m','file')
                    error('Mcmc packages not installed. Run fetchMcmcPackages and try again');
                end
            end
            
            model.ssfun = @(p,z) sum((p(1)*curvefun(p(3:end),X)+p(2)-z).^2);
            model.priorfun = @(p,x) sum(opts.precision.*p(:).^2);
            
            G = curvefun(p0(3:end),X);
            sigma2 = mean((G-z).^2);
            
            simparams.par0 = p0;
            simparams.sigma2 = sigma2;
            simparams.n0 = 1;
            simparams.n  = numel(z);
            simparams.bounds = [lb';ub'];

            mcmcopts.nsimu = opts.mcmcIters;
            mcmcopts.qcov = eye(length(p0))*sigma2*100;
            mcmcopts.adaptint = 10;
            
            nmin = 5000;
            params = [];
            s2samples = [];
            while size(params,1) < nmin
                %Run dram
                [diagnostics,paramsn,s2samplesn] = dramrun(model,z,simparams,mcmcopts);
                params = [params;paramsn];
                s2samples = [s2samples;s2samplesn];
                
                %Run coda and look at convergence diagnostics
                r = coda(params);
                %Get the number of iterations coda says I should plus a
                %little padding
                nmin = max(cell2mat({r.n}));
                
                mcmcopts.nsimu = ceil(nmin-size(params,1)+500);
                mcmcopts.qcov = diagnostics.qcov;
                simparams.par0 = params(end,:)';
            end
            
            maxburn = max(cell2mat({r.nburn}));
            
            params = params(maxburn+1+500:end,:);
            s2samples = s2samples(maxburn+1+500:end,:);
            
            params = [params';sqrt(s2samples)'];
            varnames{end+1} = 'sigmanoise';
            results.diagnostics = diagnostics;
            iscircular = [iscircular;false];
        otherwise
            error('Unknown method of obtaining error bars: %s',opts.errorbars);
    end
    
    if size(params,2) == 1
        for ii = 1:length(varnames)
            results.(varnames{ii}) = params(ii);
        end
    else
        results = getSummaryStatistics(params',varnames,iscircular);
        results.samples = params;
        params = median(params');
        
        Gp = 0;
        for ii = 1:size(results.samples,2)
            p = results.samples(:,ii);
            Gp = Gp + curvefun(p(3:end),X)*p(1) + p(2);
        end
        results.Gp = Gp/size(results.samples,2);
    end
    results.G = curvefun(params(3:end),X)*params(1) + params(2);
    results.sse = sum((z-results.G).^2);
    results.sse0 = sum((z-mean(z)).^2);
    results.r2 = 1-results.sse/results.sse0;



function p = getMLestimate(curvefun,display,p0,X,z,lb,ub,mask)
    opts = optimset('Display',display,'Jacobian','on');
    params = lsqcurvefit(@(p,X) dodiff(curvefun,p,p0,X,mask),p0(mask),X,z,lb(mask),ub(mask),opts);
    p = p0;
    p(mask) = params;


function [eta,deta] = dodiff(curvefun,p,p0,x,mask)
%Do differentiation of curve/surf wrt parameters
    ps = p0;
    ps(mask) = p;
    p = ps;
    a = p(1);
    b = p(2);
    p = p(3:end);
    eta = curvefun(p,x);
    deta = zeros(length(eta),length(p)+2);
    delta = 1e-5;
    for ii = 1:length(p)
        dp = zeros(size(p));
        dp(ii) = delta;
        deta(:,ii+2) = (curvefun(p+dp,x)-eta)/delta;
    end
    deta(:,2) = 1;
    deta(:,1) = eta;
    deta = deta(:,mask);
    
function results = getSummaryStatistics(samples,names,iscircular)
    p = [0.025,.05,.25,.5,.75,.95,.975];
    for ii = 1:size(samples,2)
        if iscircular(ii) > 0
            x = mean(cos(samples(:,ii)));
            y = mean(sin(samples(:,ii)));
            
            r = sqrt(y.^2+x.^2);
            results.means.(names{ii}) = atan2(mean(y),mean(x));
            results.stds.(names{ii})  = sqrt(-2*log(r));
        else
            results.means.(names{ii}) = mean(samples(:,ii));
            results.stds.(names{ii})  =  std(samples(:,ii));
            results.quantiles.(names{ii}) = quantile(samples(:,ii),p);
        end
    end
    
    results.quantiles.key = p;

