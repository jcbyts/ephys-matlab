function [h,el] = plotellipse(mu, covmat, r, varargin);
% h = plotellipse(mu, covmat, stdev, varargin);
%
% Plots an ellipse of constant probability (i.e. a contour of
% constant deviation stdev) for a given bivariate gaussian distr.
%
% Inputs:
%    mu = column vector with the mean of the distr.
%    covmat = 2x2 covariance matrix
%    stdev = the standard deviation contour we wish to plot (e.g. 1, 2, .2, etc)
%    varargin = additional style arguments for line plot
%
%  Sample call:  (plots 1-standard deviation ellipse for Gaussian centered at 1,2):
%     h = plotellipse([1;2], [1 .9; .9 1], 1, 'b--', 'linewidth', 2);
%
% Output: 
%    h = handle to plotted line
%    el = 100x2 matrix where the two columns provide the x and y values of the 
%         ellipse

[U, S, V] = svd(covmat);  % singular value decomposition of covariance
thet = [0:(2*pi)/99:(2*pi+.0001)]; %  angular variable theta

% make a circle, stretch by the sqrt of the eigenvalues, then
% rotate to align with eigenvectors
el = [r*U*sqrt(S)*[cos(thet); sin(thet)]]'; 

% add the mean offset for x and y
el(:,1) = el(:,1)+mu(1);  % x component
el(:,2) = el(:,2)+mu(2);  % y component

% Plot it
h = plot(el(:,1), el(:,2), varargin{:}); 
