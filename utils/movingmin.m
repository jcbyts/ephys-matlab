function y = movingmin(x, n, blksz, dim);
%MOVINGMIN Find the minimum in a sliding window.
%   Y = MOVINGMIN(X,K) finds the moving minimum (over N points) along the
%   first non-singleton dimension of X.
%
%   At each point we compute the minimum value over a sliding window
%   starting at the current point and looking forward over a total
%   of N points, e.g., for N = 3, 
%
%     X = [1 3 5 4 6 3 5 2 3 10]
%     Y = [1 3 4 3 3 2 2 2 3 10]
%
%   Y = MOVINGMIN(X,N,BLKSZ) computes the output in blocks of BLKSZ
%   samples. Specifying BLKSZ << LENGTH(X) is useful on memory limited
%   systems (MOVINGMIN uses a working matrix of N x BLKSZ elements). By
%   default, BLKSZ == LENGTH(X) - this yields the fastest execution.
%
%   For matrices and N-D arrays, Y = MOVINGMIN(X,N,[],DIM) or 
%   Y = MOVINGMIN(X,N,BLKSZ,DIM) operates along the dimension DIM.
%
% $Id: movingmin.m,v 1.3 2009-04-16 02:08:24 shaunc Exp $

narginchk(1,4);

if nargin < 2,
  n = [];
end

if nargin < 3,
  blksz = [];
end

if nargin < 4,
  dim = [];
end

if ~isempty(dim) && dim > ndims(x)
	error('The specified dimension exceeds the dimensions of X.')
end

% reshape x into the right dimension.
if isempty(dim),
	% by default we work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else,
	% make dim in the first (row) dimension
	perm = [dim,1:dim-1,dim+1:ndims(x)];
	x = permute(x,perm);
end

% verify that the block size is valid...
sz = size(x);
if isempty(blksz),
	blksz = sz(1); % sz(1) is the number of rows of x (default)
else,
	blksz = blksz(:);
end

% initialize y with the correct dimension
y = zeros(sz);

% call movmin (see below)
for i = 1:prod(sz(2:end)),
	y(:,i) = movmin(x(:,i),n,blksz);
end

% revert y to the original shape of x
if isempty(dim),
	y = shiftdim(y, -nshifts);
else,
	y = ipermute(y,perm);
end


%--------------------------------------------------------------------------
function y = movmin(x,n,blksz),
% note: one dimensional moving minimum

nx = length(x);

X = [x; x(nx)*ones(n-1,1)];
y = zeros(nx,1);

% work in blocks to save memory
indr = (0:n-1)';
indc = 1:nx;
for i = 1:blksz:nx,
  ind = indc(ones(1,n),i:min(i+blksz-1,nx)) + ...
        indr(:,ones(1,min(i+blksz-1,nx)-i+1));
  xx = reshape(X(ind),n,min(i+blksz-1,nx)-i+1);
  y(i:min(i+blksz-1,nx)) = min(xx);
end
