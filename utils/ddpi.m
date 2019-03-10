function data = ddpi(fname)
%set(0,'DefaultFigureWindowStyle','docked')
% close all;
% clear all;

global P1_ROI_SIZE;
global P4_ROI_SIZE;
global THRESHOLD;
global IMG_WIDTH;
global IMG_HEIGHT;

P1_ROI_SIZE = 128;
P4_ROI_SIZE = 32;
THRESHOLD = 128;
IMG_WIDTH = 2040;
IMG_HEIGHT = 1088;

if nargin < 1
    fname = 'test.avi';
end

fprintf('Opening [%s]\n', fname)

[~, fn, ~] = fileparts(fname);
outfname = [fn 'tracersm.mat'];
if exist(outfname, 'file')
    tmp = load(outfname);
    data = tmp.data;
    return
end

v = VideoReader(fname);

tic
fprintf('Running tracking analysis...')
[p1x, p1y, p1offx, p1offy, p4x, p4y, p4offx, p4offy] = track(v);

fprintf('[%02.2fs]\n', toc)

data.x1 = p1x + p1offx;
data.x4 = p4x + p4offx;
data.y1 = p1y + p1offy;
data.y4 = p4y + p4offy;


save(outfname, 'data');
end

function [p1x, p1y, p1offx, p1offy, p4x, p4y, p4offx, p4offy] = track(video)
global P1_ROI_SIZE;
global P4_ROI_SIZE;
global THRESHOLD;

p1x = [];
p1y = [];
p1offx = [];
p1offy = [];

p4x = [];
p4y = [];
p4offx = [];
p4offy = [];

frame = 1;
while hasFrame(video)
    img = readFrame(video);
    img = gpuArray(img(:, :, 1));
    
    fimg = imgaussfilt(double(img), 1);
    [limg, nLabels] = bwlabel(fimg > THRESHOLD, 8);
    
    if (nLabels == 2)
        Bbox = findBlob(limg, nLabels);
        
        [sx, ex, sy, ey] = ROI(Bbox(1, :), P1_ROI_SIZE);
        p1Img = fimg(sy:ey, sx:ex);
        [p1x(frame), p1y(frame)] = radialcenter(p1Img);
%         [p1x(frame), p1y(frame)] = centerOfMass(p1Img, 200);
        p1offx(frame) = sx;
        p1offy(frame) = sy;
        
        [sx, ex, sy, ey] = ROI(Bbox(2, :), P4_ROI_SIZE);
        p4Img = fimg(sy:ey, sx:ex);
        [p4x(frame), p4y(frame)] = radialcenter(p4Img);
%         [p4x(frame), p4y(frame)] = centerOfMassWithResample(p4Img, 125, 10);
        p4offx(frame) = sx;
        p4offy(frame) = sy;
        
    else
        p1x(frame) = -1;
        p4x(frame) = -1;
        p1y(frame) = -1;
        p4y(frame) = -1;
        p1offx(frame) = -1;
        p1offy(frame) = -1;
        p4offx(frame) = -1;
        p4offy(frame) = -1;
    end
    frame = frame + 1;
end
end

function [sx, ex, sy, ey] = ROI(Bbox, roi_size)
global IMG_WIDTH;
global IMG_HEIGHT;

SIZE = int32(roi_size / 2);
bx = int32((Bbox(1) + Bbox(3)) / 2);
by = int32((Bbox(2) + Bbox(4)) / 2);
sx = max(1, bx - SIZE);
ex = min(IMG_WIDTH, bx + SIZE - 1);
sy = max(1, by - SIZE);
ey = min(IMG_HEIGHT, by + SIZE - 1);
end

function [cx, cy] = centerOfMass(img, threshold)
[h, w] = size(img);
img(img < threshold) = 0;
xIdx = repmat(gpuArray(double(1:w)), h, 1);
yIdx = repmat(gpuArray(double(1:h)'), 1, w);
s = gather(sum(img(:)));
sx = img .* xIdx;
sy = img .* yIdx;

cx = gather(sum(sx(:))) / s;
cy = gather(sum(sy(:))) / s;
end

function [cx, cy] = centerOfMassWithResample(img, threshold, RESAMPLE)
[h, w] = size(img);

F = griddedInterpolant(gather(double(img)));
xq = (0:1/RESAMPLE:w)';
yq = (0:1/RESAMPLE:h)';
img = F({xq, yq});

[h, w] = size(img);
img(img < threshold) = 0;
xIdx = repmat(double(1:w), h, 1);
yIdx = repmat(double(1:h)', 1, w);
s = sum(img(:));
sx = img .* xIdx;
sy = img .* yIdx;

cx = sum(sx(:)) / s  / RESAMPLE;
cy = sum(sy(:)) / s / RESAMPLE;
end

function Bbox = findBlob(limg, nLabels)
global IMG_HEIGHT;
Bbox = zeros(nLabels, 4);
Areas = zeros(nLabels, 1);
for i=1:nLabels
    idx = int32(find(limg == i));
    xIdx = unique(gather(idivide(idx, IMG_HEIGHT)));
    yIdx = unique(gather(mod(idx, IMG_HEIGHT)));
    
    Bbox(i, :) = [min(xIdx), min(yIdx), max(xIdx), max(yIdx)];
    width = Bbox(i, 3) - Bbox(i, 1);
    height = Bbox(i, 4) - Bbox(i, 2);
    Areas(i) = width * height;
end
[~, I] = sort(Areas, 'descend');
Bbox = Bbox(I, :);
end



function [xc yc sigma] = radialcenter(I)
I = gather(double(I));
% Number of grid points
[Ny Nx] = size(I);

% grid coordinates are -n:n, where Nx (or Ny) = 2*n+1
% grid midpoint coordinates are -n+0.5:n-0.5;
% The two lines below replace
%    xm = repmat(-(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5,Ny-1,1);
% and are faster (by a factor of >15 !)
% -- the idea is taken from the repmat source code
xm_onerow = -(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5;
xm = xm_onerow(ones(Ny-1, 1), :);
% similarly replacing
%    ym = repmat((-(Ny-1)/2.0+0.5:(Ny-1)/2.0-0.5)', 1, Nx-1);
ym_onecol = (-(Ny-1)/2.0+0.5:(Ny-1)/2.0-0.5)';  % Note that y increases "downward"
ym = ym_onecol(:,ones(Nx-1,1));

% Calculate derivatives along 45-degree shifted coordinates (u and v)
% Note that y increases "downward" (increasing row number) -- we'll deal
% with this when calculating "m" below.
dIdu = I(1:Ny-1,2:Nx)-I(2:Ny,1:Nx-1);
dIdv = I(1:Ny-1,1:Nx-1)-I(2:Ny,2:Nx);

% Smoothing --
h = ones(3)/9;  % simple 3x3 averaging filter
fdu = conv2(dIdu, h, 'same');
fdv = conv2(dIdv, h, 'same');
dImag2 = fdu.*fdu + fdv.*fdv; % gradient magnitude, squared

% Slope of the gradient .  Note that we need a 45 degree rotation of
% the u,v components to express the slope in the x-y coordinate system.
% The negative sign "flips" the array to account for y increasing
% "downward"
m = -(fdv + fdu) ./ (fdu-fdv);

% *Very* rarely, m might be NaN if (fdv + fdu) and (fdv - fdu) are both
% zero.  In this case, replace with the un-smoothed gradient.
NNanm = sum(isnan(m(:)));
if NNanm > 0
    unsmoothm = (dIdv + dIdu) ./ (dIdu-dIdv);
    m(isnan(m))=unsmoothm(isnan(m));
end
% If it's still NaN, replace with zero. (Very unlikely.)
NNanm = sum(isnan(m(:)));
if NNanm > 0
    m(isnan(m))=0;
end

%
% Almost as rarely, an element of m can be infinite if the smoothed u and v
% derivatives are identical.  To avoid NaNs later, replace these with some
% large number -- 10x the largest non-infinite slope.  The sign of the
% infinity doesn't matter
try
    m(isinf(m))=10*max(m(~isinf(m)));
catch
    % if this fails, it's because all the elements are infinite.  Replace
    % with the unsmoothed derivative.  There's probably a more elegant way
    % to do this.
    m = (dIdv + dIdu) ./ (dIdu-dIdv);
end


% Shorthand "b", which also happens to be the
% y intercept of the line of slope m that goes through each grid midpoint
b = ym - m.*xm;

% Weighting: weight by square of gradient magnitude and inverse
% distance to gradient intensity centroid.
sdI2 = sum(dImag2(:));
xcentroid = sum(sum(dImag2.*xm))/sdI2;
ycentroid = sum(sum(dImag2.*ym))/sdI2;
w  = dImag2./sqrt((xm-xcentroid).*(xm-xcentroid)+(ym-ycentroid).*(ym-ycentroid));

% least-squares minimization to determine the translated coordinate
% system origin (xc, yc) such that lines y = mx+b have
% the minimal total distance^2 to the origin:
% See function lsradialcenterfit (below)
[xc yc] = lsradialcenterfit(m, b, w);



%%
% Return output relative to upper left coordinate
xc = xc + (Nx+1)/2.0;
yc = yc + (Ny+1)/2.0;

% A rough measure of the particle width.
% Not at all connected to center determination, but may be useful for tracking applications;
% could eliminate for (very slightly) greater speed
Isub = I - min(I(:));
[px,py] = meshgrid(1:Nx,1:Ny);
xoffset = px - xc;
yoffset = py - yc;
r2 = xoffset.*xoffset + yoffset.*yoffset;
sigma = sqrt(sum(sum(Isub.*r2))/sum(Isub(:)))/2;  % second moment is 2*Gaussian width

%%

end

function [xc yc] = lsradialcenterfit(m, b, w)
% least squares solution to determine the radial symmetry center

% inputs m, b, w are defined on a grid
% w are the weights for each point
wm2p1 = w./(m.*m+1);
sw  = sum(sum(wm2p1));
smmw = sum(sum(m.*m.*wm2p1));
smw  = sum(sum(m.*wm2p1));
smbw = sum(sum(m.*b.*wm2p1));
sbw  = sum(sum(b.*wm2p1));
det = smw*smw - smmw*sw;
xc = (smbw*sw - smw*sbw)/det;    % relative to image center
yc = (smbw*smw - smmw*sbw)/det; % relative to image center

end