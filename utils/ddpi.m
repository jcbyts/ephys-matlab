function data = ddpi(fname, overwrite, debug, eyelinkIllumination)
%set(0,'DefaultFigureWindowStyle','docked')
% close all;
% clear all;

if nargin < 4
    eyelinkIllumination = false; % default assumption is no eyelink illumination
end

if nargin < 3 || isempty(debug)
    debug = false;
end

if nargin < 2 || isempty(overwrite)
    overwrite = false;
end


global DEBUG;
global P1_ROI_SIZE;
global P4_ROI_SIZE;
global THRESHOLD;
global PUPILTHRESH;
global IMG_WIDTH;
global IMG_HEIGHT;
global EYELINK;

EYELINK = eyelinkIllumination;
DEBUG = debug;
P1_ROI_SIZE = 30;
P4_ROI_SIZE = 20;
THRESHOLD = 75;
PUPILTHRESH = 2;


if nargin < 1
    fname = 'test.avi';
end

fprintf('Opening [%s]\n', fname)

[~, fn, ~] = fileparts(fname);
outfname = [fn 'tracersm.mat'];
if exist(outfname, 'file') && ~overwrite
    tmp = load(outfname);
    data = tmp.data;
    return
end

v = VideoReader(fname);

IMG_WIDTH = v.Width;
IMG_HEIGHT = v.Height;

tic
fprintf('Running tracking analysis...')
[p1x, p1y, p1offx, p1offy, p4x, p4y, p4offx, p4offy, pupx, pupy, time] = track(v);

fprintf('[%02.2fs]\n', toc)

data.x1 = p1x + p1offx;
data.x4 = p4x + p4offx;
data.y1 = p1y + p1offy;
data.y4 = p4y + p4offy;
data.xp = pupx;
data.yp = pupy;
data.time = time;



save(outfname, 'data');
end

function [p1x, p1y, p1offx, p1offy, p4x, p4y, p4offx, p4offy, pupx, pupy, time] = track(video)
global P1_ROI_SIZE;
global P4_ROI_SIZE;
global THRESHOLD;
global PUPILTHRESH;
global DEBUG;
global IMG_HEIGHT;
global IMG_WIDTH;
global EYELINK;

p1x = [];
p1y = [];
p1offx = [];
p1offy = [];

p4x = [];
p4y = [];
p4offx = [];
p4offy = [];

pupx = [];
pupy = [];

time = [];

frame = 1;
while 1
    if DEBUG
        disp([video.CurrentTime video.Duration])
    end
    
    if ~hasFrame(video)
        break
    end
    
    % store video file time
    time(frame) = video.CurrentTime;
    
    % read frame -> store on GPU
    img = readFrame(video);
    img = gpuArray(img(:, :, 1));
    
    % --- detect pupil ROI
    
    % threshold and smooth
    fimg = imgaussfilt(double(img<PUPILTHRESH), 2);
    
    % find pupil ROI
    pimg = bwlabel(fimg > 0, 8);
    pimg(pimg==mode(pimg(:,1))) = 0; % remove any big edge artifacts
    pimg(pimg==mode(pimg(:,end))) = 0; % remove any big edge artifacts
    pid = mode(pimg(pimg~=0));
    pimg = pimg==pid;
    xy = contour(pimg, [1 1]); % get contour around pupil
    label = bwlabel(xy(1,:)>1);
    pcontour = (xy(:,label==1));
    
    x = sum(pimg)>0;
    y = sum(pimg,2)>0;
    ixPupil = double(y(:))*double(x(:))';   
    [ii,jj] = find(ixPupil);
    ii = unique(ii);
    jj = unique(jj);
    
    if isempty(ii) % no pupil
        p1x(frame) = -1;
        p4x(frame) = -1;
        p1y(frame) = -1;
        p4y(frame) = -1;
        p1offx(frame) = -1;
        p1offy(frame) = -1;
        p4offx(frame) = -1;
        p4offy(frame) = -1;
        pupx(frame) = -1;
        pupy(frame) = -1;
        continue
    end
    
    % analysis restricted to pupil ROI
    img = img(ii,jj); % only analyze blobs in pupil
    
    % contour of the pupil outline
    pcontour(1,:) = pcontour(1,:) - gather(jj(1));
    pcontour(2,:) = pcontour(2,:) - gather(ii(1));
    
    % new image size
    [IMG_HEIGHT, IMG_WIDTH] = size(img);
    
    
    % --- high threshold to find putative purkinje images
    fimg = imgaussfilt(double(img), 1);
    [limg, nLabels] = bwlabel(fimg > THRESHOLD, 8);
    
    if (nLabels >= 2) % need at least two images for this to work
        
        try % incase of badness put analyses in a try-catch statement
            
            % --- find all blobs in ROI
            [Bbox, Areas] = findBlob(limg, nLabels);
            
            if DEBUG % plot blobs
                figure(1); clf
                subplot(3,1,1:2)
                imagesc(img); hold on
                for ii = 1:size(Bbox,1)
                    plot(Bbox(ii,[1 1 3 3 1]), Bbox(ii,[2 4 4 2 2]), 'y')
                end
                title('Detected Blobs')
                drawnow
            end
        
            % Edge Cases 1: weird aspect ratios
            
            % flag blobs with really weird aspect ratios (easy to remove)
            boxWidth = (Bbox(:,3)-Bbox(:,1));
            boxHeight = (Bbox(:,4)-Bbox(:,2));
            aspect = boxWidth ./ boxHeight;
            weirdshape = abs(aspect - 1) >= .5 & Areas > 10;
        
            boxCtr = [Bbox(:,1) + boxWidth/2, Bbox(:,2) + boxHeight/2];
        
            if DEBUG % plot removed blobs
                for ii = find(weirdshape)'
                    plot(Bbox(ii,[1 1 3 3 1]), Bbox(ii,[2 4 4 2 2]), 'b')
                end
                drawnow
            end
            
            Bbox(weirdshape,:) = [];
            Areas(weirdshape,:) = [];
            boxCtr(weirdshape,:) = [];
        
            % --- fit ellipse to the pupil 
            fit = fit_ellipse(pcontour(1,:)', pcontour(2,:)');
            
            % get polygon for pupil
            phi = fit.phi;
            cosphi = cos(phi);
            sinphi = sin(phi);
            R = [cosphi sinphi; -sinphi cosphi];
            thetas = linspace(0, 2*pi, 100);
            ex = fit.X0 + fit.a * cos(thetas);
            ey = fit.Y0 + fit.b * sin(thetas);
            epoly = R * [ex; ey];
                
            % store pupil center
            pupx(frame) = fit.X0_in;
            pupy(frame) = fit.Y0_in;
        
            % only analyze blobs within the pupil ellipse
            in = find(inpolygon(boxCtr(:,1), boxCtr(:,2),epoly(1,:), epoly(2,:)));
        if DEBUG
            plot(epoly(1,:), epoly(2,:), 'r')
            plot(boxCtr(:,1), boxCtr(:,2), 'c.')
            for ii = in(:)'
                plot(Bbox(ii,[1 1 3 3 1]), Bbox(ii,[2 4 4 2 2]), 'm')
            end
        end
        
        if numel(in) < 2
            error
        end
        
        % P1 is the biggest blob in the pupil (findBlob is sorted by area)
        p1ix = in(1);
        
        % p4 is the smallest blob in the pupil
        p4ix = in(end);
        
        % Edge cases
        if EYELINK
            if numel(in) > 2 % more than two blobs in pupil
                bigblobs = Areas(in) > 100;
                bigblobs = in(bigblobs);
                if numel(bigblobs) > 1 % more than one big blob
                    [~, ix] = max(Bbox(bigblobs,1)); % find the bib blob that is furthest to the right
                    p1ix = bigblobs(ix);
                end
            end
        end
        
        [sx, ex, sy, ey] = ROI(Bbox(p1ix, :), P1_ROI_SIZE);
        
        p1Img = fimg(sy:ey, sx:ex);
        [p1x(frame), p1y(frame)] = radialcenter(p1Img);
%         [p1x(frame), p1y(frame)] = centerOfMass(p1Img, 200);
        p1offx(frame) = sx;
        p1offy(frame) = sy;
        
        
        if sum(Areas(in) < 60) > 1
            [~, p4ix] = max(Bbox(in,1),[],1);
            p4ix = in(p4ix);
        end
     
        [sx, ex, sy, ey] = ROI(Bbox(p4ix, :), P4_ROI_SIZE);
        p4Img = fimg(sy:ey, sx:ex);
        [p4x(frame), p4y(frame)] = radialcenter(p4Img);
%         [p4x(frame), p4y(frame)] = centerOfMassWithResample(p4Img, 125, 10);
        p4offx(frame) = sx;
        p4offy(frame) = sy;
         
        
        if DEBUG
           
            subplot(3,1,1:2)
            plot(epoly(1,:), epoly(2,:), 'r')
            plot(Bbox(p1ix,[1 1 3 3 1]), Bbox(p1ix,[2 4 4 2 2]), 'r')
            plot(Bbox(p4ix,[1 1 3 3 1]), Bbox(p4ix,[2 4 4 2 2]), 'g')
            
            plot(p1x(frame)+p1offx(frame)-1, p1y(frame)+p1offy(frame)-1, '+r', 'MarkerSize', 10);
            plot(p4x(frame)+p4offx(frame)-1, p4y(frame)+p4offy(frame)-1, '+g', 'MarkerSize', 10);
            plot(pupx(frame),pupy(frame), 'm+')
            title(frame)
            subplot(3,2,5)
            imagesc(p1Img); hold on
            plot(p1x(frame), p1y(frame), '+r', 'MarkerSize', 10)
            subplot(3,2,6)
            imagesc(p4Img); hold on
            plot(p4x(frame), p4y(frame), '+g', 'MarkerSize', 10)
            drawnow
        end
        
        catch
        
        p1x(frame) = -1;
        p4x(frame) = -1;
        p1y(frame) = -1;
        p4y(frame) = -1;
        p1offx(frame) = -1;
        p1offy(frame) = -1;
        p4offx(frame) = -1;
        p4offy(frame) = -1;
        pupx(frame) = -1;
        pupy(frame) = -1;
        end
        
    else
        p1x(frame) = -1;
        p4x(frame) = -1;
        p1y(frame) = -1;
        p4y(frame) = -1;
        p1offx(frame) = -1;
        p1offy(frame) = -1;
        p4offx(frame) = -1;
        p4offy(frame) = -1;
        pupx(frame) = -1;
        pupy(frame) = -1;
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

function [Bbox, Areas] = findBlob(limg, nLabels)
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
[Areas, I] = sort(Areas, 'descend');
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
