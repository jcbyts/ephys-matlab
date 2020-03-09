function data = ddpi_offline(fname, varargin)
% dDPI analysis of video file
% Inputs:
%   fname <string> - full path to movie file (*.avi, *.mp4)
%   
% Optional Arguments (as argument pairs):
%   'overwrite'     <logical> - overwrite existing analysis? (default: false)
%   'debug'         <logical> - debugger output (default: false)
%   'p1ROI'         <scalar>  - size of P1 ROI (default: 30 pixels)
%   'p4ROI'         <scalar>  - size of P4 ROI (default: 20 pixels)
%   'fitPupil'      <logical> - fit pupil with ellipse and restrict
%                               analysis to within pupil (default: false)
%   'pupilThresh'   <scalar>  - threshold to detect pupil (default: 2)
%   'piThresh'      <scalar>  - threshold to detect purkinje images (default: 100)
%   'eyelinkOn'     <logical> - was the eyelink illuminator on? (default: false)
%   'savevideo'     <logical> - save a video of the analysis? (default: false)
%   'method'        <string>  - method for center tracking ('rsm' or 'com')
%                               rsm - radial symmetric mean (default)
%                               com - center of mass (less robust to noise)
%   'manualROID'    <array>   - manually specified ROI [x1, y1, x2, y2] 
%
%  Outputs:
%   data <struct>
%       time    [1 x nFrames] - time in the video file at each frame
%       x1      [1 x nFrames] - x position of 1st purkinje image
%       x4      [1 x nFrames] - x position of 4th purkinje image
%       y1      [1 x nFrames] - y position of 1st purkinje image
%       y4      [1 x nFrames] - y position of 4st purkinje image
%       xpup    [1 x nFrames] - x position of pupil
%       ypup    [1 x nFrames] - y position of pupil
%
% Example Call:
%   track = ddpi_offline('test.avi');
%
%   track = ddpi_offline('test.avi', 'overwrite', true, 'debug', true, ...
%       'piThresh', 128);

ip = inputParser(); % parse optional arguments
ip.addParameter('overwrite', false)
ip.addParameter('debug', false)
ip.addParameter('p1ROI', 30)
ip.addParameter('p4ROI', 20)
ip.addParameter('pupilThresh', 2)
ip.addParameter('piThresh', 100)
ip.addParameter('eyelinkOn', false)
ip.addParameter('saveVideo', false)
ip.addParameter('fitPupil', false)
ip.addParameter('method', 'rsm', @(x) strcmp(x, {'rsm', 'com'}))
ip.addParameter('manualROI', [])
ip.parse(varargin{:})

global DEBUG;
global P1_ROI_SIZE;
global P4_ROI_SIZE;
global THRESHOLD;
global PUPILTHRESH;
global IMG_WIDTH;
global IMG_HEIGHT;
global EYELINK;
global METHOD;
global SAVEVIDEO;
global FITPUPIL;
global MANUALROI;

EYELINK     = ip.Results.eyelinkOn;
DEBUG       = ip.Results.debug;
P1_ROI_SIZE = ip.Results.p1ROI;
P4_ROI_SIZE = ip.Results.p4ROI;
THRESHOLD   = ip.Results.piThresh;
PUPILTHRESH = ip.Results.pupilThresh;
METHOD      = ip.Results.method;
SAVEVIDEO   = ip.Results.saveVideo;
FITPUPIL    = ip.Results.fitPupil;
MANUALROI   = ip.Results.manualROI;

if nargin < 1
    fname = uigetfile(fullfile(pwd,'*.avi;*.mp4;*.mov'), 'Select video file');
end

% track options
optnames = fieldnames(ip.Results);
optvals = cellfun(@(x) num2str(ip.Results.(x), '%1.0f'), optnames, 'uni', 0);
optsstr = [cellfun(@(x) x(1:2), optnames, 'uni', 0) optvals]';
optsstr = sprintf('_%s%s', optsstr{:});

[~, fn, ~] = fileparts(fname);
outfname = [fn optsstr 'trace.mat']; % analysis file name
videoout = [fn optsstr 'trace.mp4']; % video file name
    
if exist(outfname, 'file') && ~overwrite && ( (SAVEVIDEO && exist(videoout, 'file')) || ~SAVEVIDEO)
    tmp = load(outfname);
    data = tmp.data;
    return
end

vidout = [];

if SAVEVIDEO
    vidout = VideoWriter(videoout);
    vidout.Quality = 100; % high quality output
    open(vidout)
end

fprintf('Opening [%s]\n', fn)

v = VideoReader(fname);

IMG_WIDTH  = v.Width;
IMG_HEIGHT = v.Height;

if isempty(MANUALROI)
    MANUALROI = [0 0 IMG_WIDTH IMG_HEIGHT];
end

tic
fprintf('Running tracking analysis...')
[p1x, p1y, p1offx, p1offy, p4x, p4y, p4offx, p4offy, pupx, pupy, time] = track(v, vidout);
fprintf('[%02.2fs]\n', toc)

if SAVEVIDEO
    close(vidout)
end

% format output
data.time = time;
data.x1 = p1x + p1offx;
data.x4 = p4x + p4offx;
data.y1 = p1y + p1offy;
data.y4 = p4y + p4offy;
data.xp = pupx;
data.yp = pupy;

save(outfname, 'data');
end

function [p1x, p1y, p1offx, p1offy, p4x, p4y, p4offx, p4offy, pupx, pupy, time] = track(video, vidout)
% main tracking analysis
% Steps:
% 1) Read video frame
% 2) detect pupil ROI / fit pupil (optional)
% 3) blob analysis to detect putative purkinje images
% 4) track P1 and P4 with radial symmetric mean or center of mass

global P1_ROI_SIZE;
global P4_ROI_SIZE;
global THRESHOLD;
global PUPILTHRESH;
global DEBUG;
global IMG_HEIGHT;
global IMG_WIDTH;
global EYELINK;
global METHOD;
global SAVEVIDEO;
global FITPUPIL;
global MANUALROI;

if nargin < 2 || isempty(vidout)
    vidout = [];
    assert(~SAVEVIDEO, 'save video is on, but there is no video writer object')
end

    
% this is not ideal - should preallocate if possible
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

frame = 1; % frame iterator

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
    
    if FITPUPIL
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
    end % if fitpupil
    
    % --- high threshold to find putative purkinje images
    fimg = imgaussfilt(double(img), 1);
    [limg, nLabels] = bwlabel(fimg > THRESHOLD, 8);
    
    if (nLabels >= 2) % need at least two images for this to work
        
        try % incase of badness put analyses in a try-catch statement
            
            % --- find all blobs in ROI
            [Bbox, Areas] = findBlob(limg, nLabels);
            
            if DEBUG || SAVEVIDEO % plot blobs
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
            
            Bbox(weirdshape,:)   = [];
            Areas(weirdshape,:)  = [];
            boxCtr(weirdshape,:) = [];
        
            % --- fit ellipse to the pupil 
            if FITPUPIL
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
            else
                epoly = [MANUALROI([1 1 3 3 1]); MANUALROI([2 4 4 2 2])];
                pupx(frame) = -1;
                pupy(frame) = -1;
            end
        
            % only analyze blobs within the pupil ellipse (or manual ROI)
            in = find(inpolygon(boxCtr(:,1), boxCtr(:,2),epoly(1,:), epoly(2,:)));
            
            if DEBUG || SAVEVIDEO % plot pupil ellipse and blob centers
                plot(epoly(1,:), epoly(2,:), 'r')
                plot(boxCtr(:,1), boxCtr(:,2), 'c.')
                for ii = in(:)'
                    plot(Bbox(ii,[1 1 3 3 1]), Bbox(ii,[2 4 4 2 2]), 'm')
                end
                drawnow
            end
        
            if numel(in) < 2 % less than 2 putative purkinje images
                error % break the try statement
            end
        
            % P1 is the biggest blob in the pupil (findBlob is sorted by area)
            p1ix = in(1);
        
            % p4 is the smallest blob in the pupil
            p4ix = in(end);
        
            % Edge cases 2: the eyelink illuminator
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
            switch METHOD
                case 'rsm'
                    [p1x(frame), p1y(frame)] = radialcenter(p1Img);
                case 'com'
                    [p1x(frame), p1y(frame)] = centerOfMass(p1Img, 200);
            end
            
            p1offx(frame) = sx;
            p1offy(frame) = sy;
        
            % Edge case 3: lots of small blobs
            % take the one that is furthest to the right within the ROI
            if sum(Areas(in) < 60) > 1
                [~, ix] = max(Bbox(in,1),[],1);
                p4ix = in(ix);
            end
     
            [sx, ex, sy, ey] = ROI(Bbox(p4ix, :), P4_ROI_SIZE);
            p4Img = fimg(sy:ey, sx:ex);
            
            switch METHOD
                case 'rsm'
                    [p4x(frame), p4y(frame)] = radialcenter(p4Img);
                case 'com'
                    [p4x(frame), p4y(frame)] = centerOfMassWithResample(p4Img, 125, 10);
            end
            
            p4offx(frame) = sx;
            p4offy(frame) = sy;
         
        
            if DEBUG || SAVEVIDEO
                % plot the output of the analysis
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
            
            if SAVEVIDEO
                currFrame = getframe(gcf);
                writeVideo(vidout, currFrame);
            end
            
        catch
            
            p1x(frame) = -1; %#ok<*AGROW>
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



function [xc, yc, sigma] = radialcenter(I)
I = gather(double(I));
% Number of grid points
[Ny, Nx] = size(I);

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
[xc, yc] = lsradialcenterfit(m, b, w);



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

function [xc, yc] = lsradialcenterfit(m, b, w)
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


function ellipse_t = fit_ellipse( x,y,axis_handle )
%
% fit_ellipse - finds the best fit to an ellipse for the given set of points.
%
% Format:   ellipse_t = fit_ellipse( x,y,axis_handle )
%
% Input:    x,y         - a set of points in 2 column vectors. AT LEAST 5 points are needed !
%           axis_handle - optional. a handle to an axis, at which the estimated ellipse 
%                         will be drawn along with it's axes
%
% Output:   ellipse_t - structure that defines the best fit to an ellipse
%                       a           - sub axis (radius) of the X axis of the non-tilt ellipse
%                       b           - sub axis (radius) of the Y axis of the non-tilt ellipse
%                       phi         - orientation in radians of the ellipse (tilt)
%                       X0          - center at the X axis of the non-tilt ellipse
%                       Y0          - center at the Y axis of the non-tilt ellipse
%                       X0_in       - center at the X axis of the tilted ellipse
%                       Y0_in       - center at the Y axis of the tilted ellipse
%                       long_axis   - size of the long axis of the ellipse
%                       short_axis  - size of the short axis of the ellipse
%                       status      - status of detection of an ellipse
%
% Note:     if an ellipse was not detected (but a parabola or hyperbola), then
%           an empty structure is returned

% =====================================================================================
%                  Ellipse Fit using Least Squares criterion
% =====================================================================================
% We will try to fit the best ellipse to the given measurements. the mathematical
% representation of use will be the CONIC Equation of the Ellipse which is:
% 
%    Ellipse = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0
%   
% The fit-estimation method of use is the Least Squares method (without any weights)
% The estimator is extracted from the following equations:
%
%    g(x,y;A) := a*x^2 + b*x*y + c*y^2 + d*x + e*y = f
%
%    where:
%       A   - is the vector of parameters to be estimated (a,b,c,d,e)
%       x,y - is a single measurement
%
% We will define the cost function to be:
%
%   Cost(A) := (g_c(x_c,y_c;A)-f_c)'*(g_c(x_c,y_c;A)-f_c)
%            = (X*A+f_c)'*(X*A+f_c) 
%            = A'*X'*X*A + 2*f_c'*X*A + N*f^2
%
%   where:
%       g_c(x_c,y_c;A) - vector function of ALL the measurements
%                        each element of g_c() is g(x,y;A)
%       X              - a matrix of the form: [x_c.^2, x_c.*y_c, y_c.^2, x_c, y_c ]
%       f_c            - is actually defined as ones(length(f),1)*f
%
% Derivation of the Cost function with respect to the vector of parameters "A" yields:
%
%   A'*X'*X = -f_c'*X = -f*ones(1,length(f_c))*X = -f*sum(X)
%
% Which yields the estimator:
%
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%       |  A_least_squares = -f*sum(X)/(X'*X) ->(normalize by -f) = sum(X)/(X'*X)  |
%       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% (We will normalize the variables by (-f) since "f" is unknown and can be accounted for later on)
%  
% NOW, all that is left to do is to extract the parameters from the Conic Equation.
% We will deal the vector A into the variables: (A,B,C,D,E) and assume F = -1;
%
%    Recall the conic representation of an ellipse:
% 
%       A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
% 
% We will check if the ellipse has a tilt (=orientation). The orientation is present
% if the coefficient of the term "x*y" is not zero. If so, we first need to remove the
% tilt of the ellipse.
%
% If the parameter "B" is not equal to zero, then we have an orientation (tilt) to the ellipse.
% we will remove the tilt of the ellipse so as to remain with a conic representation of an 
% ellipse without a tilt, for which the math is more simple:
%
% Non tilt conic rep.:  A`*x^2 + C`*y^2 + D`*x + E`*y + F` = 0
%
% We will remove the orientation using the following substitution:
%   
%   Replace x with cx+sy and y with -sx+cy such that the conic representation is:
%   
%   A(cx+sy)^2 + B(cx+sy)(-sx+cy) + C(-sx+cy)^2 + D(cx+sy) + E(-sx+cy) + F = 0
%
%   where:      c = cos(phi)    ,   s = sin(phi)
%
%   and simplify...
%
%       x^2(A*c^2 - Bcs + Cs^2) + xy(2A*cs +(c^2-s^2)B -2Ccs) + ...
%           y^2(As^2 + Bcs + Cc^2) + x(Dc-Es) + y(Ds+Ec) + F = 0
%
%   The orientation is easily found by the condition of (B_new=0) which results in:
% 
%   2A*cs +(c^2-s^2)B -2Ccs = 0  ==> phi = 1/2 * atan( b/(c-a) )
%   
%   Now the constants   c=cos(phi)  and  s=sin(phi)  can be found, and from them
%   all the other constants A`,C`,D`,E` can be found.
%
%   A` = A*c^2 - B*c*s + C*s^2                  D` = D*c-E*s
%   B` = 2*A*c*s +(c^2-s^2)*B -2*C*c*s = 0      E` = D*s+E*c 
%   C` = A*s^2 + B*c*s + C*c^2
%
% Next, we want the representation of the non-tilted ellipse to be as:
%
%       Ellipse = ( (X-X0)/a )^2 + ( (Y-Y0)/b )^2 = 1
%
%       where:  (X0,Y0) is the center of the ellipse
%               a,b     are the ellipse "radiuses" (or sub-axis)
%
% Using a square completion method we will define:
%       
%       F`` = -F` + (D`^2)/(4*A`) + (E`^2)/(4*C`)
%
%       Such that:    a`*(X-X0)^2 = A`(X^2 + X*D`/A` + (D`/(2*A`))^2 )
%                     c`*(Y-Y0)^2 = C`(Y^2 + Y*E`/C` + (E`/(2*C`))^2 )
%
%       which yields the transformations:
%       
%           X0  =   -D`/(2*A`)
%           Y0  =   -E`/(2*C`)
%           a   =   sqrt( abs( F``/A` ) )
%           b   =   sqrt( abs( F``/C` ) )
%
% And finally we can define the remaining parameters:
%
%   long_axis   = 2 * max( a,b )
%   short_axis  = 2 * min( a,b )
%   Orientation = phi
%
%

% initialize
orientation_tolerance = 1e-3;

% empty warning stack
warning( '' );

% prepare vectors, must be column vectors
x = x(:);
y = y(:);

% remove bias of the ellipse - to make matrix inversion more accurate. (will be added later on).
mean_x = mean(x);
mean_y = mean(y);
x = x-mean_x;
y = y-mean_y;

% the estimation for the conic equation of the ellipse
X = [x.^2, x.*y, y.^2, x, y ];
a = sum(X)/(X'*X);

% check for warnings
if ~isempty( lastwarn )
    disp( 'stopped because of a warning regarding matrix inversion' );
    ellipse_t = [];
    return
end

% extract parameters from the conic equation
[a,b,c,d,e] = deal( a(1),a(2),a(3),a(4),a(5) );

% remove the orientation from the ellipse
if ( min(abs(b/a),abs(b/c)) > orientation_tolerance )
    
    orientation_rad = 1/2 * atan( b/(c-a) );
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    [a,b,c,d,e] = deal(...
        a*cos_phi^2 - b*cos_phi*sin_phi + c*sin_phi^2,...
        0,...
        a*sin_phi^2 + b*cos_phi*sin_phi + c*cos_phi^2,...
        d*cos_phi - e*sin_phi,...
        d*sin_phi + e*cos_phi );
    [mean_x,mean_y] = deal( ...
        cos_phi*mean_x - sin_phi*mean_y,...
        sin_phi*mean_x + cos_phi*mean_y );
else
    orientation_rad = 0;
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
end

% check if conic equation represents an ellipse
test = a*c;
switch (1)
case (test>0),  status = '';
case (test==0), status = 'Parabola found';  warning( 'fit_ellipse: Did not locate an ellipse' );
case (test<0),  status = 'Hyperbola found'; warning( 'fit_ellipse: Did not locate an ellipse' );
end

% if we found an ellipse return it's data
if (test>0)
    
    % make sure coefficients are positive as required
    if (a<0), [a,c,d,e] = deal( -a,-c,-d,-e ); end
    
    % final ellipse parameters
    X0          = mean_x - d/2/a;
    Y0          = mean_y - e/2/c;
    F           = 1 + (d^2)/(4*a) + (e^2)/(4*c);
    [a,b]       = deal( sqrt( F/a ),sqrt( F/c ) );    
    long_axis   = 2*max(a,b);
    short_axis  = 2*min(a,b);

    % rotate the axes backwards to find the center point of the original TILTED ellipse
    R           = [ cos_phi sin_phi; -sin_phi cos_phi ];
    P_in        = R * [X0;Y0];
    X0_in       = P_in(1);
    Y0_in       = P_in(2);
    
    % pack ellipse into a structure
    ellipse_t = struct( ...
        'a',a,...
        'b',b,...
        'phi',orientation_rad,...
        'X0',X0,...
        'Y0',Y0,...
        'X0_in',X0_in,...
        'Y0_in',Y0_in,...
        'long_axis',long_axis,...
        'short_axis',short_axis,...
        'status','' );
else
    % report an empty structure
    ellipse_t = struct( ...
        'a',[],...
        'b',[],...
        'phi',[],...
        'X0',[],...
        'Y0',[],...
        'X0_in',[],...
        'Y0_in',[],...
        'long_axis',[],...
        'short_axis',[],...
        'status',status );
end

% check if we need to plot an ellipse with it's axes.
if (nargin>2) && ~isempty( axis_handle ) && (test>0)
    
    % rotation matrix to rotate the axes with respect to an angle phi
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    
    % the axes
    ver_line        = [ [X0 X0]; Y0+b*[-1 1] ];
    horz_line       = [ X0+a*[-1 1]; [Y0 Y0] ];
    new_ver_line    = R*ver_line;
    new_horz_line   = R*horz_line;
    
    % the ellipse
    theta_r         = linspace(0,2*pi);
    ellipse_x_r     = X0 + a*cos( theta_r );
    ellipse_y_r     = Y0 + b*sin( theta_r );
    rotated_ellipse = R * [ellipse_x_r;ellipse_y_r];
    
    % draw
    hold_state = get( axis_handle,'NextPlot' );
    set( axis_handle,'NextPlot','add' );
    plot( new_ver_line(1,:),new_ver_line(2,:),'r' );
    plot( new_horz_line(1,:),new_horz_line(2,:),'r' );
    plot( rotated_ellipse(1,:),rotated_ellipse(2,:),'r' );
    set( axis_handle,'NextPlot',hold_state );
end

end
