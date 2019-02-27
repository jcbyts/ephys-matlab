function patchFill3D(x,y,color,xoffset,yoffset,vert,varargin)
% patchFill3D(x,y,color,xoffset,yoffset,vert,varargin)
ip = inputParser;
ip.addParameter('FaceAlpha', .5);
ip.parse(varargin{:})

if nargin < 3 || isempty(color)
    color = 'k';
end

if nargin < 4 || isempty(xoffset)
    xoffset = 0;
end

if nargin < 5 || isempty(yoffset)
    yoffset = 0;
end

if nargin < 6 || isempty(vert)
    vert = true;
end

% plotArgs = varargin;

xHist = y(:);
yHist = x(:);

[yHist, ind] = sort(yHist);
xHist = xHist(ind);
yPatch = [yHist(1); yHist(:); yHist(end)];
xPatch = [0; xHist(:); 0];
if ip.Results.FaceAlpha==0
    lw = 2;
else
    lw = .5;
end

if vert
    fill3(zeros(size(xPatch)), xPatch + xoffset, yPatch + yoffset,  color, 'facecolor', color, 'edgecolor', 'none', 'facealpha', ip.Results.FaceAlpha);
    hold on;
    % outline thicker, on top
    plot3(zeros(size(xHist)), xHist + xoffset, yHist + yoffset, '-', 'color', color, 'linewidth', lw);
else
    fill3(zeros(size(xPatch)) + xoffset, yPatch + yoffset, xPatch,  color, 'facecolor', color, 'edgecolor', 'none', 'facealpha', ip.Results.FaceAlpha); hold on
    plot3(zeros(size(xHist)) + xoffset, yHist + yoffset, xHist, '-', 'color', color, 'linewidth', lw);
end