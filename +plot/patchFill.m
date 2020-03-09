function patchFill(x,y,color,xoffset,yoffset,vert,varargin)

ip = inputParser;
ip.addParameter('FaceAlpha', .5);
ip.addParameter('LineWidth', 1);

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


lw = ip.Results.LineWidth;


if vert
    patch(xPatch + xoffset, yPatch + yoffset,  color, 'facecolor', color, 'edgecolor', 'none', 'facealpha', ip.Results.FaceAlpha);
    hold on;
    % outline thicker, on top
    if lw > 0
        plot(xHist + xoffset, yHist + yoffset, '-', 'color', color, 'linewidth', lw);
    end
else
    patch(yPatch + yoffset, xPatch + xoffset,  color, 'facecolor', color, 'edgecolor', 'none', 'facealpha', ip.Results.FaceAlpha); hold on
    if lw > 0
        plot(yHist + yoffset, xHist + xoffset, '-', 'color', color, 'linewidth', lw);
    end
end