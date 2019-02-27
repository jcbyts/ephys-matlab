function clr = getNextPlotColor(ax)
% clr = getNextPlotColor(ax)

if nargin < 1
    ax = gca;
end

colorOrder = get(ax, 'ColorOrder');
next = mod(length(get(ax, 'Children')), size(colorOrder, 1))+1;

clr = colorOrder(next,:);