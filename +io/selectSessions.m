function [dataDir, val] = selectSessions(dataDir)

if nargin<1
    dataDir = uigetdir(pwd, 'Select Data Directory');
end

fl = dir(dataDir);
fl(1:2) = []; % remove relative directories
fl(~arrayfun(@(x) x.isdir, fl)) = []; % remove non-directory items

l = io.listHolder;

fig = uifigure('Name', 'Select Sessions. Close figure to continue', 'Position', [100 500 350 400]); % 'CloseRequestFcn', @(fig, event) myclosereq(fig,lbx)
lbx = uilistbox(fig, 'Position', [10 10 300 350], 'Items', {fl.name}, 'MultiSelect', 'on', 'ValueChangedFcn', @(lbx, event) selectionChanged(lbx,l)); %

uiwait(fig)
val = l.list;
% 
function selectionChanged(lbx,l)
l.list = lbx.Value;


