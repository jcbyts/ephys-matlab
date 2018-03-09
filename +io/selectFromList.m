function val = selectFromList(list)
% list is a cell array of strings

if nargin == 0
    warning('need to pass in a cell array of strings to select')
    val = [];
    return
end

l = io.listHolder;

fig = uifigure('Name', 'Click to select. Close figure to continue', 'Position', [100 500 350 400]); % 'CloseRequestFcn', @(fig, event) myclosereq(fig,lbx)
lbx = uilistbox(fig, 'Position', [10 10 300 350], 'Items', list, 'MultiSelect', 'on', 'ValueChangedFcn', @(lbx, event) selectionChanged(lbx,l)); %

uiwait(fig)
val = l.list;
% 
function selectionChanged(lbx,l)
l.list = lbx.Value;


