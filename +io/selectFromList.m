function val = selectFromList(list)
% list is a cell array of strings

if nargin == 0
    warning('need to pass in a cell array of strings to select')
    val = [];
    return
end

l = io.listHolder;

fig = uifigure('Name', 'Electrode Select GUI', 'Position', [100 500 350 400]); % 'CloseRequestFcn', @(fig, event) myclosereq(fig,lbx)
fig.Color = [1 1 1];
lbl = uitextarea(fig, 'Position', [10 330 300 60], 'Value', {'Select an electrode. Close figure to confirm. If you need to create a new configuration, edit harwarde.electrodeFactory'});
lbl.BackgroundColor = [1 1 1];
lbl.HorizontalAlignment = 'center';
% text(fig, 
lbx = uilistbox(fig, 'Position', [10 10 300 330], 'Items', list, 'MultiSelect', 'on', 'ValueChangedFcn', @(lbx, event) selectionChanged(lbx,l)); %

uiwait(fig)
val = l.list;
% 
function selectionChanged(lbx,l)
l.list = lbx.Value;


