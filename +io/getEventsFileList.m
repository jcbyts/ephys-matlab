function evList = getEventsFileList(sessionInfo)
% evList = getEventsFileList(sessionInfo)

nSessions = numel(sessionInfo);

isOps = isfield(sessionInfo, 'root');

if isOps
    pathField = 'root';
else
    pathField = 'path';
end

evList = {};

for i = 1:nSessions
   fl = dir(fullfile(sessionInfo(i).(pathField), 'all_channels*'));
   
   evList = [evList; arrayfun(@(x) fullfile(x.folder, x.name), fl, 'UniformOutput', false)]; %#ok<AGROW>
   
end