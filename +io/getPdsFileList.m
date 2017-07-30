function pdsList = getPdsFileList(sessionInfo)

nSessions = numel(sessionInfo);

isOps = isfield(sessionInfo, 'root');

if isOps
    pathField = 'root';
else
    pathField = 'path';
end

pdsList = {};

for i = 1:nSessions
   fl = dir(fullfile(sessionInfo(i).(pathField), '*.PDS'));
   
   pdsList = [pdsList; arrayfun(@(x) fullfile(x.folder, x.name), fl, 'UniformOutput', false)]; %#ok<AGROW>
   
end