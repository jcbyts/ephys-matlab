function pdsList = getPdsFileList(sessionInfo)
% GET PDS FILE LIST finds all PDS files in a session or sessions
% Inputs:
%   sessionInfo@struct - struct(array) of session structs from
%                        io.loadSession(oepath)
% Output:
%   pdsList@cell - array of strings, paths to PDS files
% Example call:
%   pdsList = io.getPdsFileList(sessionInfo)

% 2017.08.14    jly     wrote it

% --- check if more than one session combined
nSessions = numel(sessionInfo);

% --- check if it is an ops struct (io.loadOps)
isOps = isfield(sessionInfo, 'root');

if isOps
    pathField = 'root';
else
    pathField = 'path';
end

% --- check that appropriate version of MATLAB is used
[~,d] = version();
assert(datenum(d) > datenum(2016, 0, 0), 'getPdsList: requires matlab 2016a or greater. dir() output contains folder subfield')

% --- build list of PDS files
pdsList = {};

for i = 1:nSessions
   fl = dir(fullfile(sessionInfo(i).(pathField), '*.PDS'));
   
   pdsList = [pdsList; arrayfun(@(x) fullfile(x.folder, x.name), fl, 'UniformOutput', false)]; %#ok<AGROW>
   
end