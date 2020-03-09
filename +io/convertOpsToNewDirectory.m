function ops = convertOpsToNewDirectory(ops_old, DIRECTORY)
% convertOpsToNewDirectory takes an ops struct and renames all paths to
% point to a new directory. It's useful when copying data from one location
% to another.
% ops_new = convertOpsToNewDirectory(ops_old, DIRECTORY)
%   where DIRECTORY is the path to a new location where the data have been
%   copied

pathfields = {'chanMap', 'fbinary', 'fproc', 'root'};

ops = ops_old;

% where is the end of the new directory
newPath = regexp(DIRECTORY, '\\', 'split');

for iField = 1:numel(pathfields)
    
    opspath = ops.(pathfields{iField});
    
    % match the end of the new directory to existing ops paths and fix the
    % head
    opsParts = regexp(opspath, '\\', 'split');
    pathend = find(strcmp(newPath{end}, opsParts));
    if isempty(pathend)
        disp('no path match')
        pathend = find(cellfun(@(x) strcmp(x(1), '_'), opsParts)) - 1;
    end
    
    pathargs = [newPath, opsParts(pathend+1:end)];
    newopspath = fullfile(pathargs{:});
    
    ops.(pathfields{iField}) = newopspath;
    
end
    