function ops = loadOps(fpath)
% LOADOPS loads the ops struct that is required for KiloSort (and holds directory information about the raw data)
% Inputs:
%   fpath@string    - path to session
% Outputs:
%   ops@struct      - ops struct. built by io.oe2dat --> used for KiloSort
% Example Call:
%   oepath = 'C:\Data\Ellie_2017-08-09_13-04-23_ShankD15MT6';
%   ops = io.loadOps(oepath);

if isa(fpath, 'table')
    fpath = fullfile(getpref('EPHYS', 'SERVER_DATA'), fpath.Directory{1});
end

warning('off') % suppress warning while loading
shankDirs = dir(fullfile(fpath, '_shank*'));

if isempty(shankDirs)
    ops = [];
    return
end

for iShank = 1:numel(shankDirs)
    
    tmp = loadops_helper(fullfile(fpath, shankDirs(iShank).name));
    ops(iShank) = tmp; %#ok<AGROW>
end
warning('on') % warnings back on


function tmp = loadops_helper(fpath)
    tmp = load(fullfile(fpath, 'ops.mat'));
    % parfor is a special function so matlab doesn't let you load it as a
    % variable. Correct for this here. In theory we should've called it
    % something else, but this is what Kilosort uses.
    if isfield(tmp, 'parfor_')
        tmp.parfor = tmp.parfor_;
        tmp = rmfield(tmp, 'parfor_');
    else
        tmp.parfor = false;
    end