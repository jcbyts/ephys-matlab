function info = loadEphysInfo(fpath)
% LOADEPHYSINFO loads the meta data required for reconstructing the timestamps/voltage of raw ephys data
% Inputs:
%   fpath@string    - path to session
% Outputs:
%   info@struct      - info struct. built by io.oe2dat --> meta-data required to make sense of raw ephys samples
% Example Call:
%   oepath = 'C:\Data\Ellie_2017-08-09_13-04-23_ShankD15MT6';
%   info   = io.loadEphysInfo(oepath);

% 2017.08.14    jly     wrote it

warning('off') % turn off warnings (some variable names might cause flags)
if isstruct(fpath) % it's an ops
    fpath = fpath.root;
end

shankDirs = dir(fullfile(fpath, '_shank*'));

if isempty(shankDirs)
    fname = fullfile(fpath, 'ephys_info.mat');
    if exist(fname, 'file')
        info = load(fname);
    else
        info = [];
    end
    return
end

for iShank = 1:numel(shankDirs)
    tmp = load(fullfile(fpath, shankDirs(iShank).name, 'ephys_info.mat'));

    info(iShank) = tmp; %#ok<AGROW>
end

warning('on') % turn warnings back on