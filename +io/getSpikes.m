function sp = getSpikes(sess, tag)
% GET SPIKES gets spike times and cluster info
% Inputs:
%   session@struct - session struct
% Output:
%   sp@cell     - cell array of spike structs
% Example call:
%   sp = io.getSpikes(sess)

% 2017.08.14    jly     wrote it
if nargin < 2
    tag = 'Kilo';
end

ephys_dirs = dir(fullfile(sess.path, '_shank*'));
nDirs = numel(ephys_dirs);

sp = {};
for i = 1:nDirs
    fspike = fullfile(sess.path, ephys_dirs(i).name, sprintf('sp-%s.mat', tag));
    if exist(fspike, 'file')
        sp{i} = load(fspike);
    end
end