function ops = defaultOps(ops)

if nargin < 1 || isempty(ops)
    ops = struct();
end

[~,thisComputer]=system('hostname');
switch thisComputer(1:end-1)
    case 'Gravedigger'
        ops.GPU         = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)
        ops.parfor      = 1; % whether to use parfor to accelerate some parts of the algorithm
        ops.verbose     = 1; % whether to print command line progress
        ops.showfigures = 0; % whether to plot figures during optimization
    otherwise
        ops.GPU         = 0; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)
        ops.parfor      = 0; % whether to use parfor to accelerate some parts of the algorithm
        ops.verbose     = 1; % whether to print command line progress
        ops.showfigures = 1; % whether to plot figures during optimization
end

ops.NchanTOT = ops.Nchan;

info = load(fullfile(ops.root, 'ephys_info.mat'));