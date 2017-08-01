function runKiloSort(ops, varargin)

ip = inputParser();
ip.addOptional('merge', false)
ip.addOptional('GPU', true)
ip.addOptional('parfor', true)
ip.addOptional('verbose', true)
ip.addOptional('showfigures', true)

ip.parse(varargin{:});

if nargin < 2
    merge = false;
end

fprintf('-----------------------------------------------------------------\n')
fprintf('-----------------------------------------------------------------\n')
fprintf('-----------------------------------------------------------------\n')
fprintf('Running Kilosort\n')

if isfield(ops, 'parfor_')
    ops.parfor = ops.parfor_;
    ops = rmfield(ops, 'parfor_');
end

merge = ip.Results.merge;
ops.GPU         = ip.Results.GPU;
ops.parfor      = ip.Results.parfor;
ops.verbose     = ip.Results.verbose;
ops.showfigures = ip.Results.showfigures;

if ops.GPU
    fprintf('Using GPU for faster processing\n')
    fprintf('Opening gpuDevice. If this is the first time this was run during \nthis matlab session, this can be very slow\n')
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

% [path_to, file, ext]=fileparts(ops.fbinary);
% fOut=fullfile(path_to, [file '_hp' ext]);
% 
% % if ~exist(fOut, 'file')
% %     ops = removeArtifacts(ops);
% % end
% 
% ops.fbinary = fOut;

% main spike-sorting routine
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% saving
fprintf('saving matlab results file\n')
save(fullfile(ops.root,  'rez.mat'), 'rez', 'ops', '-v7.3');

if merge
    fprintf('Attemptying Automerge\n')
    rez                = merge_posthoc2(rez);
end

fprintf('saving python files for Phy\n')

% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);