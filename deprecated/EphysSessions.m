function EphysSessions
% EPHYS SESSIONS concatenates multiples sessions together


dataPath = io.concatenateSessions();

%%
dataPath = 'C:\Data\Ellie_3-sessions_2017-07-24_2017-07-26\';


%% KiloSort

ops = load(fullfile(dataPath, 'ops.mat'));
ops.parfor = true;

%% artifact detection goes here?

% additional preprocessing steps before running kilosort
% 1) save high-pass filtered version of the binary file for visualization
%    in Phy
% 2) detect large events that span multiple channels and remove them
ops=removeArtifacts(ops);

%% check that the timing lines up
% checkRawBinaryMatch(ops, 10e3)

%% main spike-sorting routine
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)


%% saving
% save matlab results file
fprintf('saving matlab results file\n')
save(fullfile(ops.root,  'rez.mat'), 'rez', 'ops', '-v7.3');

% rez                = merge_posthoc2(rez);
fprintf('saving python files for Phy\n')
% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);