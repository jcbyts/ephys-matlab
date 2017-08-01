function KiloAutomerge(ops)
% KiloAutomerge(ops)

load(fullfile(ops.root,  'rez.mat'))

fprintf('Attempting Automerge\n')
rez                = merge_posthoc2(rez); %#ok<NODEF>

fprintf('saving python files for Phy\n')

% save python results file for Phy
rezToPhy(rez, ops.root);