function KiloAutomerge(ops)
% KiloAutomerge(ops)

tmp = load(fullfile(ops.root,  'rez.mat'));
rez = tmp.rez;

fprintf('Attempting Automerge\n')
rez                = merge_posthoc3(rez);

fprintf('saving python files for Phy\n')

% save python results file for Phy
rezToPhy(rez, ops.root);