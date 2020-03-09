function KiloUndoAutomerge(ops)

load(fullfile(ops.root,  'rez.mat'))

fprintf('saving python files for Phy\n')

% save python results file for Phy
rezToPhy(rez, ops.root);