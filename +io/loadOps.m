function ops = loadOps(fpath)

warning('off')
ops = load(fullfile(fpath, 'ops.mat'));

if isfield(ops, 'parfor_')
    ops.parfor = ops.parfor_;
    ops = rmfield(ops, 'parfor_');
end
