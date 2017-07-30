function [sessionInfo, ops, info] = loadSession(oepath)
% [sessionInfo, ops, info] = loadSession(oepath)

ops = io.loadOps(oepath);

info = load(fullfile(ops.root, 'ephys_info.mat'));

load(fullfile(oepath, 'session_info.mat'))