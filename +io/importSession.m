function importSession(thisSession)
% this is a bottom up import of a session. Assumes nothing has been run
% yet. It will abort if it finds that some

SERVER_DATA_DIR = getpref('EPHYS', 'SERVER_DATA');
disp(thisSession)

if isnan(thisSession.Electrode{1})
    disp('No electrode data. skipping.')
	return
end

shank = hardware.electrodeFactory(thisSession.Electrode{1});

oepath = fullfile(SERVER_DATA_DIR, thisSession.Directory{1});

% check if import has already been run
derived = dir(fullfile(oepath, '_*'));
try
    for iDir = 1:numel(derived)
        ops_old = io.loadOps(oepath);
        
        for i = 1:numel(ops_old)
            o = ops_old(i);
            
            o = io.convertOpsToNewDirectory(o, oepath);
            
            save(fullfile(o.root, 'ops.mat'), '-v7.3', '-struct', 'o')
        end
    end
end

io.oe2dat(oepath, shank, 'overwrite', false, 'verbose', true);

thisSession.oe2dat = true;
% --- load session: This is central to all operations
sess = io.loadSession(oepath);



overwrite = false; % if it breaks, run again with true
PDS = io.getPds(sess, overwrite);


io.getEdf(sess, PDS, overwrite);


%% Step 4: Import LFP
% Loop over sessions and import the local field potential
% This is really slow because of the line-noise fitting and correction
ops = io.loadOps(oepath);

plotIt = false;
for i = 1:numel(ops)
    [~, ~, info] = io.getLFP(ops(i),overwrite,plotIt);
end

if info.phaseCorrection
    thisSession.LfpPhaseCorrection = true;
end

io.writeMeta(thisSession)
