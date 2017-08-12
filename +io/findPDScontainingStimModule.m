function hasStim = findPDScontainingStimModule(PDS, stim)
% hasStim = findPDScontainingStimModule(PDS, stim)

hasStim = false(numel(PDS),1);
for i = 1:numel(PDS)
    if isfield(PDS{i}.initialParametersMerged, stim)
        ixcond = cellfun(@(x) isfield(x, stim), PDS{i}.conditions);
        ixdata = cellfun(@(x) isfield(x, stim), PDS{i}.data);
        if any(ixcond)
            ixuse = cellfun(@(x) isfield(x.(stim), 'use'), PDS{i}.conditions(ixcond));
            inds  = find(ixcond);
            hasStim(i) = any(cellfun(@(x) x.(stim).use, PDS{i}.conditions(inds(ixuse)))) | PDS{i}.initialParametersMerged.(stim).use;
        elseif any(ixdata)
            ixuse = cellfun(@(x) isfield(x.(stim), 'use'), PDS{i}.data(ixdata));
            inds  = find(ixdata);
            hasStim(i) = any(cellfun(@(x) x.(stim).use, PDS{i}.conditions(inds(ixuse)))) | PDS{i}.initialParametersMerged.(stim).use;
        end
    end
end