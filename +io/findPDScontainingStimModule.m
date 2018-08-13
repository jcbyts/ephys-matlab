function hasStim = findPDScontainingStimModule(PDS, stim)
% FIND PDS CONTATINING STIM MODULE
% Takes in a cell aray of PDS data structs and returns which one contain a particular stimulus module stim
% Inputs:
%   PDS@cell    - cell array of PDS 
%   stim@string - name of module to recover
% Output:
%   hasStim@logical - boolean array of which PDS have the stim
% Example call:
%   hasStim = findPDScontainingStimModule(PDS, stim)

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
            if isempty(inds) || ~any(ixuse)
                continue
            end
            hasStim(i) = any(cellfun(@(x) x.(stim).use, PDS{i}.conditions(inds(ixuse)))) | PDS{i}.initialParametersMerged.(stim).use;
        end
    end
end