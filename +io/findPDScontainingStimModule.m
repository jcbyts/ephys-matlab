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
    
    trial = pds.getPdsTrialData(PDS{i});
    
    if isfield(trial, stim)
        hasStim(i) = any(arrayfun(@(x) x.(stim).use, trial));
    end
end