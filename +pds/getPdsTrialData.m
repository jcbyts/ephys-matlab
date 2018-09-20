function data = getPdsTrialData(PDS)
% getPdsTrialData merges the conditions and data and outputs the parameters
% and data for each trial. 
% The hierarchy of parameters is like this:
% data -> conditions -> initialParameters 

if ~isempty(PDS.conditions)    
    A = cellfun(@(x) mergeStruct(PDS.initialParametersMerged, x), PDS.conditions, 'uni', 1);
else
    A = repmat({PDS.initialParametersMerged}, 1, numel(PDS.data));
end

data = cellfun(@(x,y) mergeStruct(x, y), A, PDS.data, 'uni', 1);