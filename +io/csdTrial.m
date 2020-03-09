function csdTrial = csdTrial(PDS, varargin)
% build CSD trial struct
% Inputs:
%   PDS@cell - array of PDS structs
% Output:
%   csdTrial@struct - array of trials
% Example call:
%   csdTrial = csdTrial(PDS)


stim = 'csdFlash';

hasStim = io.findPDScontainingStimModule(PDS, stim);


csdTrial = struct();
trialNum = 0;

for i = find(hasStim(:)')
    
    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 

    stimTrials = find(trialIx);
    
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
    
        kTrial = trialNum + j;
        
        csdTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
        csdTrial(kTrial).start      = csdTrial(kTrial).frameTimes(1);
        csdTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - csdTrial(kTrial).start;
    
        csdTrial(kTrial).on         = PDS{i}.data{thisTrial}.(stim).on;
        csdTrial(kTrial).onset      = csdTrial(kTrial).frameTimes(diff(csdTrial(kTrial).on)==1);
    end
    
    trialNum = kTrial;
    
end
