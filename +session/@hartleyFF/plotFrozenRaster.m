function plotFrozenRaster(h, s, clustId)

if nargin < 3
    clustId = s.cids;
end
nUnits = numel(clustId);

hartleyTrial = h.trial;
frozenTrials = find([h.trial.frozenSequence]);
if any(frozenTrials)
    tmp = arrayfun(@(x) x.frameTimes(1:x.frozenSequenceLength:end)', hartleyTrial(frozenTrials), 'UniformOutput', false)';
    tmp = cellfun(@(x) x(:), tmp(:), 'uni', false);
    sequenceStarts = cell2mat(tmp);
    
    cmap = lines(nUnits);
    
    for kUnit = 1:nUnits
        st = s.st(s.clu==clustId(kUnit));
        
        seqLength = hartleyTrial(frozenTrials(1)).frozenSequenceLength;
        
        [spcnt, bcenters] = pdsa.binSpTimes(st, sequenceStarts, [0 seqLength/120], 1e-3);
        
        [i, j] = find(spcnt);
        
        uIx = s.cids==clustId(kUnit);
        
        if isfield(s, 'cgs') && s.cgs(uIx)==3
            plot(bcenters(j), i-numel(sequenceStarts)*(kUnit-1), '.', 'Color', cmap(kUnit,:)); hold on
        else
            plot(bcenters(j), i-numel(sequenceStarts)*(kUnit-1), '.', 'Color', repmat(.5, 1, 3)); hold on
        end
        
    end
end

end