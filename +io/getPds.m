function PDS = getPds(sessionInfo)

pdsList = io.getPdsFileList(sessionInfo);

evList = io.getEventsFileList(sessionInfo);

nPdsFiles = numel(pdsList);

PDS = cell(nPdsFiles,1);
fprintf('Found [%d] PDS files\n', nPdsFiles)
fprintf('Aligning with OE clock now\n')
for kPdsFile = 1:nPdsFiles
    
    tmp = load(pdsList{kPdsFile}, '-mat');
    if ~isempty(tmp.PDS.functionHandles)
        fhlist = fieldnames(tmp.PDS.functionHandles);
        for i = 1:numel(fhlist)
            if isa(tmp.PDS.functionHandles.(fhlist{i}), 'matlab.ui.Figure')
                close(tmp.PDS.functionHandles.(fhlist{i}))
            end
        end
    end
    [OE2PTBfit, ~,PTB2OE, maxreconstructionerror] = io.sync2OeClock(tmp.PDS, evList);
    
    tmp.PDS.PTB2OE    = PTB2OE;
    tmp.PDS.OE2PTBfit = OE2PTBfit;
    tmp.PDS.maxreconstructionerror = maxreconstructionerror;
    
    if isempty(maxreconstructionerror)
        fprintf('%d) No Ephys Data\n', kPdsFile);
    else
        fprintf('%d) Strobe times aligned. Max reconstruction error is %2.3f ms\n', kPdsFile, maxreconstructionerror*1e3)
    end
    
    PDS{kPdsFile} = tmp.PDS;
end

noEphys = cellfun(@(x) isempty(x.maxreconstructionerror), PDS);

PDS(noEphys) = [];
