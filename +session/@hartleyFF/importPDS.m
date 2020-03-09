function [trial, display, trialIdx] = importPDS(PDS)

pdsDate = PDS.initialParametersMerged.session.initTime;
if isfield(PDS.initialParametersMerged.git, 'pep')
    
    if any(strfind(PDS.initialParametersMerged.git.pep.status, 'branch cleanup'))
        
        if pdsDate > datenum(2018, 02, 01)
            [trial, display, trialIdx] = session.hartleyFF.importPDS_v2(PDS);
        else
            error('unknown version')
        end
        
    else
        warning('hartleyFF: git tracking failed. running import version 2')
        try
            [trial, display, trialIdx] = session.hartleyFF.importPDS_v2(PDS);
        catch
            error('version 2 import failed')
        end
    end
else
    [trial, display, trialIdx] = session.hartleyFF.importPDS_v1(PDS);
end

end