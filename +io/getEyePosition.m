function eyepos = getEyePosition(PDS, kTrial)


cm = PDS.initialParametersMerged.eyelink.calibration_matrix;

eyeIdx = PDS.initialParametersMerged.eyelink.eyeIdx;
useRaw = PDS.initialParametersMerged.eyelink.useRawData;

% --- Get eye position from the eyelink buffer
sampleNames=PDS.data{kTrial}.eyelink.sampleIds;

if useRaw
    ixEyePos=(cellfun(@(x) any(strfind(x, 'EyeRaw')), sampleNames)); % index to find Eye position
    ptbElXY   = PDS.data{kTrial}.eyelink.samples(ixEyePos,:);
    XY = (cm(:,:,eyeIdx)*[ptbElXY; ones(1,size(ptbElXY,2))])';
else 
	ixEyePos=(cellfun(@(x) any(strfind(x, 'EyeX')), sampleNames)) | (cellfun(@(x) any(strfind(x, 'EyeY')), sampleNames)); % index to find Eye position
    XY = PDS.data{kTrial}.eyelink.samples(ixEyePos,:)';
end
ixPupil=cellfun(@(x) any(strfind(x, 'Pupil')), sampleNames);



ptbElTime = PDS.data{kTrial}.eyelink.samples(1,:);
       
        
time = ptbElTime/1e3-PDS.data{kTrial}.timing.eyelinkStartTime(2);



pupil = PDS.data{kTrial}.eyelink.samples(ixPupil,:)';

eyepos = [time(:) XY pupil(:)];

