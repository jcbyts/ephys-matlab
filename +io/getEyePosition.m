function eyepos = getEyePosition(PDS, kTrial)

if PDS.initialParametersMerged.eyelink.use == 1
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
        
%***************************************************************
%*****% ADD by Shanna to take in arrington info*****************
%***************************************************************
elseif PDS.initialParametersMerged.arrington.use == 1 
           
    cm = PDS.initialParametersMerged.arrington.calibration_matrix;
    eyeIdx = PDS.initialParametersMerged.arrington.eyeIdx;
    %********* I don't know what to do about this here, useRaw or not
    useRaw = PDS.initialParametersMerged.eyelink.useRawData;
    
    k = kTrial;
    iFrame = PDS.data{k}.iFrame;
    times = (PDS.data{k}.timing.flipTimes(1,:) - PDS.data{k}.timing.arringtonStartTime(1))';
    eye_xy = PDS.data{k}.faceforage.eyes(1:iFrame,1:2);
    pupil = zeros(iFrame,1);   % don't have pupil from arrington
    if (useRaw == 1)
               eye_xy = (cm(:,:,1)*[eye_xy'; ones(1,size(eye_xy',2))])';
    end
    eyepos = [times eye_xy pupil];
end
