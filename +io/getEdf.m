function [data, timestamps, elInfo] = getEdf(ops, PDS, overwrite)
% [data, timestamps, elInfo] = getEdf(ops, PDS)

if nargin < 3
    overwrite = false;
    if nargin < 2
        PDS = io.getPds(ops);
    end
end

elInfo.timestamps = [];
elInfo.fragments  = [];
elInfo.sampleRate = [];
elInfo.dateNum    = [];
elInfo.fields     = {'EyeX', 'EyeY', 'PupilArea'};
elInfo.bitDeg     = 1e-3;

feye    = fullfile(ops.root, 'eyepos.dat');
felinfo = fullfile(ops.root, 'eye_info.mat');

if exist(feye, 'file') && ~overwrite
   elInfo = load(felinfo);
   
   fid = fopen(feye, 'r');
   fseek(fid, 0, 'eof');
   filesize = ftell(fid);
   fseek(fid, 0, 'bof');
   
   nFields = numel(elInfo.fields);
   nTotSamps = filesize/nFields/2;
   
   data = double(fread(fid, [nFields nTotSamps], '*uint16'));
   data(1:2,:) = data(1:2,:)*elInfo.bitDeg(1) + elInfo.bitDeg(2);
   data(1,data(3,:)==0) = nan;
   data(2,data(3,:)==0) = nan;
   
   timestamps = io.convertSamplesToTime(1:nTotSamps, elInfo.sampleRate, elInfo.timestamps(:), elInfo.fragments(:));
   
    return
    
end

fidout = fopen(feye, 'w');

nPds = numel(PDS);

for kPds = 1:nPds

    [dataRAW, ~, el_info] = getEdfData(ops, PDS{kPds});
    fwrite(fidout, dataRAW', '*uint16')
    
    elInfo.timestamps   = [elInfo.timestamps el_info.timestamps];
    elInfo.fragments    = [elInfo.fragments el_info.fragments];
    elInfo.sampleRate   = [elInfo.sampleRate el_info.sampleRate];
    elInfo.dateNum      = [elInfo.dateNum el_info.dateNum];
    elInfo.bitDeg       = el_info.bitDeg;
end

elInfo.sampleRate = unique(elInfo.sampleRate);
assert(numel(elInfo.sampleRate)==1, 'Eyelink sample rate changed throughout session. This is not supported')

fclose(fidout);
save(felinfo, '-v7.3', '-struct', 'elInfo');

fid = fopen(feye, 'r');
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid, 0, 'bof');

nFields = numel(elInfo.fields);
nTotSamps = filesize/nFields/2;

data = double(fread(fid, [nFields nTotSamps], '*uint16'));
data(1:2,:) = data(1:2,:)*elInfo.bitDeg(1) + elInfo.bitDeg(2);
data(1,data(3,:)==0) = nan;
data(2,data(3,:)==0) = nan;

timestamps = io.convertSamplesToTime(1:nTotSamps, elInfo.sampleRate, elInfo.timestamps(:), elInfo.fragments(:));




end

function [data, timestamps, info] = getEdfData(ops, PDS)



assert(PDS.initialParametersMerged.eyelink.useRawData, 'This import is only designed for raw data + calibration matrix')

% Find all times the calibration matrix changed
eyelinkChanged = find(cellfun(@(x) isfield(x, 'eyelink'), PDS.initialParameters));
pauseTrials    = find(cellfun(@(x) ~isempty(regexp(x, 'Pause', 'once')), PDS.initialParameterNames));
eyelinkChangedAfterPause = intersect(eyelinkChanged, pauseTrials);
changedParamNames = PDS.initialParameterNames(eyelinkChangedAfterPause);
trialNums = cellfun(@(x) str2double(cell2mat(regexp(x, '\d+', 'match'))), changedParamNames);

calibMatChanged = cellfun(@(x) isfield(x.eyelink, 'calibration_matrix'), PDS.initialParameters(eyelinkChangedAfterPause));

calibMat = [];
if ~isempty(PDS.initialParametersMerged.eyelink.calibration_matrix)
    calibMat = {PDS.initialParametersMerged.eyelink.calibration_matrix};
end

calibMat = [calibMat cellfun(@(x) x.eyelink.calibration_matrix, PDS.initialParameters(eyelinkChangedAfterPause), 'UniformOutput', false)];

calibMatChangeIdx = trialNums(calibMatChanged);

ElTimeCmChanged = [0 cellfun(@(x) x.timing.eyelinkStartTime(2), PDS.data(calibMatChangeIdx))]*1e3;

eyeIdx = PDS.initialParametersMerged.eyelink.eyeIdx;


[~, edfFile, ~] = fileparts(PDS.initialParametersMerged.session.file);

elFile=fullfile(ops.root, [edfFile '.edf']);
hasEdf=exist(elFile, 'file');
if hasEdf
    el=edfmex(elFile);
    elTrialStarts=arrayfun(@(x) strcmp(x.message, 'TRIALSTART'), el.FEVENT);
    elTrialTS=double(arrayfun(@(x) x.sttime, el.FEVENT(elTrialStarts)));
end

pdsTrialStarts=cellfun(@(x) x.timing.eyelinkStartTime(1), PDS.data);

EL2OEfit =[elTrialTS(:) ones(numel(elTrialTS),1)]\PDS.PTB2OE(pdsTrialStarts(:));
EL2OE  =@(x) x*EL2OEfit(1) + EL2OEfit(2); % this converts from Eyelink sample times to OE times

edfTime=double(el.FSAMPLE.time);

x=double(el.FSAMPLE.px);
y=double(el.FSAMPLE.py);

numCmChanges = numel(calibMat);

% build raw input
X=[x(eyeIdx,:); y(eyeIdx,:); ones(1,size(y,2))];

eyeXyDeg = zeros(size(X,2), 2);

for iCalib = 1:numCmChanges
    
    cm=calibMat{iCalib}(:,:,eyeIdx);
    
    if iCalib == numCmChanges
        iix = edfTime > ElTimeCmChanged(iCalib);
    else
        iix = edfTime > ElTimeCmChanged(iCalib) & edfTime < ElTimeCmChanged(iCalib+1);
    end
    
    eyePxTmp = cm*X(:,iix); % convert raw to pixels
    
    win = PDS.initialParametersMerged.display.winRect(3:4);
    % nan times when the eye was off the screen
    offScreen = (eyePxTmp(1,:)>=win(1) | eyePxTmp(1,:)<=0) | (eyePxTmp(2,:)>=win(2) | eyePxTmp(2,:)<=0);
    
    % center
    eyePxTmp = bsxfun(@minus, eyePxTmp, PDS.initialParametersMerged.display.ctr(1:2)');
    
    eyeDegTmp = pds.px2deg(eyePxTmp, PDS.initialParametersMerged.display.viewdist, PDS.initialParametersMerged.display.px2w);
    eyeDegTmp(:,offScreen) = nan;
    
    eyeXyDeg(iix, :) = eyeDegTmp';
end

dT = diff(edfTime);
ss = mode(dT);

info.sampleRate = 1e3/ss;

tmp = pds.px2deg(-PDS.initialParametersMerged.display.ctr(1:2)', PDS.initialParametersMerged.display.viewdist, PDS.initialParametersMerged.display.px2w);
tmp = min(tmp);

info.bitDeg = [1/100 tmp];

eyeXyDeg(:,2) = eyeXyDeg(:,2)*-1;

data = [int16((eyeXyDeg - tmp)*100) el.FSAMPLE.pa(eyeIdx,:)'];

breaks = find(dT ~= ss);

info.fragments = diff([0 breaks numel(edfTime)]);
timestamps = EL2OE(edfTime);
info.timestamps = timestamps([1 breaks]);
info.dateNum = PDS.initialParametersMerged.session.initTime;
end
    