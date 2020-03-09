function [data, timestamps] = loadAndPreprocess(ops, inds)

if nargout>1
    [data, timestamps] = io.loadRaw(ops, inds);
else
    data = io.loadRaw(ops, inds);
end

% notch filter
if ops.NotchFilter60
    wo = 60/(ops.fs/2);  bw = wo/35;
    [b,a] = iirnotch(wo,bw);
    

    datr = data';

    datr = filter(b, a, datr);
    datr = flipud(datr);
    datr = filter(b, a, datr);
    datr = flipud(datr);

    data = datr';
end

switch ops.softwareReferencing
    case 'commonaverage'
        data     = bsxfun(@minus, data, mean(data));
    otherwise
end

% --- bandpass filtering
if isfield(ops, 'fshigh')
    if isfield(ops, 'fslow')
        data = preprocess.highpass(data, ops.fs, ops.fshigh, ops.fslow, 3);
    else
        data = preprocess.highpass(data, ops.fs, ops.fshigh, 0, 3);
    end
end

if isfield(ops, 'artifactThresh') && ~isinf(ops.artifactThresh)
    [data, bad] = preprocess.removeChannelArtifacts(data, ops.artifactThresh, ops.artifactNchans, 50); 
end

% --- apply channel map
if exist(ops.chanMap, 'file')
    load(ops.chanMap)
    data = data(chanMap(connected),:);
end

% --- convert to milivolts
data = data * ops.bitVolts;