function [data, timestamps] = loadRaw(ops, inds, inVolts)
% [data, timestamps] = loadRaw(ops, inds)

if nargin < 3
    inVolts = false;
end

fid = fopen(ops.fbinary);

fseek(fid, 0, 'eof');
filesize = ftell(fid);

fseek(fid, 0, 'bof');

Nchan = ops.Nchan;

nTotSamp = filesize/Nchan;
if ~exist('inds', 'var') || isempty(inds)
    inds = 1:nTotSamp;
    bufferSize = [Nchan nTotSamp];
else
    fseek(fid, 2*inds(1)*Nchan, 'bof');
    bufferSize = [Nchan inds(end)-(inds(1)-1)];
end

data = double(fread(fid, bufferSize, '*int16'));

if nargout > 1 || inVolts
   
    info = load(fullfile(ops.root, 'ephys_info.mat'));
    if inVolts
        data = double(data)*info.bitVolts;
    end
    
    timestamps = io.convertSamplesToTime(inds, info.sampleRate, info.timestamps(:), info.fragments(:));
end

fclose(fid);
