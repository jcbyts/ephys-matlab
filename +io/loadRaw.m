function [data, timestamps] = loadRaw(ops, inds, inVolts, fproc)
% [data, timestamps] = loadRaw(ops, inds, inVolts, fproc)

if nargin < 4
    fproc = false;
end

if nargin < 3
    inVolts = false;
end

if fproc
    fid = fopen(ops.fproc);
else
    fid = fopen(ops.fbinary);
end

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

data = fread(fid, bufferSize, '*int16');

if nargout > 1 || inVolts
   
    info = load(fullfile(ops.root, 'ephys_info.mat'));
    if inVolts
        data = double(data)*info.bitVolts;
        
        
%         if isfield(info, 'artifacts')
%             for iCh = 1:size(data,1)
%                 artifacts = info.artifacts - inds(1);
%                 artifacts = artifacts(artifacts > 0 & artifacts < (inds(end) - inds(1)));
%                 data(iCh,artifacts) = nan;
%             end
%         end
    end
    
    timestamps = io.convertSamplesToTime(inds, info.sampleRate, info.timestamps(:), info.fragments(:));
end

fclose(fid);
