function W = getSpikeWaveformsFromOps(ops, sp, varargin)
% PLOT WAVEFORMS
% Inputs:
%   ops@struct - ops struct
%   sp@struct  - spikes struct
% Optional Arguments (as pairs)
%   'numWaveforms'    - default 500
%   'postSpikeBuffer' - number of samples after spike
%   'preSpikeBuffer'  - number of samples before spike
%   'clusterIds'      - list of clusters to plot
%   'useMean'         - plot mean +- SD instead of all waveforms
% Output:
%   wfs - figure handle

% 2017.08.14    jly     wrote it
ip = inputParser();
ip.addOptional('numWaveforms', 500)
ip.addOptional('postSpikeBuffer', 40)
ip.addOptional('preSpikeBuffer', 10)
ip.addOptional('clusterIds', [])
ip.addOptional('useMean', true)
ip.parse(varargin{:})

% useMean = ip.Results.useMean;
% 
if ~exist(ops.fproc, 'file')
    warning('spikeWaveformsFromOps: No high-pass filtered data. Using wideband. This is bad. You should run preproccess.save_filtered_data(ops)')
    ops.fproc = ops.fbinary; % point to wideband data instead
%     useMean = false;
end

if isempty(ip.Results.clusterIds)
    if isfield(sp, 'clusterDepths')
        [~, depthIdx] = sort(sp.clusterDepths);
    else
        depthIdx = 1:numel(sp.cids);
    end
    
    clustId = sp.cids(depthIdx);
else
    clustId = ip.Results.clusterIds;
end

% open file to read raw data
fname = ops.fproc;
fid   = fopen(fname, 'r');
buffer = [ops.Nchan ip.Results.postSpikeBuffer+ip.Results.preSpikeBuffer];

if exist(ops.chanMap, 'file')
    load(ops.chanMap)
else
    chanMap = 1:ops.Nchan;
end

W = struct();
W.cids = clustId;
W.wftax = (-ip.Results.preSpikeBuffer:ip.Results.postSpikeBuffer)/ops.fs;
W.chx   = sp.xc;   
W.chy   = sp.yc;
W.wf = [];
W.sd = [];

for kClust = 1:numel(clustId)
    
    iix = sp.clu==clustId(kClust);
    n = sum(iix);
    ss = sp.ss(iix);
    ss = ss - ip.Results.preSpikeBuffer;
    
    numWaveforms = ip.Results.numWaveforms;
    
    n = min(n, numWaveforms);
    
    wfs = zeros(n, prod(buffer));
    
    for i = 1:n
        fseek(fid, 2*ss(i)*ops.Nchan, 'bof');
        data = double(fread(fid, buffer, '*int16'));
        data = data(chanMap,:)';
        
        wfs(i,:) = data(:)*ops.bitVolts;
    end
    
    
    wf = median(wfs);
   
    ci = prctile(wfs, [15.8650   84.1350]);
    sd = mean(abs(bsxfun(@minus, ci, wf)));
    
    W.wf = [W.wf wf(:)];
    W.sd = [W.sd sd(:)];
    
%     figure(1); clf
%     plot(wf); hold on
%     plot(ci(1,:))
%     plot(ci(2,:))
end

fclose(fid);