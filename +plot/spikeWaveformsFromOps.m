function fig = spikeWaveformsFromOps(ops, sp, varargin)
% PLOT WAVEFORMS
% Inputs:
%   ops@struct - ops struct
%   sp@struct  - spikes struct
% Optional Arguments (as pairs)
%   'figure'          - figure number or handle
%   'numWaveforms'    - default 500
%   'postSpikeBuffer' - numbe of samples after spike to plot
%   'preSpikeBuffer'  - numbe of samples after spike to plot
%   'clusterIds'      - list of clusters to plot
% Output:
%   fig - figure handle
% Example call:
%   fig = plot.spikeWaveformsFromOps(ops, sp, varargin)

% 2017.08.14    jly     wrote it

ip = inputParser();
ip.addOptional('figure', [])
ip.addOptional('numWaveforms', 500)
ip.addOptional('postSpikeBuffer', 40)
ip.addOptional('preSpikeBuffer', 10)
ip.addOptional('clusterIds', [])
ip.parse(varargin{:})

if isempty(ip.Results.figure)
    fig = gca;
else
    fig = figure(ip.Results.figure);
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

cmap = lines;

% open file to read raw data
fname = ops.fbinary;
fid   = fopen(fname, 'r');
buffer = [ops.Nchan ip.Results.postSpikeBuffer+ip.Results.preSpikeBuffer];

if exist(ops.chanMap, 'file')
    load(ops.chanMap)
else
    chanMap = 1:ops.Nchan;
end

for kClust = 1:numel(clustId)
    
iix = sp.clu==clustId(kClust);
n = sum(iix);
ss = sp.ss(iix);

for i = ss(1:ceil((n/ip.Results.numWaveforms)):n)'
       
    fseek(fid, (i-10)*2*ops.Nchan, 'bof');
    data = double(fread(fid, buffer, '*int16'));
    data = data(chanMap,:)';
    data = bsxfun(@minus, data, mean(data([1:10 (buffer(2)-10):buffer(2)],:)));
    if isfield(sp, 'yc')
        wf = bsxfun(@plus, data*ops.bitVolts, 5*flipud(sp.yc)');
    else
        wf = data*ops.bitVolts;
    end
    plot((1:buffer(2))+kClust*buffer(2), wf, 'Color', cmap(kClust,:)); hold on
end
drawnow
end
axis tight

fclose(fid);