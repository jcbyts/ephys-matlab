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
%   'useMean'         - plot mean +- SD instead of all waveforms
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
ip.addOptional('useMean', true)
ip.addOptional('yscale', 150)
ip.parse(varargin{:})

useMean = ip.Results.useMean;

if ~exist(ops.fproc, 'file')
    warning('spikeWaveformsFromOps: No high-pass filtered data. Using wideband')
    ops.fproc = ops.fbinary; % point to wideband data instead
    useMean = false;
end
    
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
fname = ops.fproc;
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
    xax = (1:buffer(2))+kClust*buffer(2);
    
    numWaveforms = ip.Results.numWaveforms;
    
    if isfield(sp, 'wfs')
        xax = sp.wftax + kClust*diff(sp.wftax([1 end]));
        offset = ip.Results.yscale*reshape(repmat((1:numel(sp.xc)), numel(xax), 1), [], 1)';
        wfs = reshape(bsxfun(@plus, sp.wfs(iix,:), offset), sum(iix), numel(xax), numel(sp.xc));
        
        if useMean
            mwf = squeeze(mean(wfs));
            plot(xax, mwf, 'Color', cmap(kClust,:)); hold on
            sd  = squeeze(std(wfs));
            
            plot(xax, mwf+sd, ':', 'Color', cmap(kClust,:));
            plot(xax, mwf-sd, ':', 'Color', cmap(kClust,:));
            text(xax(2), max(mwf(:) + sd(:)) + ip.Results.yscale/2, sprintf('unit: %d', kClust), 'Color', cmap(kClust,:))
        else
            
            
        end
        
    else
        
        
        if useMean
            wfs = zeros(numWaveforms,numel(xax), numel(sp.xc));
        end
        
        wctr = 1;
        for i = ss(1:ceil((n/numWaveforms)):n)'
            
            fseek(fid, (i-10)*2*ops.Nchan, 'bof');
            data = double(fread(fid, buffer, '*int16'));
            data = data(chanMap,:)';
            data = bsxfun(@minus, data, mean(data([1:10 (buffer(2)-10):buffer(2)],:)));
            if isfield(sp, 'yc')
                wf = bsxfun(@plus, data*ops.bitVolts, 5*flipud(sp.yc)');
            else
                wf = data*ops.bitVolts;
            end
            
            if useMean
                wfs(wctr,:,:) = wf;
                wctr = wctr + 1;
            else
                plot(xax, wf, 'Color', cmap(kClust,:)); hold on
            end
        end
        
        if useMean
            mwf = squeeze(mean(wfs));
            plot(xax, mwf, 'Color', cmap(kClust,:)); hold on
            
            %         th = mode(diff(mwf(1,:)))/2;
            sd  = squeeze(std(wfs));
            %         sd(sd > th) = nan;
            plot(xax, mwf+sd, ':', 'Color', cmap(kClust,:));
            plot(xax, mwf-sd, ':', 'Color', cmap(kClust,:));
        end
        
        drawnow
    end
    
end
axis tight

fclose(fid);