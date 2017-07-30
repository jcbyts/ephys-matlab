function ops = oe2dat(oepath)

% --- declare constants
NUM_HEADER_BYTES    = 1024; % size of header for each analog file
SAMPLES_PER_RECORD  = 1024;
BLOCK_BYTES         = 2070; % see load_open_ephys_data_faster

fephys  = fullfile(oepath, 'ephys.dat');
fops    = fullfile(oepath, 'ops.mat');
finfo   = fullfile(oepath, 'ephys_info.mat');

% --- get channel info
fl  = dir(fullfile(oepath, '*CH*.continuous'));
str = 'CH(\d+)';
chlist = cellfun(@(x) cell2mat(regexp(x, str, 'match')) , {fl.name}, 'UniformOutput', false);

channels    = unique(cellfun(@(x) str2double(x(3:end)), chlist));
numChannels = numel(channels);

fs = cell(numChannels,1);
for j = 1:numChannels
    tmp = dir(fullfile(oepath, sprintf('*CH%d.continuous',j) ));
    % if separate files are saved by open ephys gui
    tmp_ = dir(fullfile(oepath, sprintf('*CH%d_*.continuous', j) ));
    
    fs{j} = [tmp tmp_];
end

% --- convert raw data
fidout      = fopen(fephys, 'w');

nblocks = cellfun(@(x) numel(x), fs);
if numel(unique(nblocks))>1
    error('different number of blocks for different channels!')
end

fid = cell(numChannels, 1);

nBlocks     = unique(nblocks);
nSamples    = SAMPLES_PER_RECORD;  % fixed to 1024 for now!
fprintf('Found [%d] total blocks\n', nBlocks)

% --- initialize meta data
info = struct('format', [], ...
    'version', [], ...
    'dateNum', [], ...
    'sampleRate', [], ...
    'bitVolts', [], ...
    'timestamps', [], ...
    'fragments', []);

tic
for k = 1:nBlocks
    fprintf('Block %d\n', k)
    
    for j = 1:numChannels
        fid{j}             = fopen(fullfile(fs{j}(k).folder, fs{j}(k).name));
        
        % --- Check that OE version is correct
        fseek(fid{j},0,'bof');
        hdr = fread(fid{j}, NUM_HEADER_BYTES, 'char*1');
        info_ = getHeader(hdr);
        assert(info_.header.version == 0.4, 'open ephys version has changed. You must now check that everything still works. Sorry, Jake')
        
        % save header info
        if j==1 && k==1
            info.format     = info_.header.format;
            info.version    = info_.header.version;
            info.sampleRate = info_.header.sampleRate;
            info.bitVolts   = info_.header.bitVolts;
        else
            assert(strcmp(info.format,info_.header.format), 'FORMAT does not match across blocks')
            assert(info.version == info_.header.version, 'VERSION does not match across blocks')
            assert(info.sampleRate == info_.header.sampleRate, 'SAMPLERATE does not match across blocks')
            assert(info.bitVolts == info_.header.bitVolts, 'BITVOLTS does not match across blocks')
        end
        
        % --- get timestamps from this file
        if j==1
            % get file size
            fseek(fid{j},0,'eof');
            filesize  = ftell(fid{j});
            numIdx = floor((filesize - NUM_HEADER_BYTES)/BLOCK_BYTES);
            
            % skip header again
            fseek(fid{j}, NUM_HEADER_BYTES, 'bof');
            
            % read timestamps
            seg = fread(fid{j}, numIdx, '1*int64', BLOCK_BYTES - 8, 'l');
            
            % save timestamps of all starts and stops
            ts_ = seg(1);
            ns_ = numel(seg)*SAMPLES_PER_RECORD;
            stops = find(diff(seg)~=SAMPLES_PER_RECORD)+1;
            
            if any(stops)
                ts_ = [ts_ seg(stops)'];
                ns_ = diff([0 stops(:)' numel(seg)])*SAMPLES_PER_RECORD;
                
%                 fprintf('Found multiple paused segments\n')
%                 figure
%                 plot(seg, '.')
%                 hold on
%                 plot(cumsum([[0 ns_(1:end-1)]; ns_], 2)/SAMPLES_PER_RECORD, [ts_; ts_])
            end
            
            nFragments = numel(ns_);
            assert(nFragments==numel(ts_), 'oe2dat: TIMESTAMPS and FRAGMENTS are different sizes')
            
            % get datetime of recording start
            dn_ =datenum(info_.header.date_created, 'dd-mmm-yyyy HHMMSS');
            
            info.dateNum = [info.dateNum repmat(dn_,1,nFragments)];
            info.timestamps = [info.timestamps ts_];
            info.fragments  = [info.fragments ns_];
            
            % --- skip header
            fseek(fid{j}, NUM_HEADER_BYTES, 0);
        end
        
        
    end
    
    
    % --- save raw data to new binary file
    nsamps = 0;
    flag = 1;
    while 1
        samples = zeros(nSamples * 1000, numChannels, 'int16');
        for j = 1:numChannels
            collectSamps    = zeros(nSamples * 1000, 1, 'int16');
            
            rawData         = fread(fid{j}, 1000 * (nSamples + 6), '1030*int16', 10, 'b');
            
            nbatches        = ceil(numel(rawData)/(nSamples+6));
            for s = 1:nbatches
                rawSamps = rawData((s-1) * (nSamples + 6) +6+ [1:nSamples]);
                collectSamps((s-1)*nSamples + [1:nSamples]) = rawSamps;
            end
            samples(:,j)         = collectSamps;
        end
        
        if nbatches<1000
            flag = 0;
        end
        if flag==0
            samples = samples(1:s*nSamples, :);
        end
        
        samples         = samples';
        fwrite(fidout, samples, 'int16');
        
        nsamps = nsamps + size(samples,2);
        
        if flag==0
            break;
        end
    end
    nSamplesBlocks(k) = nsamps; %#ok<AGROW>
    
    for j = 1:numChannels
        fclose(fid{j});
    end
    
end

fclose(fidout);

save(finfo, '-v7.3', '-struct', 'info')

fprintf('Done [%02.2fs]\n', toc)



% --- build ops
ops = struct();
ops.root    = oepath;
ops.fbinary = fullfile(ops.root, 'ephys.dat');
ops.Nchan   = numChannels;
ops.nSamplesBlocks = nSamplesBlocks;
ops.bitVolts = info.bitVolts;
ops.fs =info.sampleRate;

% --- default parameters

ops.datatype            = 'dat';  % binary ('dat', 'bin') or 'openEphys'
ops.fproc               = fullfile(oepath, 'tempWh.dat'); % residual from RAM of preprocessed data				
% define the channel map as a filename (string) or simply an array
ops.chanMap             = fullfile(oepath, 'chanMap.mat'); % make this file using createChannelMapFile.m
if ~exist(ops.chanMap, 'file')
	ops.chanMap = 1:ops.Nchan; % treated as linear probe if unavailable chanMap file
end

ops.Nfilt               = 2*(ops.Nchan+32-mod(ops.Nchan,32));  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)     		
ops.nNeighPC            = 12; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)		
ops.nNeigh              = 16; % visualization only (Phy): number of neighboring templates to retain projections of (16)		
		
% options for channel whitening		
ops.whitening           = 'noSpikes'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)		
ops.nSkipCov            = 1; % compute whitening matrix from every N-th batch (1)		
ops.whiteningRange      = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)		
		
ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info). 		

% other options for controlling the model and optimization		
ops.Nrank               = 3;    % matrix rank of spike template model (3)		
ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)		
ops.maxFR               = 20000;  % maximum number of spikes to extract per batch (20000)		
ops.NotchFilter60       = false;
ops.softwareReferencing = 'none';
% ops.fshigh              = 200;   % frequency for high pass filtering		
% ops.fslow             = 2000;   % frequency for low pass filtering (optional)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection		
ops.scaleproc           = 200;   % int16 scaling of whitened data		
ops.NT                  = 128*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory) 		
% for GPU should be multiple of 32 + ntbuff		
		
% the following options can improve/deteriorate results. 		
% when multiple values are provided for an option, the first two are beginning and ending anneal values, 		
% the third is the value used in the final pass. 		
ops.Th               = [6 12 12];    % threshold for detecting spikes on template-filtered data ([6 12 12])		
ops.lam              = [10 30 30];   % large means amplitudes are forced around the mean ([10 30 30])		
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)		
ops.momentum         = 1./[20 1000];  % start with high momentum and anneal (1./[20 1000])		
ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)		
ops.mergeT           = .1;           % upper threshold for merging (.1)		
ops.splitT           = .1;           % lower threshold for splitting (.1)		
		
% options for initializing spikes from data		
ops.initialize      = 'fromData';    %'fromData' or 'no'		
ops.spkTh           = -6;      % spike threshold in standard deviations (4)		
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])		
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])		
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])		
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)		
ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)		
		
% load predefined principal components (visualization only (Phy): used for features)		
dd                  = load('PCspikes2.mat'); % you might want to recompute this from your own data		
ops.wPCA            = dd.Wi(:,1:7);   % PCs 		
		
% options for posthoc merges (under construction)		
ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)		
ops.epu     = Inf;		

if ispc
    dmem=memory;
    ops.ForceMaxRAMforDat   = dmem.MaxPossibleArrayBytes;
else
    ops.ForceMaxRAMforDat   = 20e9; % maxim
end

save(fops, '-v7.3', '-struct', 'ops')


end

function info = getHeader(hdr)
eval(char(hdr'));
info.header = header;
end

