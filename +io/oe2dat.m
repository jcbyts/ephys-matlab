function [ops, info] = oe2dat(oepath, shanks, varargin)
% ops = oe2dat(oepath, shanks, varargin)
% Inputs:
%   oepath@string - full path to the data
%   shanks@cell   - array of electrodes used
%
% Optional Inputs (argument pairs):
%   overwrite@logical - overwrite imported data (default: false)
%   verbose@logical   - print to the command window
ops = [];

if nargin < 2
    shanks{1}.name = 'dummy';
end

numShanks  = numel(shanks);
shankPaths = cell(numShanks,1);
for i = 1:numShanks
    shankPaths{i} = sprintf('_shank%02.0f_%s', i, shanks{i}.name);
end

% --- parse optional arguments
ip = inputParser();
ip.addParameter('overwrite', false)
ip.addParameter('verbose', true)
ip.parse(varargin{:});

% -------------------------------------------------------------------------
% Setup for data conversion


% --- save session info
fsess   = fullfile(oepath, 'session_info.mat'); % meta data about the session

% --- build meta data about session (Subject, Data, Time, Note, Path)
sessionInfo = struct('subject', [], 'dateNum', [], 'date', [], 'time', [], 'note', [], 'nShanks', []);

% extract session info from directory naming convention
pat = '(?<subject>\D+)\_(?<date>[\d-]+)\_(?<time>[\d-]+)\_(?<note>[^]+)';

[~, sessionPath] = fileparts(oepath);
s = regexp(sessionPath, pat, 'names');
sessionInfo.subject = s.subject;
sessionInfo.dateNum = datenum([s.date '-' s.time], 'yyyy-mm-dd-HH-MM-SS');
sessionInfo.date    = datestr(sessionInfo.dateNum, 'yyyy-mm-dd');
sessionInfo.time    = datestr(sessionInfo.dateNum, 'HH:MM:SS');
sessionInfo.note    = s.note;
sessionInfo.path    = oepath;
sessionInfo.nShanks = numShanks;
for i = 1:numShanks
    sessionInfo.shankPaths{i} = fullfile(sessionInfo.path, shankPaths{i});
    if ~exist(sessionInfo.shankPaths{i}, 'dir')
        mkdir(sessionInfo.shankPaths{i})
    end
end

% overwrite session info??
save(fsess, 'sessionInfo')

maxChan = 0;
% --- loop over shanks and import
clear ops

for iShank = 1:numShanks
    
    % combine headstages if more than one
    if numel(shanks{iShank}.headstages) == 1
        headstage = shanks{iShank}.headstages{1};
    else
        headstage = combineHeadstages(shanks{iShank}.headstages);
    end
   
   % build channel map
   chanMap = headstage.channelMap(shanks{iShank}.channelMap) + maxChan;
   xcoords = shanks{iShank}.xcoords;
   ycoords = shanks{iShank}.ycoords;
   zcoords = shanks{iShank}.zcoords;
   
   [ops(iShank), info(iShank)] = oe2dat_helper(oepath, shankPaths{iShank}, chanMap, xcoords, ycoords, zcoords, ip.Results.verbose, ip.Results.overwrite);
   
   maxChan = maxChan + max(headstage.channelMap);
end



function [ops, info] = oe2dat_helper(oepath, shankName, chanMap, xcoords, ycoords, zcoords, verbose, overwrite)

if ~exist('verbose', 'var')
    verbose = false;
end

if ~exist('overwrite', 'var')
    overwrite = false;
end

% --- declare constants
NUM_HEADER_BYTES    = 1024; % size of header for each analog file
SAMPLES_PER_RECORD  = 1024;
BLOCK_BYTES         = 2070; % see load_open_ephys_data_faster for how this is calculated

% --- name of output files
fephys  = fullfile(oepath, shankName, 'ephys.dat');        % raw data (binary file)
fops    = fullfile(oepath, shankName, 'ops.mat');          % spike-sorting struct
finfo   = fullfile(oepath, shankName, 'ephys_info.mat');   % meta data about the binary

% --- Exit if this has already been run
if exist(fops, 'file') && ~overwrite
    ops  = load(fops);
    info = load(finfo);
    return
else
    disp('replacing old files');
end

numChannels = numel(chanMap);

fs = cell(numChannels,1);
for j = 1:numChannels
    ch = chanMap(j);
    
    if verbose
        fprintf('*CH%d.continuous\n',ch);
    end
    
    tmp = dir(fullfile(oepath, sprintf('*CH%d.continuous',ch) ));
    % if separate files are saved by open ephys gui
    tmp_ = dir(fullfile(oepath, sprintf('*CH%d_*.continuous', ch) ));
    
    fs{j} = [tmp tmp_(:)'];
end

% -------------------------------------------------------------------------
% --- convert raw data

% open binary file for writing
fidout      = fopen(fephys, 'w');

% number of recording blocks must be the same for all channels (this will
% never be violated unless one of the headstages was unplugged during
% recording)
nblocks = cellfun(@(x) numel(x), fs);
if numel(unique(nblocks))>1
    error('different number of blocks for different channels! Split the recordings into separate folders')
end

fid = cell(numChannels, 1);

nBlocks     = unique(nblocks);
nSamples    = SAMPLES_PER_RECORD;  % fixed to 1024 for now!
fprintf('Found [%d] total blocks\n', nBlocks)

% --- initialize meta data struct
info = struct('format', [], ... % format the original data was saved in
    'version', [], ... % version of the OE GUI
    'dateNum', [], ... % date/time of each recording
    'sampleRate', [], ... % sampling rate of data collection
    'bitVolts', [], ...   % multiplier to convert integers to millivolts
    'timestamps', [], ... % timestamp for each recording fragment
    'fragments', []);     % number of samples per recording fragment

tic
for k = 1:nBlocks
    fprintf('Block %d\n', k)
    
    % preallocate samples
    filename = fullfile(fs{j}(k).folder, fs{j}(k).name);
	data_seg = load_open_ephys_data_faster(filename, 'unscaledInt16');
    samples  = zeros(numChannels, numel(data_seg), 'int16');
        
    for j = 1:numChannels
        fprintf('Channel %d\n', j)
        
        filename = fullfile(fs{j}(k).folder, fs{j}(k).name);
        [data_seg, timestamps_seg, info_] = load_open_ephys_data_faster(filename, 'unscaledInt16');
        
        assert(info_.header.version == 0.4, 'open ephys version has changed. You must now check that everything still works. Sorry, Jake')
        
        % --- save header info
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
            % save timestamps of all starts and stops
            ts_ = timestamps_seg(1);
            ns_ = numel(data_seg);
            stops = find(round(diff(info_.header.sampleRate*timestamps_seg))~=1);
            
            if any(stops)
                ts_ = [ts_ timestamps_seg(stops+1)'];
                ns_ = diff([0 stops(:)' ns_]);
                
                % Uncomment this if debugging timestamps
                %                  fprintf('Found multiple paused segments\n')
                %                  figure
                %                  plot(seg, '.')
                %                  hold on
                %                  plot(cumsum([[0 ns_(1:end-1)]; ns_], 2)/SAMPLES_PER_RECORD, [ts_; ts_])
                
                
            end
            
            nFragments = numel(ns_);
            assert(nFragments==numel(ts_), 'oe2dat: TIMESTAMPS and FRAGMENTS are different sizes')
            
            % get datetime of recording start
            dn_ =datenum(info_.header.date_created, 'dd-mmm-yyyy HHMMSS');
            
            info.dateNum    = [info.dateNum repmat(dn_,1,nFragments)];
            info.timestamps = [info.timestamps ts_];
            info.fragments  = [info.fragments ns_];
        end
        
        samples(j,:) = data_seg;
    end
    
    nsamps = size(samples,2);
    fwrite(fidout, samples, 'int16');
    
    nSamplesBlocks(k) = nsamps; %#ok<AGROW>
    
end

fclose(fidout);

% convert sample timestamps to time
% info.timestamps = info.timestamps/info.sampleRate;

save(finfo, '-v7.3', '-struct', 'info')

fprintf('Done [%02.2fs]\n', toc)


% --- build ops
ops = struct();
ops.root    = fullfile(oepath, shankName);
ops.fbinary = fullfile(ops.root, 'ephys.dat');

ops.Nchan          = numChannels;
ops.nSamplesBlocks = nSamplesBlocks;
ops.bitVolts       = info.bitVolts;
ops.fs             = info.sampleRate;

% --- Parameters for KiloSort
ops.datatype            = 'dat';  % binary
ops.fproc               = fullfile(ops.root, 'tempWh.dat'); % residual from RAM of preprocessed data				
% define the channel map as a filename (string) or simply an array
ops.chanMap             = fullfile(ops.root, 'chanMap.mat'); % make this file using createChannelMapFile.m

chanMap   = 1:numChannels;
connected = true(size(xcoords));
fs        = ops.fs;

save(ops.chanMap, 'chanMap', 'xcoords', 'ycoords', 'zcoords', 'connected', 'fs');
% if ~exist(ops.chanMap, 'file')
% 	ops.chanMap = 1:ops.Nchan; % treated as linear probe if unavailable chanMap file
% end

% default number of clusters
ops.Nfilt               = 2*(ops.Nchan+32-mod(ops.Nchan,32));  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)     		
ops.nNeighPC            = 12; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)		
ops.nNeigh              = 12; % visualization only (Phy): number of neighboring templates to retain projections of (16)		
		
% options for channel whitening		
ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)		
ops.nSkipCov            = 2; % compute whitening matrix from every N-th batch (1)		
ops.whiteningRange      = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)		
		
ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info). 		

% other options for controlling the model and optimization		
ops.Nrank               = 3;        % matrix rank of spike template model (3)		
ops.nfullpasses         = 6;        % number of complete passes through data during optimization (6)		
ops.maxFR               = 500e3;    % maximum number of spikes to extract per batch (20000)		
ops.NotchFilter60       = false;
ops.softwareReferencing = 'none';
ops.fshigh              = 300;   % frequency for high pass filtering		
% ops.fslow               = 2500;   % frequency for low pass filtering (optional)
ops.artifactThresh      = 500;
ops.artifactNchans      = 12;

ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection		
ops.scaleproc           = 200;   % int16 scaling of whitened data		
ops.NT                  = 2*128*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory) 		
% for GPU should be multiple of 32 + ntbuff		
		
% the following options can improve/deteriorate results. 		
% when multiple values are provided for an option, the first two are beginning and ending anneal values, 		
% the third is the value used in the final pass. 		
ops.Th               = [4 10 10];    % threshold for detecting spikes on template-filtered data ([6 12 12])		
ops.lam              = [5 5 5];   % large means amplitudes are forced around the mean ([10 30 30])		
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)		
ops.momentum         = 1./[20 1000];  % start with high momentum and anneal (1./[20 1000])		
ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)		
ops.mergeT           = .1;           % upper threshold for merging (.1)		
ops.splitT           = .1;           % lower threshold for splitting (.1)		
		
% options for initializing spikes from data		
ops.initialize      = 'fromData';    %'fromData' or 'no'		
ops.spkTh           = -4;      % spike threshold in standard deviations (4)		
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])		
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])		
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])		
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)		
ops.nFiltMax        = 50e3;   % maximum "unique" spikes to consider (10000)		
		
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

function info = getHeader(hdr)
eval(char(hdr'));
info.header = header;

function [headstage] = combineHeadstages(headstages)

% --- build channel map info
numHeadstages = numel(headstages);

nHeadstageChannels = zeros(numHeadstages,1);
for iHeadstage = 1:numHeadstages
    nHeadstageChannels(iHeadstage) = sum(~isnan(headstages{iHeadstage}.channelMap));
end

% build new combined headstage
headstage = hardware.headstage.headstage;

chOffset = zeros(numHeadstages,1);
chCount  = zeros(numHeadstages,1);
for iHeadstage = 1:numHeadstages
    headstage.name         = [headstage.name headstages{iHeadstage}.name ','];
    headstage.manufacturer = [headstage.manufacturer headstages{iHeadstage}.manufacturer ','];
    headstage.model        = [headstage.model headstages{iHeadstage}.model ','];
    headstage.connector    = [headstage.connector headstages{iHeadstage}.connector ','];
    headstage.filter       = [headstage.filter; headstages{iHeadstage}.filter];
    headstage.samplingRate = [headstage.samplingRate headstages{iHeadstage}.samplingRate];
    headstage.gains        = [headstage.gains headstages{iHeadstage}.gains];
    tmp = max(headstage.channelMap);
    if isempty(tmp)
        chOffset(iHeadstage) = 0;
    else
        chOffset(iHeadstage) = tmp;
    end
    headstage.channelMap = [headstage.channelMap headstages{iHeadstage}.channelMap + chOffset(iHeadstage)];
    chCount(iHeadstage)  = numel(headstages{iHeadstage}.channelMap);
end

% remove comma at end
headstage.name(end) = [];
headstage.manufacturer(end) = [];
headstage.model(end) = [];
headstage.connector(end) = [];
