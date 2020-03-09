function varargout = concatenateSessions(dataPath, varargin)
% CONCATENATE SESSIONS concatenates multiples sessions together

ip = inputParser();
ip.addParameter('ChannelIndex', 1:32)
ip.addParameter('CsdAlign', true)
ip.addParameter('CsdMaxShift', 4)
ip.addParameter('CsdAnchorSession', [])

ip.parse(varargin{:})

if nargin<1
    dataPath = uigetdir(pwd, 'Select Data Directory');
end

[dataPath, sessionList] = io.selectSessions(dataPath);

% --- declare constants
NUM_HEADER_BYTES    = 1024; % size of header for each analog file
SAMPLES_PER_RECORD  = 1024;
BLOCK_BYTES         = 2070; % see load_open_ephys_data_faster

% --- build output directory and files
nSessions = numel(sessionList);

sessionInfo = repmat(struct('subject', [], 'dateNum', [], 'date', [], 'time', [], 'note', []), nSessions, 1);

% extract session info from directory naming convention
pat = '(?<subject>\D+)\_(?<date>[\d-]+)\_(?<time>[\d-]+)\_(?<note>[^]+)';
for i = 1:nSessions
   s = regexp(sessionList{i}, pat, 'names');
   sessionInfo(i).subject = s.subject;
   sessionInfo(i).dateNum = datenum([s.date '-' s.time], 'yyyy-mm-dd-HH-MM-SS');
   sessionInfo(i).date    = datestr(sessionInfo(i).dateNum, 'yyyy-mm-dd');
   sessionInfo(i).time    = datestr(sessionInfo(i).dateNum, 'HH:MM:SS');
   sessionInfo(i).note    = s.note;
   sessionInfo(i).path    = fullfile(dataPath, sessionList{i});
   
   % --- get channel info
   fl  = dir(fullfile(sessionInfo(i).path, '*CH*.continuous'));
   str = 'CH(\d+)';  
   chlist = cellfun(@(x) cell2mat(regexp(x, str, 'match')) , {fl.name}, 'UniformOutput', false);
   
   sessionInfo(i).channels = unique(cellfun(@(x) str2double(x(3:end)), chlist));
   sessionInfo(i).numChannels = numel(sessionInfo(i).channels);
   sessionInfo(i).ops         = dir(fullfile(sessionInfo(i).path, 'ops.mat'));
end

% --- Check that all sessions come from the same subject
assert(all(arrayfun(@(x) strcmp(x.subject, sessionInfo(1).subject), sessionInfo)), 'These sessions do not appear to be from the same subject!')

% --- Check that Ephys Session has been run on one of the sessions already


% sort in chronological order
[~, ind] = sort(arrayfun(@(x) x.dateNum, sessionInfo));

sessionInfo = sessionInfo(ind);

% --- print out session info to the command window
fprintf('Found [%d] sessions from subject [%s]\n', nSessions, sessionInfo(1).subject)
fprintf('Num\tDate\t\tTime\t\tChannels\tNote\n')
for i = 1:nSessions
   fprintf('%d)\t%s\t%s\t%d\t\t\t%s\n', i, sessionInfo(i).date, sessionInfo(i).time, sessionInfo(i).numChannels, sessionInfo(i).note)
end

% --- get index for channels that are common to all sessions
ChannelIndex = ip.Results.ChannelIndex;

goodChannels = sessionInfo(1).channels;
for i = 1:nSessions
    goodChannels = intersect(goodChannels, sessionInfo(i).channels);
end
goodChannels = intersect(goodChannels, ChannelIndex);

% --- Make new directory and concatenate binary file for analysis
fpath = fullfile(dataPath, sprintf('%s_%d-sessions_%s_%s_Ch%d-%d', sessionInfo(1).subject, nSessions, sessionInfo(1).date, sessionInfo(end).date, goodChannels(1), goodChannels(end)));

if ~exist(fpath, 'dir')
    mkdir(fpath)
end

save(fullfile(fpath, 'session_info.mat'), 'sessionInfo')

% --- check that at least one of the sessions has an ops
hasOps = arrayfun(@(x) ~isempty(x.ops), sessionInfo);
assert(any(hasOps), 'You must run EphysSession on at least one session before concatenating.')

ops = load(fullfile(sessionInfo(find(hasOps,1)).path, 'ops.mat'));
if isfield(ops, 'parfor_')
    ops.parfor = ops.parfor_;
    ops = rmfield(ops, 'parfor_');
end

ops.root  = fpath;
ops.Nchan = numel(goodChannels);

% -------------------------------------------------------------------------
% --- Try to allign channels using the CSD
if ip.Results.CsdAlign
    fprintf('Attempting allignment based on the CSD\n')
    
    % --- Step 1: Calculate CSD for each session
    clear CSD
    sessionIx = 1:numel(sessionInfo);
    figure; clf
    for i = 1:numel(sessionInfo)
        
        
        oepath = sessionInfo(i).path;
        CSD(i) = session.computeCSD(oepath);
        [sess, ~, ~] = io.loadSession(oepath);
        subplot(1,numel(sessionIx), i)
        imagesc(CSD(i).time, CSD(i).depths, CSD(i).csd); hold on
        plot(CSD(i).time(1:CSD(i).imRescale:end), bsxfun(@plus, CSD(i).pots, CSD(i).depths(1:CSD(i).imRescale:end)), 'k')
        colormap jet
        title(sess.date)
        
        drawnow
    end
    
    % --- Step 2: compare to anchor point
    if isempty(ip.Results.CsdAnchorSession)
        baseCsd = ceil(numel(sessionInfo)/2);
    else
        baseCsd = ip.Results.CsdAnchorSession;
    end
    
    shiftSize = nan(numel(CSD),1);

    for kCsd = 1:numel(CSD)

        sz = size(CSD(baseCsd).csd);
        maxShift = ip.Results.CsdMaxShift;
        
        % shift by integer channel numbers
        nChannels = ceil(sz(1)/CSD(baseCsd).imRescale);
        
        baseIdx = (maxShift+1):(nChannels-maxShift);
        
        shifts = -maxShift:maxShift;
        nShifts = numel(shifts);
        
        csd1 = CSD(baseCsd).csd(1:CSD(baseCsd).imRescale:end,:);
        csd2 = CSD(kCsd).csd(1:CSD(kCsd).imRescale:end,:);
        
        figure(1); clf
        
        for iShift = 1:nShifts
            
            
            shift = shifts(iShift);
            testIdx = baseIdx + shift;
            
            
            subplot(1,3,1)
            imagesc(csd1(baseIdx,:));
            subplot(1,3,2)
            imagesc(csd2(testIdx,:));
            subplot(1,3,3)
            
            a(iShift) = sum(sum(csd1(baseIdx,:) .* csd2(testIdx,:)));
            plot(shifts(1:iShift), a(1:iShift), '-o'); hold on
            pause(0.15)
            
        end
        [~, shiftId] = max(a);
        shiftSize(kCsd) = shifts(shiftId);
    end
    
    load(ops.chanMap)
    
    if ip.Results.CsdAlign
        goodChannels = goodChannels(ip.Results.CsdMaxShift:(numel(goodChannels)-ip.Results.CsdMaxShift));
    end
    
    ops.Nchan = numel(goodChannels);
    
    % --- get list of all files
    fs = cell(ops.Nchan,1);
    for j = 1:ops.Nchan
        
        
        fs{j} = [];
        for i = 1:nSessions
            
            % channel shifted by CSD shift size for this session
            ch = chanMap(goodChannels(j)+shiftSize(i));
            
            tmp = dir(fullfile(sessionInfo(i).path, sprintf('*CH%d.continuous',ch) ));
            % if separate files are saved by open ephys gui
            tmp_ = dir(fullfile(sessionInfo(i).path, sprintf('*CH%d_*.continuous', ch) ));
            
            fs{j} = [fs{j} tmp tmp_(:)'];
        end
    end
    
    ops.chanMap = 1:ops.Nchan;
    
    
    
else
    
    % --- get list of all files
    fs = cell(ops.Nchan,1);
    for j = 1:ops.Nchan
        
        ch = goodChannels(j);
        
        fs{j} = [];
        for i = 1:nSessions
            
            tmp = dir(fullfile(sessionInfo(i).path, sprintf('*CH%d.continuous',ch) ));
            % if separate files are saved by open ephys gui
            tmp_ = dir(fullfile(sessionInfo(i).path, sprintf('*CH%d_*.continuous', ch) ));
            
            fs{j} = [fs{j} tmp tmp_(:)'];
        end
    end
    
    
end



% --- concatenate binary files
ops.fbinary = fullfile(ops.root, 'ephys.dat');
path_to=fileparts(ops.fbinary);

if isempty(path_to)
    fname       = fullfile(ops.root, sprintf('%s.dat', ops.fbinary));
else
    fname       = ops.fbinary;
end

if ~exist(fname, 'file')
    
    fprintf('Creating concatenated binary file for all sessions\n')
    
    fidout      = fopen(fname, 'w');
    
    nblocks = cellfun(@(x) numel(x), fs);
    if numel(unique(nblocks))>1
        error('different number of blocks for different channels!')
    end
    
    nBlocks     = unique(nblocks);
    nSamples    = SAMPLES_PER_RECORD;  % fixed to 1024 for now!
    
    fid = cell(ops.Nchan, 1);
    
    fprintf('Found [%d] total blocks\n', nBlocks)
    
    % --- initiali
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
        
        for j = 1:ops.Nchan
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
                assert(nFragments==numel(ts_), 'TIMESTAMPS and FRAGMENTS are different sizes')
                
                % get datetime of recording start
                dn_ =datenum(info_.header.date_created, 'dd-mmm-yyyy HHMMSS');
                
                info.dateNum = [info.dateNum repmat(dn_,1,nFragments)];
                info.timestamps = [info.timestamps ts_];
                info.fragments  = [info.fragments ns_];
                
                % --- skip header
                fseek(fid{j}, NUM_HEADER_BYTES, 'bof');
            end
            
            
        end
        
        
        % --- save raw data to new binary file
        nsamps = 0;
        flag = 1;
        while 1
            samples = zeros(nSamples * 1000, ops.Nchan, 'int16');
            for j = 1:ops.Nchan
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
        ops.nSamplesBlocks(k) = nsamps;
        
        for j = 1:ops.Nchan
            fclose(fid{j});
        end
        
    end
    
    fclose(fidout);
    
    save(fullfile(fpath, 'ephys_info.mat'), '-v7.3', '-struct', 'info')
    
    fprintf('Done [%02.2fs]\n', toc)
    
end

ops.bitVolts = info.bitVolts;
ops.fs =info.sampleRate;

if ip.Results.CsdAlign
   
    chanMap   = 1:ops.Nchan;
    n         = numel(chanMap);
    connected = true(n, 1); % connected(1:2) = 0;
    xcoords   = zeros(1,n);
    ycoords   = 50*(1:n);
    
    % Often, multi-shank probes or tetrodes will be organized into groups of
    % channels that cannot possibly share spikes with the rest of the probe. This helps
    % the algorithm discard noisy templates shared across groups. In
    % this case, we set kcoords to indicate which group the channel belongs to.
    % In our case all channels are on the same shank in a single group so we
    % assign them all to group 1.
    
    kcoords = ones(1,n);
    fs = ops.fs;
    ops.chanMap = fullfile(ops.root, 'chanMap.mat');
    save(ops.chanMap, 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs')
    

 
    
    
end

save(fullfile(fpath, 'ops.mat'), '-v7.3', '-struct', 'ops')

if nargout > 0
    varargout{1} = fpath;
end


end

function info = getHeader(hdr)
eval(char(hdr'));
info.header = header;
end

