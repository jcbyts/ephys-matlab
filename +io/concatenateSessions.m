function varargout = concatenateSessions(dataPath)
% CONCATENATE SESSIONS concatenates multiples sessions together

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
goodChannels = sessionInfo(1).channels;
for i = 1:nSessions
    goodChannels = intersect(goodChannels, sessionInfo(i).channels);
end

% --- Make new directory and concatenate binary file for analysis
fpath = fullfile(dataPath, sprintf('%s_%d-sessions_%s_%s', sessionInfo(1).subject, nSessions, sessionInfo(1).date, sessionInfo(end).date));

save(fullfile(fpath, 'session_info.mat'), 'sessionInfo')

if ~exist(fpath, 'dir')
    mkdir(fpath)
end

% check that at least one of the sessions has an ops
hasOps = arrayfun(@(x) ~isempty(x.ops), sessionInfo);
assert(any(hasOps), 'You must run EphysSession on at least one session before concatenating.')

ops = load(fullfile(sessionInfo(find(hasOps,1)).path, 'ops.mat'));
if isfield(ops, 'parfor_')
    ops.parfor = ops.parfor_;
    ops = rmfield(ops, 'parfor_');
end

ops.root  = fpath;
ops.Nchan = numel(goodChannels);

% --- get list of all files
fs = cell(ops.Nchan,1);
for j = 1:ops.Nchan
    
    ch = goodChannels(j);
    
    fs{j} = [];
    for i = 1:nSessions
        
        tmp = dir(fullfile(sessionInfo(i).path, sprintf('*CH%d.continuous',ch) ));
        % if separate files are saved by open ephys gui 
        tmp_ = dir(fullfile(sessionInfo(i).path, sprintf('*CH%d_*.continuous', ch) ));
        
        fs{j} = [fs{j} tmp tmp_];
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

save(fullfile(fpath, 'ops.mat'), '-v7.3', '-struct', 'ops')

if nargout > 0
    varargout{1} = fpath;
end


end

function info = getHeader(hdr)
eval(char(hdr'));
info.header = header;
end

