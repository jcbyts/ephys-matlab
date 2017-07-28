function EphysSessions

[dataPath, sessionList] = io.selectSessions();

path_to=fileparts(ops.fbinary);

if isempty(path_to)
    fname       = fullfile(ops.root, sprintf('%s.dat', ops.fbinary));
else
    fname       = ops.fbinary;
end

% --- build output directory and files
nSessions = numel(sessionList);

sessionInfo = repmat(struct('Subject', [], 'DateNum', [], 'Date', [], 'Time', [], 'Note', []), nSessions, 1);

% extract session info from directory naming convention
pat = '(?<subject>\D+)\_(?<date>[\d-]+)\_(?<time>[\d-]+)\_(?<note>[^]+)';
for i = 1:nSessions
    s = regexp(sessionList{i}, pat, 'names');
   sessionInfo(i).Subject = s.subject;
   sessionInfo(i).DateNum = datenum([s.date '-' s.time], 'yyyy-mm-dd-HH-MM-SS');
   sessionInfo(i).Date    = datestr(sessionInfo(i).DateNum, 'yyyy.mm.dd');
   sessionInfo(i).Time    = datestr(sessionInfo(i).DateNum, 'HH:MM:SS');
   sessionInfo(i).Note    = s.note;
end

% --- Check that all sessions come from the same subject
assert(all(arrayfun(@(x) strcmp(x.Subject, sessionInfo(1).Subject), sessionInfo)), 'These sessions do not appear to be from the same subject!')

[~, ind] = sort(arrayfun(@(x) x.DateNum, sessionInfo));

sessionInfo = sessionInfo(ind);

fprintf('Found [%d] sessions from subject [%s]\n', nSessions, sessionInfo(1).Subject)
fprintf('Num\tDate\t\tTime\t\tNote\n')
for i = 1:nSessions
   fprintf('%d)\t%s\t%s\t%s\n', i, sessionInfo(i).Date, sessionInfo(i).Time, sessionInfo(i).Note)
end

%% TO DO: build new directory and concatenate binary file for analysis

fidout      = fopen(fname, 'w');
%
clear fs

for j = 1:ops.Nchan
    tmp = dir(fullfile(ops.root, sprintf('*CH%d.continuous',j) ));
    tmp_ = dir(fullfile(ops.root, sprintf('*CH%d_*.continuous', j) )); % if separate files are saved by open ephys gui
    
    fs{j} = [tmp tmp_];
end


nblocks = cellfun(@(x) numel(x), fs);
if numel(unique(nblocks))>1
   error('different number of blocks for different channels!') 
end
%
nBlocks     = unique(nblocks);
nSamples    = 1024;  % fixed to 1024 for now!

fid = cell(ops.Nchan, 1);

tic
for k = 1:nBlocks
    for j = 1:ops.Nchan
        fid{j}             = fopen(fullfile(ops.root, fs{j}(k).name));
        % discard header information
        fseek(fid{j}, 1024, 0);
    end
    %
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

toc