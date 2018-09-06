function save_filtered_data(ops)
% fproc

if exist(ops.fproc, 'file')==2
    disp('file already exists')
    return
end


% build filter
if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
    [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
else
    [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
end

% load Raw
dataRAW = io.loadRaw(ops, [], false);
dataRAW = double(dataRAW');

% common average reference
dataRAW = bsxfun(@minus, dataRAW, mean(dataRAW,2));
    
% filter
datr = filter(b1, a1, dataRAW);
datr = flipud(datr);
datr = filter(b1, a1, datr);
datr = flipud(datr);

% convert to a 16-bit integer
datr = int16(datr');

% open file for writing
fid = fopen(ops.fproc, 'w');

% write to disk
fwrite(fid, datr, 'int16');
    
% cleaup
fclose(fid);
fprintf('Done\n')