function [plot_buffer, ops] = preprocessVops(plot_buffer, ops)
% preprocessVops applies the filtering / referencing / channel map
% specified in ops

if exist(fullfile(ops.root, ops.chanMap), 'file')
    error('io.preprocessVops: Implement loading of channel map here')   
end

assert(size(plot_buffer,2)>size(plot_buffer,1), 'data is probably transposed')

% apply channel map
if numel(ops.chanMap) < ops.Nchan
    newChans = setdiff(1:ops.Nchan, ops.chanMap);
    ops.chanMap = [ops.chanMap(:); newChans(:)];
end
plot_buffer     = plot_buffer(ops.chanMap, :);

% notch filter
if ops.NotchFilter60
    wo = 60/(ops.fs/2);  bw = wo/35;
    [b,a] = iirnotch(wo,bw);
    
%     fs = ops.fs; fo = 60;  q = 35; bw = (fo/(fs/2))/q;
%     [b,a] = iircomb(fs/fo, bw, 'notch');

    datr = plot_buffer';

    datr = filter(b, a, datr);
    datr = flipud(datr);
    datr = filter(b, a, datr);
    datr = flipud(datr);

    plot_buffer = datr';
    
    
    wo = 180/(ops.fs/2);  bw = wo/35;
    [b,a] = iirnotch(wo,bw);
    
    datr = plot_buffer';

    datr = filter(b, a, datr);
    datr = flipud(datr);
    datr = filter(b, a, datr);
    datr = flipud(datr);

    plot_buffer = datr';
end

switch ops.softwareReferencing
    case 'commonaverage'
        plot_buffer     = bsxfun(@minus, plot_buffer, mean(plot_buffer));
    otherwise
end

% --- bandpass filtering
if isfield(ops, 'fshigh')
    
    % setup highpass or bandpass filter
    if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
        [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
    else
        [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
    end
    
    % filter
    datr = plot_buffer';
    
    datr = filter(b1, a1, datr);
	datr = flipud(datr);
	datr = filter(b1, a1, datr);
	datr = flipud(datr);
    
    if ~isinf(ops.artifactThresh)
        sizeThresh = ops.artifactThresh;
        chanThresh = ops.artifactNchans;
        
        bad=unique(bsxfun(@plus, find(sum(abs(datr)>sizeThresh,2)>chanThresh), -1:20));
        bad(bad<1)=[];
        bad(bad>size(datr,1))=[];
        
        datr(bad,:)=0;
    end

    plot_buffer = datr';
end

% --- convert to milivolts
for i = 1:ops.Nchan
    plot_buffer(i,:) = ops.bitVolts*plot_buffer(i,:);
end