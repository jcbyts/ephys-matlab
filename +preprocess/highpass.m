function datr = highpass(data, fs, fshigh, fslow, order)

if nargin < 5
    order = 3;
    if nargin < 4
        fslow = 0;
    if nargin < 3
        fshigh = 300;
        if nargin < 2
            error('preprocess.highpass: need to input sampling rate')
        end
    end
    end
end

if fslow<fs/2
    [b1, a1] = butter(order, [fshigh/fs,fslow/fs]*2, 'bandpass');
else
    [b1, a1] = butter(order, fshigh/fs*2, 'high');
end


sz = size(data);
if sz(2) > sz(1)
	trpose = true;
    datr = data';
else
    trpose = false;
    datr = data;
end
% filter
datr = filter(b1, a1, datr);
datr = flipud(datr);
datr = filter(b1, a1, datr);
datr = flipud(datr);

if trpose
    datr = datr';
end