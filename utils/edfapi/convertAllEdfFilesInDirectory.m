
fl = dir('*edf');

for i = 1:numel(fl)
    fname = fl(i).name;
    edf = edfmex(fname);
    save(strrep(fname, '.edf', '.mat'), '-struct', 'edf');
end