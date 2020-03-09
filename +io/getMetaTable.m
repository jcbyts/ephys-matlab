function [data, meta_file] = getMetaTable()
% % Get Meta Table reads in the experimental meta data and returns a table
% meta = getMetaTable()

SERVER_DATA_DIR = getpref('EPHYS', 'SERVER_DATA');

meta_file = fullfile(SERVER_DATA_DIR, 'meta_data.xls');

[~, ~, alldata] = xlsread(meta_file);
    
data = cell2table(alldata(2:end,:), 'VariableNames', alldata(1,:));