
SERVER_DATA_DIR = getpref('EPHYS', 'SERVER_DATA');
LOCAL_DATA_DIR  = getpref('EPHYS', 'LOCAL_DATA');


fldlist = dir(LOCAL_DATA_DIR);
fldlist(1:2) = []; % remove . and ..
% only keep directories
fldlist = fldlist(arrayfun(@(x) x.isdir, fldlist));

% loop over directories
nSessions = numel(fldlist);
for kDir = 1:nSessions
    
    all_data = dir(fullfile(LOCAL_DATA_DIR, fldlist(kDir).name));
    all_data(1:2) = []; % remove ., ..
    derived_data = dir(fullfile(LOCAL_DATA_DIR, fldlist(kDir).name, '_*'));
    all_ = arrayfun(@(x) x.name, all_data, 'UniformOutput', false);
    derived_ = arrayfun(@(x) x.name, derived_data, 'UniformOutput', false);
    
    delete_list = setdiff(all_, derived_);
    
    % make sure you don't delete any mat or m files
    mfiles = dir(fullfile(LOCAL_DATA_DIR, fldlist(kDir).name, '*.m*'));
    m_ = arrayfun(@(x) x.name, mfiles, 'UniformOutput', false);
    delete_list = setdiff(delete_list, m_);
    
    if isempty(derived_data)
        disp('no derived data in this folder. skipping')
        continue
    end
    
    for i = 1:numel(derived_data)
        local_  = fullfile(LOCAL_DATA_DIR, fldlist(kDir).name, derived_{i});
        server_ = fullfile(SERVER_DATA_DIR, fldlist(kDir).name, derived_{i});
        
        if ~exist(server_, 'dir')
            fprintf('[%s] copying to [%s]\n', derived_{i}, fullfile(SERVER_DATA_DIR, fldlist(kDir).name));
            copyfile(local_, server_);
            fprintf('Done\n');
        end
        
    end
    
    for j = 1:numel(delete_list)
        
        local_  = fullfile(LOCAL_DATA_DIR, fldlist(kDir).name, delete_list{j});
        server_ = fullfile(SERVER_DATA_DIR, fldlist(kDir).name, delete_list{j});
        
        if exist(server_, 'file')
            fprintf('deleting [%s]\n', local_)
            delete(local_)
        else
            status = copyfile(local_, server_);
            if status
                delete(local_)
            end
        end
        
    end
    
    
%     if input('continue?')~=1
%         break
%     end
end


