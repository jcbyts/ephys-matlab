
server_directory = 'Z:\Data\PLDAPS\Ellie\';
local_directory  = 'C:\Data';


fldlist = dir(local_directory);
fldlist(1:2) = []; % remove . and ..
% only keep directories
fldlist = fldlist(arrayfun(@(x) x.isdir, fldlist));

% loop over directories
nSessions = numel(fldlist);
for kDir = 1:nSessions
    
    all_data = dir(fullfile(local_directory, fldlist(kDir).name));
    all_data(1:2) = [];
    derived_data = dir(fullfile(local_directory, fldlist(kDir).name, '_*'));
    all_ = arrayfun(@(x) x.name, all_data, 'UniformOutput', false);
    derived_ = arrayfun(@(x) x.name, derived_data, 'UniformOutput', false);
    
    delete_list = setdiff(all_, derived_);
    
    % make sure you don't delete any mat or m files
    mfiles = dir(fullfile(local_directory, fldlist(kDir).name, '*.m*'));
    m_ = arrayfun(@(x) x.name, mfiles, 'UniformOutput', false);
    delete_list = setdiff(delete_list, m_);
    
    if isempty(derived_data)
        disp('no derived data in this folder. skipping')
        continue
    end
    
    for i = 1:numel(derived_data)
        local_  = fullfile(local_directory, fldlist(kDir).name, derived_{i});
        server_ = fullfile(server_directory, fldlist(kDir).name, derived_{i});
        
        if ~exist(server_, 'dir')
            fprintf('[%s] copying to [%s]\n', derived_{i}, fullfile(server_directory, fldlist(kDir).name));
            copyfile(local_, server_);
            fprintf('Done\n');
        end
        
    end
    
    for j = 1:numel(delete_list)
        
        local_  = fullfile(local_directory, fldlist(kDir).name, delete_list{j});
        server_ = fullfile(server_directory, fldlist(kDir).name, delete_list{j});
        
        if exist(server_, 'file')
            fprintf('deleting [%s]\n', local_)
            delete(local_)
        else
            copyfile(local_, server_)
        end
        
    end
    
    
%     if input('continue?')~=1
%         break
%     end
end


