
% --- List of paths that are required
reposDir = 'C:\Users\Jake\Repos';
dirs{1} = fullfile(reposDir, 'KiloSort');
dirs{2} = fullfile(reposDir, 'npy-matlab');
dirs{3} = fullfile(reposDir, 'spikes');
dirs{4} = fullfile(reposDir, 'sortingQuality');

reposDir = 'C:\Users\Jake\Dropbox\MatlabCode\Repos\';
dirs{5}  = fullfile(reposDir, 'scalablerf'); % for receptive field estimation
dirs{6}  = fullfile(reposDir, 'pdstools');   % useful matlab functions for ephys
dirs{7}  = fullfile(reposDir, 'ncclabcode'); % general useful matlab functions
dirs{8}  = fullfile(reposDir, 'PLDAPS');     % stimulus generation
dirs{9} = fullfile(reposDir, 'edfmex');      % mex file for reading edf
% dirs{10}  = fullfile(reposDir, 'PLDAPStools'); % tools for interacting with PDLAPS

% warning off % turn warnings off
for j=1:length(dirs)
    a=genpath(dirs{j});
    if any(strfind(computer, 'WIN'))
        b=textscan(a,'%s','delimiter',';');
    else
        b=textscan(a,'%s','delimiter',':');
    end
    b=b{1};
    b(~cellfun(@isempty,strfind(b,'.git')))=[];
    addpath(b{:})
    fprintf('[%s] added to the path\n', dirs{j})
end

% warning on