function addEphysMatlab
pathto = fileparts(mfilename('fullpath'));

% --- List of paths that are required
% reposDir = 'C:\Users\Jake\Repos';
reposDir = getpref('ephysmatlab', 'repos');
kDir = 1;
dirs{kDir} = fullfile(reposDir, 'KiloSort');
kDir = kDir + 1;
dirs{kDir} = fullfile(reposDir, 'npy-matlab');
kDir = kDir + 1;
dirs{kDir} = fullfile(reposDir, 'spikes');
kDir = kDir + 1;
dirs{kDir} = fullfile(reposDir, 'sortingQuality');
kDir = kDir + 1;
dirs{kDir} = fullfile(reposDir, 'analysis-tools'); % open ephys tools



% kDir = kDir + 1;
% dirs{kDir}  = fullfile(reposDir, 'scalablerf'); % for receptive field estimation
% kDir = kDir + 1;
% dirs{kDir}  = fullfile(reposDir, 'pdstools');   % useful matlab functions for ephys
% kDir = kDir + 1;
% dirs{kDir}  = fullfile(reposDir, 'ncclabcode'); % general useful matlab functions
kDir = kDir + 1;
dirs{kDir}  = fullfile(reposDir, 'PLDAPS');     % stimulus generation
% kDir = kDir + 1;
% dirs{kDir} = fullfile(reposDir, 'edfmex');      % mex file for reading edf
kDir = kDir + 1;
dirs{kDir}  = fullfile(reposDir, 'pds-stimuli');
% dirs{10}  = fullfile(reposDir, 'PLDAPStools'); % tools for interacting with PDLAPS

% warning off % turn warnings off
for j=1:length(dirs)
    try
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
    catch
        warning('[%s] not added to the path\n', dirs{j})
    end
end

addpath(pathto)
addpath(fullfile(pathto, 'utils'))
addpath(fullfile(pathto, 'utils', 'edfapi'))
% warning on