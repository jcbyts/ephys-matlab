%% Import Script
% This script shows how to import a session from pldaps using the full
% pipeline.

% Run this first-> 
addEphysMatlab

% run import session to start the process. It will have you select a
% folder, and will import the raw ephys, apply the channel map, and filter
% the LFP. Finally, it will add the session to the meta table, which you
% can then update.
%io.importSession();

io.importSession([], 'overwrite', false);

%% load the session from the meta table

meta = io.getExperimentsAnd(); % get all experiments meta data

%thisSession = meta(146,:); % choose the exact session you want to look at from the meta table 
thisSession = meta(end,:); % if you just want the last one
disp(thisSession)

%% Import the stimuli

thisSession = io.importStimulusProtocols(thisSession);

%% run artifact rejection

preprocess.removeArtifactsManual(thisSession);

%% sort with kilosort

preprocess.runKiloSort(thisSession, 'merge', true)

%% Ready to run phy? Run this cell and then follow the prompt

commandwindow
fprintf('Open the START menu and select Anaconda prompt\n')

fprintf('Copy and paste in the following\n')

fprintf('activate phy\n')
fprintf('phy template-gui params.py\n')

fprintf('Then cd into the directory that represents your data\n')

%% 

meta = io.getExperimentsAnd(); % get all experiments meta data

%thisSession = meta(146,:); % choose the exact session you want to look at from the meta table 
thisSession = meta(end,:); % if you just want the last one
disp(thisSession)

sp = io.getSpikesFromKilo(thisSession);

%% sort with threshold
preprocess.runSingleChannelSpikeSortThreshold(thisSession)

%% sort with mixture of gaussians
preprocess.runSingleChannelSpikeSortMog(thisSession)
