function import_function(oepath, shank)
%% IMPORT SCRIPT
% This script shows you how to import an ephys session collected with
% PLDAPS / Open-Ephys / Eyelink
% In the future it will support other eyetrackers, but for the time being,
% it expects that every session will contain three types of files
% 1) Open-Ephys files (.continuous, .spikes, .events) - These are binary
%           files. They contain the raw electrophysiological data and meta 
%           data about the session.
% 2) PLDAPS files (*.PDS) - These are .mat files created by PLDAPS. They
%           contain the stimulus and behavioral data, as well as the meta
%           data required to synchronize the Ephys and Display computer
%           clocks
% 3) Eyelink files (*.edf) - These are proprietrary filetpe from SR
%           research for storying eyelink data.

% cd C:\Users\Jake\Repos\ephys-matlab\
% addEphysMatlab
%% Step 0: setup for import

% Make sure you have already copied the data from the server to the local 
% disk or the entire import will be unbearably slow

% NOTE: I have this set up for a session where the Shank2 probe is plugged
%       into headstage 1, and a single electode is plugged into headstage 2

% --- choose the session directory you want to analyze
% % oepath = uigetdir();
% 
% % --- setup the hardware you used
% clear shank
% 
% % setup shank2 (plugged into first headstage)
% shank{1} = hardware.electrode.Shank2;
% % specify what type of headstage was used
% shank{1}.headstages{1} = hardware.headstage.intan_RHD2132;
% shank{1}.name = 'V1';
% 
% % For the single electrodes specify a custom channel map using numbers
% % relative to the headstage start. eg., if using ch36 plugged into
% % headstage 2, this should be channel 4 (relative to the start of
% % headstage 2)
% 
% % % list single electrode channels (this can be a vector if > 1 electrode used)
% % chanMap = 4;
% % shank{2} = hardware.electrode.customChannelMap(chanMap);
% % shank{2}.name = 'MtBurrHoleMapping';
% % If you chose hardware.electrode.customChannelMap(chNum), the channel map
% % will be chanMap. That's it. No specifying headstage necessary.


%% Step 1: convert the raw ephys data
% we represent all of our ephys data as a binary file (integers only) and a
% *.mat file that specifies how to recover timestamps / voltage from the
% binary file. 

% --- The ops struct
% Most of the information about the session, probe, and parameters for 
% spike-sortinf using KiloSort are contained in the ops struct. It is
% required for KiloSort, but I added some extra fields that are useful for
% our book-keeping

% --- convert raw data (this can be slow)
io.oe2dat(oepath, shank);


%% Step 2: Spike Sort Array Shanks with Kilo
ops = io.loadOps(oepath);

% which "Shanks" correspond to arrays?
shanksToSortWithKilo = find(cellfun(@(x) numel(x.channelMap)>16, shank));

% Sorting channels that are not part of a single array will result in
% errors

for i = 1:numel(shanksToSortWithKilo)
    iShank = shanksToSortWithKilo(i);
    
    % run Kilosort
    preprocess.runKiloSort(ops(iShank), 'GPU', true, 'parfor', true, 'verbose', true, 'showfigures', false, 'merge', false);

    % run KiloAutomerge
    preprocess.KiloAutomerge(ops(iShank));
    
    
    fprintf('--------------------------------------------------------------\n')
    fprintf('--------------------------------------------------------------\n')
    fprintf('For manual stage using Phy run "Anaconda Prompt" from the START menu\n')
    fprintf('Enter the following lines into the prompt:\n\n')
    fprintf('activate phy\n')
    fprintf('cd %s\n', ops(iShank).root)
    fprintf('phy template-gui params.py\n')
    commandwindow
    pause(10)
    fprintf('\n\n')
    fprintf('For keyboard shortcuts and sorting instructions, go to:\n')
    fprintf('http://phy-contrib.readthedocs.io/en/latest/template-gui/#keyboard-shortcuts\n')
    
    fprintf('importing spikes to matlab\n')
    io.getSpikesFromKilo(ops(iShank));
end


%% Step 3: Spike sort single electrodes (if they exist)
singleTrodes = find(cellfun(@(x) numel(x.channelMap) < 4, shank));
if any(singleTrodes)
    clear sp
    for i = 1:numel(singleTrodes)
        iShank = singleTrodes(i);

        sp(i) = preprocess.runSingleChannelSpikeSortMog(ops(iShank));
    end
end

%% Step 4: Import LFP
% Loop over sessions and import the local field potential
% This is really slow because of the line-noise fitting and correction
ops = io.loadOps(oepath);

overwrite = false;
plotIt = false;
for i = 1:numel(ops)
    io.getLFP(ops(i),overwrite,plotIt);
end

%% Step 5: Import PLDAPS / eyelink files

% --- load session: This is central to all operations
[sess, ops, info] = io.loadSession(oepath);

overwrite = true; % if it breaks, run again with true 
PDS = io.getPds(sess, overwrite);


[data, timestamps, elInfo] = io.getEdf(sess, PDS, overwrite);

