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

%% Step 0: setup for import

% Make sure you have already copied the data from the server to the local 
% disk or the entire import will be unbearably slow

% NOTE: I have this set up for a session where the Shank2 probe is plugged
%       into headstage 1, and a single electode is plugged into headstage 2

% --- choose the session directory you want to analyze
oepath = uigetdir();

% --- setup the hardware you used
headstages = {hardware.headstage.intan_RHD2132, hardware.headstage.intan_RHD2132};

% list single electrode channels (this can be a vector if > 1 electrode used)
singleElectrodeChannelNumber = 4; % this is relative to the start of that headstage (not absolute channel #)

electrodes = {hardware.electrode.Shank2, hardware.electrode.SingleElectrodes(singleElectrodeChannelNumber)};

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
ops = io.oe2dat(oepath, 'headstages', headstages, 'electrodes', electrodes);

%%
% TODO: make Ephys session take in a path
EphysSession

%% 
ops = io.loadOps(oepath);

preprocess.runKiloSort(ops)
%%
preprocess.KiloAutomerge(ops)

%%

[session, ops, info] = io.loadSession(oepath);

[lfp, lfpTime, lfpInfo] = io.getLFP(ops);

PDS = io.getPds(session);
%%
[data, timestamps, elInfo] = io.getEdf(ops, PDS, 1);