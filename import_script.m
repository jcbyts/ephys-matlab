%% Import Script
% This script shows how to import a session from pldaps using the full
% pipeline.

% run import session to start the process. It will have you select a
% folder, and will import the raw ephys, apply the channel map, and filter
% the LFP. Finally, it will add the session to the meta table, which you
% can then update.
thisSession = io.importSession();

%% Next steps:

% 1. spike sorting
% 2. stimulus import


