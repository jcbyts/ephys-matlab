%% Import Script
% This script shows how to import a session from pldaps using the full
% pipeline.

% run import session to start the process. It will have you select a
% folder, and will import the raw ephys, apply the channel map, and filter
% the LFP. Finally, it will add the session to the meta table, which you
% can then update.
io.importSession();

%% load the session from the meta table

meta = io.getExperimentsAnd(); % get all experiments meta data

thisSession = meta(end,:);
disp(thisSession)

%% Import the stimuli

thisSession = io.importStimulusProtocols(thisSession);

%%

preprocess.runSingleChannelSpikeSortThreshold(thisSession)


%%

meta = io.getExperimentsAnd(); % get all experiments meta data

thisSession = meta(end,:);
disp(thisSession)

PDS = io.getPds(thisSession);

sp = io.getSpikes(thisSession, thisSession.SpikeSorting{1});

%%

% sp{1}

spatialMap = session.squareFlash(PDS);

%%
%x,y,x,y
window = [0 -5 5 0]; %[-3, 3]% degrees (window is the same in x and y -- TODO: separate)
binSize = 1; % degrees

stim = spatialMap.binSpace('window', window, 'binSize', binSize, 'correctEyePos', false);

nTimeBins = 20; % frames
Xd = spatialMap.buildDesignMatrix(stim, nTimeBins);


%% Plot spatial map for each unit
spikeTimes = sp{1}.st(sp{1}.clu==1);

sta = spatialMap.spikeTriggeredAverage(stim, Xd, spikeTimes, 10e3);


figure(1); clf
subplot(1,3,1)
imagesc(sta.fullRF)
subplot(1,3,2)
plot(sta.RFtime)
subplot(1,3,3)
imagesc(sta.RF)
%% play movie

for i = 1:numel(sta.RFtime)
    figure(2); clf
    subplot(1,2,1)
    imagesc(sta.fullRF); hold on
    plot(xlim, i*[1 1], 'r')
    subplot(1,2,2)
I = reshape(sta.fullRF(i,:), stim.size);
% plot(sta.fullRF(i,:))
imagesc(stim.xax, stim.yax, I, [min(sta.fullRF(:)) max(sta.fullRF(:))])
axis xy
drawnow
pause(0.5)
end

%%

y = histc(spikeTimes, stim.bins);
y(diff(stim.bins) > spatialMap.display.ifi*2) = 0;
y = sqrt(y);
nlags = 20;
xc = zeros(nlags*2+1,size(stim.X,2)); 
for i = 1:size(stim.X,2)
xc(:,i) = xcorr(y,stim.X(:,i), nlags);
end


%%
S = prf.fit_generic_linear_model(stim.X, y, 20, stim.size);


%%

imagesc(reshape(Xd(:,1:end-1)'*y, [], prod(stim.size)))

% clf
% imagesc(reshape(sum(reshape(Xd(:,1:end-1)'*y, [], prod(stim.size))), stim.size))

%%