% get sample dataset for richard

%% Find sessions that have movies 
meta = io.getExperimentsAnd('StimulusProtocols', 'HDmovies', 'SpikeSorting', 'Kilo');

%% load data
nSessions = size(meta, 1);

iSession = 1;

thisSession = meta(iSession,:);

PDS = io.getPds(thisSession);
spikes = io.getSpikes(thisSession);

trial = pds.getPdsTrialData(PDS);
disp('Done')

RF = rfmap.CoarseRFMap(PDS, spikes{1});

trialIndex = find(~arrayfun(@(x) isempty(x.HDmovies), trial));

goodTrials = trialIndex(arrayfun(@(x) x.HDmovies.use, trial(trialIndex)));

% movie shown on each trial
movieFilesTrial = arrayfun(@(x) x.HDmovies.moviefilename, trial(goodTrials), 'uni', 0);
movieList = unique(movieFilesTrial);

iMovie = 3;

thisMovieTrials = goodTrials(strcmp(movieFilesTrial, movieList{iMovie}));
numTrials = numel(thisMovieTrials);



PTB2OE = PDS{1}.PTB2OE;

movieFileName = movieList{iMovie};
frameTimesEphysClock = [];
timeShownMovieFile = [];

for iTrial = 1:numTrials

    thisTrial = thisMovieTrials(iTrial);

    ft = PTB2OE(trial(thisTrial).timing.flipTimes(3,:));
    ts = trial(thisTrial).HDmovies.frameShown;
    ts(isnan(ts))=[];
    assert(numel(ft)==numel(ts))
    frameTimesEphysClock = [frameTimesEphysClock ft(:)'];
    timeShownMovieFile = [timeShownMovieFile ts(:)'];
    
end

% Check RF mask
n = numel(RF);
sx = ceil(sqrt(n));
sy = round(sqrt(n));

xs = nan(n,1);
ys = nan(n,1);
sigmas = nan(n,1);
pxsiz = median(diff(RF(1).xax));
figure(1); clf
for i = 1:n
    [x,y,sigma] = radialcenter(RF(i).RF);
    subplot(sx,sy,i)
    imagesc(RF(i).RF)
    hold on
    plot(x, y, '.r')
    plot(x+[-sigma -sigma sigma sigma -sigma], y+[-sigma sigma sigma -sigma -sigma], 'r')
    set(gca, 'XTick', 1:numel(RF(i).xax), 'XTickLabel', RF(i).xax)
    drawnow
    xs(i) = interp1(RF(i).xax, x);
    ys(i) = interp1(RF(i).yax, y);
    sigmas(i) = sigma * pxsiz;
end

figure(2); clf
plot(sigmas)
RFcenterPx = pds.deg2px([xs ys]', trial(1).display.viewdist, trial(1).display.w2px);

[eyeData, timestamps] = io.getEyeData(thisSession);
eyeXYpx = pds.deg2px(eyeData(1:2,:), trial(1).display.viewdist, trial(1).display.w2px);
ix = timestamps > frameTimesEphysClock(1)-2 & timestamps < frameTimesEphysClock(end)+2;

out.movieFileName = movieFileName;
out.frameTimesEphysClock = frameTimesEphysClock;
out.timeMovieFile = timeShownMovieFile;
out.spikeTimes = spikes{1}.st;
out.spikeId = spikes{1}.clu;
out.cluterIds = spikes{1}.cids;
out.eyeSampleTimes = timestamps(ix);
out.eyeXYpx = eyeXYpx(:,ix);
out.viewdist = trial(1).display.viewdist;
out.screenwidth = trial(1).display.widthcm;
out.w2px = trial(1).display.w2px;
out.px2w = trial(1).display.px2w;

out.RFcenterPx = RFcenterPx;
out.RFmaskSize = sigmas'/2*trial(1).display.ppd;

disp('saving')
save([thisSession.Directory{1} '.mat'], '-v7', '-struct', 'out')
disp('done')