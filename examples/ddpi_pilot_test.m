
%% Find sessions that have movies 
meta = io.getExperimentsAnd('Subject', 'Harold');

%% load data
nSessions = size(meta, 1);

iSession = 1;

thisSession = meta(end,:);

PDS = io.getPds(thisSession);


trial = pds.getPdsTrialData(PDS);
[eyedata, timestamps, elinfo] = io.getEyeData(thisSession, PDS);

disp('Done')


%% Run offline dDPI analysis on video files
ddpidir = fullfile(getpref('EPHYS', 'SERVER_DATA'), thisSession.Directory{1});
mvl = dir(fullfile(ddpidir, 'ddpi*.avi'));
overwrite = true;
debug = false;
eyelinkOn = true;
clear track
for im = 1:numel(mvl)
    track(im) = ddpi_offline(fullfile(ddpidir, mvl(im).name), 'overwrite', overwrite, ...
        'debug', debug, ...
        'eyelinkOn', eyelinkOn, ...
        'fitPupil', true, 'pupilThresh', 2, ...
        'piThresh', 75);
end


% convert dDPI filenames to clocktime

%%
figure(1); clf
pat = '(?<id>\w+)\s+(?<date>[\d-]+)\s+(?<time>[\d-]+)\.';
info = arrayfun(@(x) regexp(x.name, pat, 'names'), mvl);
ddpidt = arrayfun(@(x) datenum([x.date '-' x.time], 'yy-mm-dd-HH-MM-SS'), info);
pdst0 = trial(1).session.initTime;
offsetdt = ddpidt - pdst0;

time = [];
eyex = [];
eyey = [];

for im = 1:numel(mvl)
    elapsedtime = datestr(offsetdt(im), 'HHMMSS');
    elapsedsec = str2double(elapsedtime(1:2))*60*60 + str2double(elapsedtime(3:4))*60 + str2double(elapsedtime(5:6));
    
    time = [time track(im).time + elapsedsec];
    bad = track(im).x1==-1;
    x = track(im).x4 - track(im).x1;
    y = track(im).y4 - track(im).y1;
    x(bad) = nan;
    y(bad) = nan;
    eyex = [eyex x];
    eyey = [eyey y];
end

pdstime = time + trial(1).session.experimentStart;
oetime = trial(1).PTB2OE(pdstime);
eyex = zscore(eyex);
eyey = zscore(eyey);
figure(1); clf
plot(oetime, eyex, '-'); hold on
% plot(oetime, eyey, '-')
x = eyedata(1,:);
x(x < -10) = nan;
ix = ~isnan(x);
plot(timestamps(ix), zscore(eyedata(1,ix)))





%% side question: is there a manifold?
x1 = [track(:).x1];
x4 = [track(:).x4];
y1 = [track(:).y1];
y4 = [track(:).y4];
xp = [track(:).xp];
yp = [track(:).yp];
xp(xp < 100) = nan;
figure(1); clf
plot3(x1, x4, xp, 'ow', 'MarkerFaceColor', 'b', 'MarkerSize', 4)
xlabel('x1'); ylabel('x4'); zlabel('xp')
% figure(2); clf
% plot3(y1, x4, yp, '.')

%%


%%





pdsdt = arrayfun(@(x) datenum(2000 + x.unique_number(1), x.unique_number(2), x.unique_number(3), x.unique_number(4), x.unique_number(5), x.unique_number(6)), trial);





[~, closestTrial] = min((pdsdt - ddpidt).^2);

figure(1); clf
plot(pdsdt, 1, '.k')
hold on
plot(pdsdt(closestTrial), 1, 'or')
plot([1 1]*ddpidt, ylim, 'r')

%%
x = track(im).x4 - track(im).x1;
y = track(im).y4 - track(im).y1;
t0 = ddpidt;
plot(track(im).time, x, track(im).time, y)


%%

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