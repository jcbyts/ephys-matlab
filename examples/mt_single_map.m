

% oepath = fullfile('C:\Data', 'Ellie_2017-08-09_13-04-23_ShankD15MT6');
oepath = 'C:\Data\Ellie_2017-08-11_11-18-43_Shank2D17MT7';
ops = io.oe2dat(oepath);

%%
[session, ops, info] = io.loadSession(oepath);

PDS = io.getPds(session);

[data, timestamps, elInfo] = io.getEdf(ops, PDS);


%%
sp = struct();
[wf, timestamps, info] = load_open_ephys_data_faster(fullfile(oepath, 'SE0.spikes'));
sp.st = timestamps;


%%

[wf, timestamps, info] = load_open_ephys_data_faster(fullfile(oepath, '110_CH36.continuous'));
[b1, a1] = butter(3, 300/30000*2, 'high');
    
   
wf = filter(b1, a1, wf);
wf = flipud(wf);
wf = filter(b1, a1, wf);
wf = flipud(wf);

wf(abs(wf) > 500) = 0;
figure(1); clf
plot(timestamps, wf)

sp.st = timestamps(wf<-100 )/30e3;

%%

stim = 'DotMapping';

ppd  = PDS{1}.initialParametersMerged.display.ppd;
ifi  = PDS{1}.initialParametersMerged.display.ifi;

hasStim = io.findPDScontainingStimModule(PDS, stim);
hasFixation = cellfun(@(x) strncmp(x.initialParametersMerged.pldaps.trialFunction, 'stimuli.fixflash', 16), PDS);

hasStim = hasStim(:); % & ~hasFixation(:);

mtmapTrial = struct();
trialNum = 0;

for i = find(hasStim(:)')
    
    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 
    
    stimTrials = find(trialIx);
    
    if isempty(stimTrials)
        continue
    end
    
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
        try
            kTrial = trialNum + j;
            
            % --- Timing
            mtmapTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
            mtmapTrial(kTrial).start      = mtmapTrial(kTrial).frameTimes(1);
            mtmapTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - mtmapTrial(kTrial).start;
            
            % --- Eye position
            eyepos = io.getEyePosition(PDS{i}, thisTrial);
            mtmapTrial(kTrial).eyeSampleTime = eyepos(:,1);
            mtmapTrial(kTrial).eyeXPx        = eyepos(:,2);
            mtmapTrial(kTrial).eyeYPx        = eyepos(:,3);
            mtmapTrial(kTrial).pupilArea     = eyepos(:,4);
            
            results = pdsa.detectSaccades(eyepos(:,1)', eyepos(:,2:3)'./ppd, 'verbose', false);
            
            
            
            iix = mtmapTrial(kTrial).eyeXPx < 200 | mtmapTrial(kTrial).eyeXPx > 1800;
            iiy = mtmapTrial(kTrial).eyeYPx < 100 | mtmapTrial(kTrial).eyeXPx > 1000;
            bad = iix | iiy;
            
            mtmapTrial(kTrial).eyeXPx(bad) = nan;
            mtmapTrial(kTrial).eyeYPx(bad) = nan;
            mtmapTrial(kTrial).pupilArea(bad) = nan;
            
            % find eye position on each frame
            t = mtmapTrial(kTrial).frameTimes - mtmapTrial(kTrial).start;
            tdiff = abs(bsxfun(@minus, mtmapTrial(kTrial).eyeSampleTime, t(:)')) < mean(diff(mtmapTrial(kTrial).eyeSampleTime));
            [irow,icol] = find(diff(tdiff)==-1);
            
            mtmapTrial(kTrial).saccades = false(size(mtmapTrial(kTrial).frameTimes));
            mtmapTrial(kTrial).saccadeTimes = results(1,:)+mtmapTrial(kTrial).start;
            
            % 5. StartX
% 6. StartY
% 7. EndX
% 8. EndY
            sacDx = results(7,:) - results(5,:);
            sacDy = results(8,:) - results(6,:);
            [th, rho] = cart2pol(sacDx, sacDy);
            mtmapTrial(kTrial).saccadeDirection = th/pi*180;
            mtmapTrial(kTrial).saccadeAmp       = rho;
            
            for iSaccade = 1:size(results,2)
                ix = t > (results(1, iSaccade)-0.05) & t < (results(2,iSaccade)+0.05);
                mtmapTrial(kTrial).saccades(ix) = true;
            end
            
            mtmapTrial(kTrial).eyePosAtFrame = [mtmapTrial(kTrial).eyeXPx(irow) mtmapTrial(kTrial).eyeYPx(irow)];
            mtmapTrial(kTrial).pupilAtFrame = mtmapTrial(kTrial).pupilArea(irow);
            
            nApert = numel(mtmapTrial(kTrial).frameTimes);
            mtmapTrial(kTrial).xpos       = PDS{i}.data{thisTrial}.(stim).x(1:nApert,:);
            mtmapTrial(kTrial).ypos       = PDS{i}.data{thisTrial}.(stim).y(1:nApert,:);
            mtmapTrial(kTrial).on         = PDS{i}.data{thisTrial}.(stim).on(1:nApert,:);
            
            mtmapTrial(kTrial).xpos(~mtmapTrial(kTrial).on) = nan;
            mtmapTrial(kTrial).ypos(~mtmapTrial(kTrial).on) = nan;
            
            mtmapTrial(kTrial).direction  = PDS{i}.data{thisTrial}.(stim).direction(1:nApert,:);
            mtmapTrial(kTrial).speed      = PDS{i}.data{thisTrial}.(stim).speed(1:nApert,:) / ifi / ppd;
            
            mtmapTrial(kTrial).dx      = mtmapTrial(kTrial).speed .* cosd(mtmapTrial(kTrial).direction);
            mtmapTrial(kTrial).dy      = mtmapTrial(kTrial).speed .* sind(mtmapTrial(kTrial).direction);
                        
            mtmapTrial(kTrial).xposEye    = bsxfun(@minus, mtmapTrial(kTrial).xpos, mtmapTrial(kTrial).eyePosAtFrame(:,1));
            mtmapTrial(kTrial).yposEye    = -bsxfun(@minus, mtmapTrial(kTrial).ypos, mtmapTrial(kTrial).eyePosAtFrame(:,2));
            
            mtmapTrial(kTrial).apertureSize = [PDS{i}.data{thisTrial}.(stim).hDots.maxRadius]/ppd;
        end
        
    end
    
    trialNum = kTrial;
end

%% plot sample trial
kTrial = 1;
figure(1); clf;
subplot(1,3,1)
plot(mtmapTrial(kTrial).xpos, mtmapTrial(kTrial).ypos, '.k')
hold on
plot(mtmapTrial(kTrial).eyePosAtFrame(:,1), mtmapTrial(kTrial).eyePosAtFrame(:,2), '.r')

subplot(1,3,2)
plot(mtmapTrial(kTrial).eyePosAtFrame(:,1), 'r'); hold on
plot(mtmapTrial(kTrial).xpos, 'k.')

subplot(1,3,3)
plot(mtmapTrial(kTrial).xposEye/ppd, mtmapTrial(kTrial).yposEye/ppd, '.k')

quiver(mtmapTrial(kTrial).xposEye/ppd, mtmapTrial(kTrial).yposEye/ppd, mtmapTrial(kTrial).dx/5, mtmapTrial(kTrial).dy/5, 'AutoScale', 'off')

%%

% build grid
xwin  = [-5 5];
ywin  = [-5 5];

binSize = 1;
xax = xwin(1):binSize:xwin(2);
yax = ywin(1):binSize:ywin(2);

[xgrid,ygrid]=meshgrid(xax, yax);

nTotalFrames = sum(arrayfun(@(x) numel(x.frameTimes), mtmapTrial));

sz = size(xgrid);

dx     = zeros(nTotalFrames, prod(sz));
dy     = zeros(nTotalFrames, prod(sz));
spcnt  = zeros(nTotalFrames, 1);



iFrame = 0;

nTrials = numel(mtmapTrial);
for kTrial = 1:nTrials
    fprintf('%d / %d \n', kTrial, nTrials)
    frameIdx = iFrame + (1:numel(mtmapTrial(kTrial).frameTimes));
    
    
    x = mtmapTrial(kTrial).xposEye/ppd;
    y = mtmapTrial(kTrial).yposEye/ppd;
    
    nGrid  = prod(sz);
    nFrame = size(x,1);
    nApert = size(x,2);
    s = unique(mtmapTrial(kTrial).apertureSize);
    
    for iDot = 1:nApert
        xd = x(:, iDot) - xgrid(:)';
        yd = y(:, iDot) - ygrid(:)';
        r = sqrt(xd.^2 + yd.^2);
        r(isnan(r)) = inf;
        idx = r < s;
        [i, j] = find(idx);
%         ind = find(idx);
        
        ind = sub2ind([nFrame nGrid], i, j);
        dxtmp = zeros(nFrame, nGrid);
        dytmp = zeros(nFrame, nGrid);
        dxtmp(ind) = mtmapTrial(kTrial).dx(i,iDot);
        dytmp(ind) = mtmapTrial(kTrial).dy(i,iDot);
        
        dx(frameIdx,:) = dx(frameIdx,:) + dxtmp;
        dy(frameIdx,:) = dy(frameIdx,:) + dytmp;
        
    end
    
    spcnt(frameIdx) = histc(sp.st, mtmapTrial(kTrial).frameTimes);
    
    iFrame = iFrame + nFrame;
    
end

%% do some regreesion
nkt = ceil(.2/ifi);
Xd = rfmap.makeStimRowsDense([dx dy], nkt);

Xd = [Xd ones(nTotalFrames,1)];
%%
saccades  = cell2mat(arrayfun(@(x) x.saccades(:), mtmapTrial, 'UniformOutput', false)');

% C = (Xd(~saccades,:)'*Xd(~saccades,:));

%%
saccades = false(size(saccades));
sta = Xd(~saccades,:)'*spcnt(~saccades);

% sta = (C + 10e7*speye(nkt*nGrid*2+1))\sta;

figure(2); clf
plot(sta(1:end-1))

I = reshape(sta(1:end-1), nkt, []);
subplot(1,2,1)
imagesc(I)

subplot(1,2,2)
imagesc(reshape(mean(Xd(~saccades,1:end-1)), nkt, []));

figure(1); clf
[u,s,v] = svd(I);
timeK  = u(:,1);
spaceK = v(:,1);
signFlip = sign(sum(timeK));

timeK = timeK * signFlip;
spaceK = spaceK * signFlip;
subplot(1,2,1)
quiver(xgrid(:), ygrid(:), spaceK(1:nGrid), spaceK(nGrid+1:end))
subplot(1,2,2)
plot( (1:nkt) * ifi, timeK)

%%

    
    





%%

saccades  = cell2mat(arrayfun(@(x) x.saccades(:), mtmapTrial, 'UniformOutput', false)');
pupil  = cell2mat(arrayfun(@(x) x.pupilAtFrame(:), mtmapTrial, 'UniformOutput', false)');
pupil(isnan(pupil)) = 0;
flipTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), mtmapTrial, 'UniformOutput', false)');
y = histc(sp.st, flipTimes);
y(y>30) = 0;
y(diff(flipTimes) > ifi) = 0;

figure(2); clf
subplot(2,1,1)
plot(flipTimes, saccades*max(y), 'k'); hold on
plot(flipTimes, y)
subplot(2,1,2)
plot(xcorr(saccades, y, 100))
%%
figure(1); clf
ev = [mtmapTrial.saccadeTimes];

th = wrapTo360([mtmapTrial.saccadeDirection]);

thBins = 0:90:360;
cmap = hsv(numel(thBins));
for i = 1:numel(thBins)-1
    
    ii = th > thBins(i) & th < thBins(i+1);
    [sta, sd, bc] = pdsa.eventPsth(sp.st, ev(ii), [-.3 .3], .01);
    
    plot(bc, sta, 'Color', cmap(i,:)); hold on
    plot(bc, sta+sd, '--', 'Color', cmap(i,:));
    plot(bc, sta-sd, '--', 'Color', cmap(i,:));
    
end
xlabel('Time from saccade')
ylabel('Firing Rate')
axis tight
ax2 = axes();
for i = 1:numel(thBins)-1
quiver(ax2, 0, 0, cosd(mean(thBins(i + [0 1]))), sind(mean(thBins(i + [0 1]))), 'Color', cmap(i,:), 'AutoScale', 'off'); hold on
axis off

end
xlim([-1 10])
ylim([-10 1])


%% concatenate trials

flipTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), mtmapTrial, 'UniformOutput', false)');
xposRel   = cell2mat(arrayfun(@(x) x.xposEye, mtmapTrial, 'UniformOutput', false)');
yposRel   = cell2mat(arrayfun(@(x) x.yposEye, mtmapTrial, 'UniformOutput', false)');
direction = cell2mat(arrayfun(@(x) x.direction, mtmapTrial, 'UniformOutput', false)');
speed     = cell2mat(arrayfun(@(x) x.speed, mtmapTrial, 'UniformOutput', false)');
on        = cell2mat(arrayfun(@(x) x.on, mtmapTrial, 'UniformOutput', false)');
saccades  = cell2mat(arrayfun(@(x) x.saccades(:), mtmapTrial, 'UniformOutput', false)');
pupil     = cell2mat(arrayfun(@(x) x.pupilAtFrame(:), mtmapTrial, 'UniformOutput', false)');
validTimes = ~(any(isnan(xposRel),2) | any(isnan(yposRel),2));
% validTimes = validTimes & (pupil > 600) & ~saccades;

y = histc(sp.st, flipTimes);
y(y>20) = 0;
y(diff(flipTimes) > ifi) = 0;

% validTimes2 = false(size(validTimes));
% k = -1:22;
% % k = k+12;
% k = k;
% sacIx = find(diff(saccades)==1)+k;
% sacIx = unique(sacIx);
% sacIx(sacIx<1 | sacIx > numel(validTimes) ) = [];
% validTimes2(sacIx) = true;

validTimes = validTimes & ~saccades;
% validTimes = validTimes & validTimes2;

% validTimes = validTimes & flipTimes >5.5e3;
% sum(validTimes)


% saccades = saccades(validTimes);
% flipTimes = flipTimes(validTimes);
% xposRel   = xposRel(validTimes, :);
% yposRel   = yposRel(validTimes, :);
% on        = on(validTimes,:);
% direction = direction(validTimes,:);
% speed     = speed(validTimes,:);
% y = y(validTimes);
% 
[flipTimes, id]=sort(flipTimes);

assert(all(id==sort(id)), 'why are flip times out of order?')




% 
figure(1); clf; 
subplot(211)
plot(y)
subplot(212)
[a,lags] = xcorr(saccades, y, 100);
plot(lags, a)


% 
% xposRel = xposRel(id,:);
% yposRel = yposRel(id,:);
% on = on(id,:);
% direction = direction(id,:);
% speed = speed(id,:);
% saccades = saccades(id);
% assume cosine tuning for direction and try to get a spatial/direction map

dx = speed .* cosd(direction);
dy = speed .* sind(direction);

xwin  = [-4 4]*ppd;
ywin  = [-10 0]*ppd;

binSize = 2*ppd; %lvls(kLevel)*37.4400;
xax = xwin(1):binSize:xwin(2);
yax = ywin(1):binSize:ywin(2);

[xx,yy]=meshgrid(xax, yax);

xax = xax/ppd;
yax = yax/ppd;

sz = size(xx); %ceil([diff(ywin)/binSize diff(xwin)/binSize]);

binfun = @(x) (x==0) + ceil(x/binSize);

xgrid = binfun(xposRel - xwin(1));
ygrid = binfun(yposRel - ywin(1));


xgrid(xgrid<1) = nan;
ygrid(ygrid<1) = nan;

xgrid(xgrid>sz(2)) = nan;
ygrid(ygrid>sz(1)) = nan;

gridpos = sub2ind(sz, ygrid, xgrid);

ix = ~isnan(gridpos);

[frameNumber,gaussianNumber] = find(ix);

Xdx = sparse(frameNumber, gridpos(ix), dx(ix), numel(flipTimes), prod(sz));
Xdy = sparse(frameNumber, gridpos(ix), dy(ix), numel(flipTimes), prod(sz));
X   = sparse(frameNumber, gridpos(ix), on(ix), numel(flipTimes), prod(sz));
% X2 = [Xdx Xdy];

% X = [X ones(size(X,1), 1)];
% X2 = [X2 ones(size(X,1), 1)];

nkt =30;
Xd = rfmap.makeStimRowsDense(X, nkt);

%%

X = [Xd ones(size(Xd,1),1)];

figure(3); clf
%%
ts = 1;
nApert = numel(ts);
y0 = smooth(y, 1)-mean(y);
for t0 = 1:nApert
    subplot(1,nApert,t0)
    t = ts(t0);
    ys = [y(t:end); y(1:(t-1))]; 

%     sta = mean(X);
    sta = (X(validTimes,:)'*ys(validTimes));
%     sta = (X'*X + 1e1*eye(size(X,2)))\sta;
    subplot(1,2,1)
    I = reshape(sta(1:end-1), nkt, []);
    imagesc(I)
    subplot(1,2,2)
    plot(mean(I,2))
    
end

%%
%     sta = (X'*X + 1e-6*eye(size(X,2)))\sta;
%     if t0>1
%     sta = mean(X);
%     end
%     plot(sta(1:end-1))
%     sta(end) = [];
    if numel(sta) == (prod(sz)+1)
        imagesc(xax, yax, reshape(sta(1:end-1), [sz(1) sz(2)]), .3*[-1 1]);
        colormap(jet)
        axis xy
        grid on
        set(gca, 'gridcolor', 'y')
    else
        I = reshape(sta(1:end-1), [sz(1) sz(2)*2]);
%         I = I(1:sz(1),1:sz(2)) + I(1:sz(1),sz(2)+1:end);
        quiver(xx(:)/ppd, yy(:)/ppd, sta(1:prod(sz)), sta(prod(sz)+1:end-1))
%         plot(I(:))
        axis tight
        drawnow
% imagesc(I, .2*[-1 1])
    end
%     axis off
    title( round(1e3*t0*ifi))
end

% try for time
nkt = 60;
Xd = rfmap.makeStimRowsDense(X*sta, nkt);

% Xd = [Xd ones(size(Xd, 1), 1)];

%
nt = numel(y);
t0 = 10;
iit = t0:nt;

timeSTA = (Xd(iit,:)'*Xd(iit,:) + 10e3*eye(nkt))\(Xd(iit,:)'*(y(1:nt-(t0-1))-mean(y(:))));
tax = fliplr(((1:nkt)-t0 )*ifi);

figure(2); clf
subplot(1,3,1)
plot(tax, timeSTA)
axis tight

% timeSTA = flipud(timeSTA);

fdx = filter(timeSTA, 1, dx);
fdy = filter(timeSTA, 1, dy);

Xdx = sparse(frameNumber, gridpos(ix), fdx(ix), numel(flipTimes), prod(sz));
Xdy = sparse(frameNumber, gridpos(ix), fdy(ix), numel(flipTimes), prod(sz));
X = [Xdx Xdy];

% X   = sparse(frameNumber, gridpos(ix), filter(timeSTA, 1, on(ix)), numel(flipTimes), prod(sz));



y0 = y - mean(y);

wls = (X'*X + 1e-5*eye(size(X,2)))\(X'*y0);

subplot(1,3,2)
if numel(wls) == (prod(sz))
    I = reshape(wls(1:end), [sz(1) sz(2)]);
    imagesc(xax, yax, I)
else
    I = reshape(wls(1:end), [sz(1) sz(2)*2]);
    I1 = reshape(wls(1:prod(sz)), [sz(1) sz(2)]);
    I2 = reshape(wls(prod(sz)+1:end), [sz(1) sz(2)]);
    
    imagesc(xax, yax, abs(I1) + abs(I2))
    
%     quiver(xx(:), yy(:), wls(1:prod(sz)), wls(prod(sz)+1:end))
end

axis xy


subplot(1,3,3)
plot(X*wls, y, '.')

%% Try with full spatiotemporal


flipTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), mtmapTrial, 'UniformOutput', false)');
xposRel   = cell2mat(arrayfun(@(x) x.xposEye, mtmapTrial, 'UniformOutput', false)');
yposRel   = cell2mat(arrayfun(@(x) x.yposEye, mtmapTrial, 'UniformOutput', false)');
direction = cell2mat(arrayfun(@(x) x.direction, mtmapTrial, 'UniformOutput', false)');
speed     = cell2mat(arrayfun(@(x) x.speed, mtmapTrial, 'UniformOutput', false)');
on        = cell2mat(arrayfun(@(x) x.on, mtmapTrial, 'UniformOutput', false)');
saccades  = cell2mat(arrayfun(@(x) x.saccades(:), mtmapTrial, 'UniformOutput', false)');
pupil     = cell2mat(arrayfun(@(x) x.pupilAtFrame(:), mtmapTrial, 'UniformOutput', false)');
validTimes = ~(any(isnan(xposRel),2) | any(isnan(yposRel),2));
% validTimes = validTimes & (pupil > 600) & ~saccades;

y = histc(sp.st, flipTimes);
y(y>20) = 0;
y(diff(flipTimes) > ifi) = 0;

% validTimes2 = false(size(validTimes));
% k = -1:22;
% % k = k+12;
% k = k;
% sacIx = find(diff(saccades)==1)+k;
% sacIx = unique(sacIx);
% sacIx(sacIx<1 | sacIx > numel(validTimes) ) = [];
% validTimes2(sacIx) = true;

validTimes = validTimes & ~saccades;
% validTimes = validTimes & validTimes2;

% validTimes = validTimes & flipTimes >5.5e3;
% sum(validTimes)

% 
% saccades = saccades(validTimes);
% flipTimes = flipTimes(validTimes);
% xposRel   = xposRel(validTimes, :);
% yposRel   = yposRel(validTimes, :);
% on        = on(validTimes,:);
% direction = direction(validTimes,:);
% speed     = speed(validTimes,:);
% y = y(validTimes);
% 
[flipTimes, id]=sort(flipTimes);

assert(all(id==sort(id)), 'why are flip times out of order?')




% 
figure(1); clf; 
subplot(211)
plot(y)
subplot(212)
[a,lags] = xcorr(saccades, y, 100);
plot(lags, a)

%%
% 
% xposRel = xposRel(id,:);
% yposRel = yposRel(id,:);
% on = on(id,:);
% direction = direction(id,:);
% speed = speed(id,:);
% saccades = saccades(id);
% assume cosine tuning for direction and try to get a spatial/direction map

dx = speed .* cosd(direction);
dy = speed .* sind(direction);

xwin  = [-10 10]*ppd;
ywin  = [-10 10]*ppd;

binSize = 1*ppd; %lvls(kLevel)*37.4400;
xax = xwin(1):binSize:xwin(2);
yax = ywin(1):binSize:ywin(2);

[xx,yy]=meshgrid(xax, yax);

xax = xax/ppd;
yax = yax/ppd;

sz = size(xx); %ceil([diff(ywin)/binSize diff(xwin)/binSize]);

binfun = @(x) (x==0) + ceil(x/binSize);

xgrid = binfun(xposRel - xwin(1));
ygrid = binfun(yposRel - ywin(1));


xgrid(xgrid<1) = nan;
ygrid(ygrid<1) = nan;

xgrid(xgrid>sz(2)) = nan;
ygrid(ygrid>sz(1)) = nan;

gridpos = sub2ind(sz, ygrid, xgrid);

ix = ~isnan(gridpos);

[frameNumber,gaussianNumber] = find(ix);

Xdx = sparse(frameNumber, gridpos(ix), dx(ix), numel(flipTimes), prod(sz));
Xdy = sparse(frameNumber, gridpos(ix), dy(ix), numel(flipTimes), prod(sz));
X   = sparse(frameNumber, gridpos(ix), on(ix), numel(flipTimes), prod(sz));
X2 = [Xdx Xdy];

rfmap.makeStimRowsDense(X2, 


