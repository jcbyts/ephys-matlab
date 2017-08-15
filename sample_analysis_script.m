%% Sample Analysis Script
% This script shows an example of how to load data into the workspace
% (assuming it has been properly imported)

% --- add paths
cd C:\Users\Jake\Repos\ephys-matlab\
addEphysMatlab

%% set session directory
oepath = 'C:\Data\Ellie_2017-07-31_15-21-11_Shnkd8';

% oepath = uigetdir(); % select using GUI

%% Load relevant data into the workspace
[sess, ops, info] = io.loadSession(oepath);

PDS = io.getPds(sess);

sp = io.getSpikes(ops, info);

%% plot spike waveforms?
nShanks = numel(ops);

for iShank = 1:nShanks
    figure(iShank); clf
    if isfield(sp{iShank}, 'uQ')
        clusterIds = sp{iShank}.cids(sp{iShank}.uQ>20);
    else
        clusterIds = [];
    end
    fig(iShank) = plot.spikeWaveformsFromOps(ops(iShank), sp{iShank}, 'figure', iShank, 'clusterIds', clusterIds, 'numWaveforms', 10);
end

%% detect saccades
% This is quick and dirty. We should probably replace this with something
% more robust. I doubt we want to score every saccade in the GUI, but we
% might want to score enough of them to know that we trust the algorithm.

% load Eye position
[data, timestamps, elInfo] = io.getEdf(sess, PDS);

% --- remove bad samples
ix = any(data(1:2,:)==elInfo.bitDeg(2));
data(1,ix) = nan;
data(2,ix) = nan;

% --- detect those saccades
ix = ~any(isnan(data));
[saccades] = pdsa.detectSaccades(timestamps(ix), data(1:2,ix), 'verbose', false);

% --- plot output
figure(1); clf
plot(timestamps, data(1,:), 'k'); hold on
plot(timestamps, data(2,:), 'Color', repmat(.5, 1, 3));


tbins = timestamps(1:3e3:end);
cnt = histc(saccades(1,:), tbins);
[~, id] = max(cnt);

plot([saccades(1,:); saccades(1,:)], ylim, 'r--')
xlim(tbins(id + [0 1]))
xlabel('seconds')
ylabel('degrees')
legend({'x', 'y', 'saccade start'})

figure(2); clf
plot(saccades(3,:)/elInfo.sampleRate, saccades(4,:), '.')
xlabel('duration (s)')
ylabel('amplitude (deg)')
% lfpInfo.timestamps = lfpInfo.timestamps-lfpInfo.timestamps(1);
%% Flash / Saccade triggered CSD

% pick a shank and load LFP
iShank = 1;
[lfp, lfpTime, lfpInfo] = io.getLFP(ops(iShank));

% --- flash-triggered
csdTrial = io.csdTrial(PDS);
eventTimes = [csdTrial.onset];

figure(1); clf
subplot(1,2,1)
plot.CsdBasic(lfp, eventTimes, lfpInfo)
title('Flash Triggered')
xlabel('ms')
ylabel('depth')

% --- saccade-triggered
eventTimes = saccades(1,:);
% find saccades that happened more than 200 ms after the last saccade
% (uncontaminated)
ixLong = find([0 diff(eventTimes)>.2]);

subplot(1,2,2)
plot.CsdBasic(lfp, eventTimes(ixLong), lfpInfo)
title('Saccade Triggered')
xlabel('ms')
ylabel('depth')

% Look for saccade-direction tuning in the LFP
th = cart2pol(saccades(7,:) - saccades(5,:), saccades(8,:) - saccades(6,:))'/pi*180;

figure(2); clf
et = eventTimes(ixLong)';
ev = io.convertTimeToSamples(et, lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));
th = th(ixLong(:));
thBins = 0:30:360;
ch0 = (1:32)*40;
cmap = flipud(hsv(numel(thBins)));
for i = 1:numel(thBins)-1
    
    ii = th > thBins(i) & th < thBins(i+1);
    [sta,~,bc] = pdsa.eventTriggeredAverage(lfp, ev(ii), [-300 400]);
    
    plot(bc, bsxfun(@plus, sta, ch0), 'Color', cmap(i,:)); hold on
%     plot(bc, sta+sd, '--', 'Color', cmap(i,:));
%     plot(bc, sta-sd, '--', 'Color', cmap(i,:));
    
end
xlabel('Time from saccade')

axis tight
ax2 = axes();
for i = 1:numel(thBins)-1
quiver(ax2, 0, 0, cosd(mean(thBins(i + [0 1]))), sind(mean(thBins(i + [0 1]))), 'Color', cmap(i,:), 'AutoScale', 'off'); hold on
axis off

end
xlim([-1 10])
ylim([-10 1])
title('Saccade-triggered LFP')

figure(3); clf
s = sp{iShank};
goodUnits = s.uQ > 15;
udepths = s.clusterDepths(goodUnits);
uId = s.cids(goodUnits);
[~, depthId] = sort(udepths);
clustId = uId(depthId);
nUnits = numel(clustId); % including

ax = pdsa.tight_subplot(nUnits, 1, 0.01,  0.01);
for kUnit = 1:nUnits
   st = s.st(s.clu == clustId(kUnit));
   
   set(gcf, 'currentaxes', ax(kUnit))
   for i = 1:numel(thBins)-1
       
       ii = th > thBins(i) & th < thBins(i+1);
       [sta, sd, bc] = pdsa.eventPsth(st, et(ii), [-.3 .3], .01);
       
       plot(bc, sta, 'Color', cmap(i,:)); hold on
       plot(bc, sta+sd, '--', 'Color', cmap(i,:));
       plot(bc, sta-sd, '--', 'Color', cmap(i,:));
       
   end
   plot([0 0], ylim, 'k')
   axis tight
   axis off

end



%%
for i = 1:nUnits
    set(gcf, 'currentaxes', ax(i))
    xlim([-10 50])
end



%%
stim = 'hartley';

hasStim = io.findPDScontainingStimModule(PDS, stim);

hartleyTrial = struct();
trialNum = 0;

for i = find(hasStim(:)')
    
    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 
    
    stimTrials = find(trialIx);
    
    if isempty(stimTrials)
        continue
    end
        
    kxs=PDS{i}.data{stimTrials(1)}.hartley.kxs;
    kys=PDS{i}.data{stimTrials(1)}.hartley.kys;
    
   
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
    
        kTrial = trialNum + j;
        
        hartleyTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
        hartleyTrial(kTrial).start      = hartleyTrial(kTrial).frameTimes(1);
        hartleyTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - hartleyTrial(kTrial).start;
    
        if isfield(PDS{i}.conditions{thisTrial}, stim)
            if isfield(PDS{i}.conditions{thisTrial}.(stim), 'setupRNG')
                if strcmp(PDS{i}.conditions{thisTrial}.(stim).setupRNG, 'frozenSequence')
                    hartleyTrial(kTrial).frozenSequence = true;
                    hartleyTrial(kTrial).frozenSequenceLength = PDS{i}.conditions{thisTrial}.(stim).sequenceLength;
                else
                    hartleyTrial(kTrial).frozenSequence = false;
                    hartleyTrial(kTrial).frozenSequenceLength = nan;
                end
            
            else
                hartleyTrial(kTrial).frozenSequence = false;
                hartleyTrial(kTrial).frozenSequenceLength = nan;
            end
            
        else
                hartleyTrial(kTrial).frozenSequence = false;
                hartleyTrial(kTrial).frozenSequenceLength = nan;
        end
        
        hartleyTrial(kTrial).kx         = PDS{i}.data{thisTrial}.(stim).kx;
        hartleyTrial(kTrial).ky         = PDS{i}.data{thisTrial}.(stim).ky;
        hartleyTrial(kTrial).on         = PDS{i}.data{thisTrial}.(stim).on;
        
        eyepos = io.getEyePosition(PDS{i}, thisTrial);
        hartleyTrial(kTrial).eyeSampleTime = eyepos(:,1);
        hartleyTrial(kTrial).eyeXPx        = eyepos(:,2);
        hartleyTrial(kTrial).eyeYPx        = eyepos(:,3);
        hartleyTrial(kTrial).pupilArea     = eyepos(:,4);
    end
    
    trialNum = kTrial;
    
end

%%
flipTimes = cell2mat(arrayfun(@(x) x.frameTimes, hartleyTrial, 'UniformOutput', false))';
kx        = cell2mat(arrayfun(@(x) x.kx', hartleyTrial, 'UniformOutput', false))';
ky        = cell2mat(arrayfun(@(x) x.ky', hartleyTrial, 'UniformOutput', false))';
on        = ~(isnan(kx) | isnan(ky));

kxs = unique(kx(on));
kys = unique(ky(on));

x=arrayfun(@(x) find(x==kxs), kx(on));
y=arrayfun(@(x) find(x==kys), ky(on));

ind=sub2ind([numel(kys) numel(kxs)], y, x);


%% Build design matrix

t0 = flipTimes(1);
binsize = PDS{1}.initialParametersMerged.display.ifi;
binfun = @(t) (t==0) + ceil(t/binsize);

stimOnset = binfun(flipTimes(on)-t0);
% stimulus
X = sparse(stimOnset, ind, ones(size(ind,1),1), max(stimOnset), max(ind));

ntk=30;
Xd=rfmap.makeStimRowsSparse(X, ntk);
Xd = [Xd ones(size(X,1),1)];

%% analyze RFs binned spikes
figure(10); clf

ax = pdsa.tight_subplot(nUnits,2,.001, .1, .1);

for kUnit = 1:nUnits
    
    
    st = s.st(s.clu==clustId(kUnit)) - t0;
    
    ss = binfun(st(st>0 & st < flipTimes(end)-t0));
    ss(ss > size(X,1)) = [];
    
    y = sparse(ss, ones(numel(ss), 1), ones(numel(ss), 1), size(X,1), 1);
    
    
    % STA
    % sta =
    ttsta = (Xd'*Xd + 10e2*speye(size(Xd,2)))\(Xd'*y);
    RF{kUnit} =  ttsta(1:end-1);
    sta =reshape(RF{kUnit}, ntk, []);
 
%     RF{kUnit} = fastASD(Xd(:,1:end-1), y-mean(y), [ntk max(ind)], .1);
%     sta =reshape(RF{kUnit}, ntk, []);
%     
    
    set(gcf, 'currentaxes', ax((kUnit-1)*2 + 1))
%     sta =reshape(sta(1:end-1), ntk, []);
%     
    
    [u,~,v] = svd(full(sta));
    u(:,1) = u(:,1) - mean(u(1:5,1));
    [~, im] = max(abs(u(:,1)));
    sflip = sign(u(im,1));
%     sflip = sign(sum(v(:,1)));
    
    
    spatialRF{kUnit} = reshape(sflip*v(:,1), numel(kxs), numel(kys));
    imagesc(kxs, kys, spatialRF{kUnit})
    colormap(gray.^2)
    % subplot(2,32,kUnit +32)
    axis off
    
    set(gcf, 'currentaxes', ax((kUnit-1)*2 + 2))
    
    plot((1:ntk)*binsize,sflip*u(:,1), 'k-')
    axis off
    axis tight
    drawnow
end

set(gcf, 'renderer', 'painters')
% set(ax(k), 'XTick', -1:.1:.5, 'XTickLabel', -1:.1:.5, 'TickDir', 'out')
set(gcf, 'PaperSize', [1 5], 'PaperPosition', [0 0 1 5])
% saveas(gcf, 'hartley_RFs', 'epsc')

%%
frozenTrials = find([hartleyTrial.frozenSequence]);
if any(frozenTrials)
   sequenceStarts = cell2mat(arrayfun(@(x) x.frameTimes(1:x.frozenSequenceLength:end)', hartleyTrial(frozenTrials), 'UniformOutput', false)');
   
   cmap = lines(nUnits);
   figure(11); clf
   
   for kUnit = 1:nUnits
    st = s.st(s.clu==clustId(kUnit));
   
    seqLength = hartleyTrial(frozenTrials(1)).frozenSequenceLength;
   
   [spcnt, bcenters] = pdsa.binSpTimes(st, sequenceStarts, [0 seqLength/120], 1e-3);
   
   [i, j] = find(spcnt);
   
   uIx = s.cids==unitIds(kUnit);
   
   if s.uQ(uIx)>20 && s.cR(uIx) < .2
    plot(j, i-numel(sequenceStarts)*(kUnit-1), '.', 'Color', cmap(kUnit,:)); hold on
   else
       plot(j, i-numel(sequenceStarts)*(kUnit-1), '.', 'Color', repmat(.5, 1, 3)); hold on
   end
   
   end
end

%% ploting


figure(1); clf
ax = pdsa.tight_subplot(32,2,.001, .1, .1);

for kUnit = 1:32
    set(gcf, 'currentaxes', ax((kUnit-1)*2 + 1))
%     sta =reshape(sta(1:end-1), ntk, []);
%     
    
    sta =reshape(RF{kUnit}, ntk, []);
    [u,s,v] = svd(full(sta));
    sflip = sign(sum(u(:,1)));
   
%     if 
         imagesc(kxs, kys, sflip*spatialRF{kUnit}, .5*[-1 1])
   
    
         
    colormap gray
    % subplot(2,32,kUnit +32)
    axis square
    axis off
    
    set(gcf, 'currentaxes', ax((kUnit-1)*2 + 2))
    
    plot(-(0:(ntk-1))*binsize,flipud(sflip*u(:,1)), 'k-')
    axis off
    axis tight
    axis xy
    drawnow
    
end
axis on
set(gca, 'YTick', '', 'XTick', [-.2 -.1 0], 'XTickLabel', [-.2 -.1 0], 'box', 'off')
set(gcf, 'renderer', 'painters')
% set(ax(k), 'XTick', -1:.1:.5, 'XTickLabel', -1:.1:.5, 'TickDir', 'out')
set(gcf, 'PaperSize', [1 5], 'PaperPosition', [0 0 1 5])
saveas(gcf, 'hartley_RFs', 'pdf')

%%
%% gaussian pyramid noise

stim = 'gaussianNoiseBlobs';

ppd  = PDS{1}.initialParametersMerged.display.ppd;
ifi  = PDS{1}.initialParametersMerged.display.ifi;

hasStim = io.findPDScontainingStimModule(PDS, stim);
hasFixation = cellfun(@(x) strncmp(x.initialParametersMerged.pldaps.trialFunction, 'stimuli.fixflash', 16), PDS);

hasStim = hasStim(:); % & ~hasFixation(:);

gpyrTrial = struct();
trialNum = 0;

for i = find(hasStim(:)')
    
    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 
    
    stimTrials = find(trialIx);
    
    if isempty(stimTrials)
        continue
    end
    
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
    
        kTrial = trialNum + j;
        
        % --- Timing
        gpyrTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
        gpyrTrial(kTrial).start      = gpyrTrial(kTrial).frameTimes(1);
        gpyrTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - gpyrTrial(kTrial).start;
        
        % --- Eye position
        eyepos = io.getEyePosition(PDS{i}, thisTrial);
        gpyrTrial(kTrial).eyeSampleTime = eyepos(:,1);
        gpyrTrial(kTrial).eyeXPx        = eyepos(:,2);
        gpyrTrial(kTrial).eyeYPx        = eyepos(:,3);
        gpyrTrial(kTrial).pupilArea     = eyepos(:,4);
        
        results = pdsa.detectSaccades(eyepos(:,1)', eyepos(:,2:3)'./ppd, 'verbose', false);
        
        
        
        iix = gpyrTrial(kTrial).eyeXPx < 200 | gpyrTrial(kTrial).eyeXPx > 1800;
        iiy = gpyrTrial(kTrial).eyeYPx < 100 | gpyrTrial(kTrial).eyeXPx > 1000;
        bad = iix | iiy;
       
        gpyrTrial(kTrial).eyeXPx(bad) = nan;
        gpyrTrial(kTrial).eyeYPx(bad) = nan;
        gpyrTrial(kTrial).pupilArea(bad) = nan;
        
        % find eye position on each frame
        t = gpyrTrial(kTrial).frameTimes - gpyrTrial(kTrial).start;
        tdiff = abs(bsxfun(@minus, gpyrTrial(kTrial).eyeSampleTime, t(:)')) < mean(diff(gpyrTrial(kTrial).eyeSampleTime));
        [irow,icol] = find(diff(tdiff)==-1);
        
        gpyrTrial(kTrial).saccades = false(size(gpyrTrial(kTrial).frameTimes));
        
        for iSaccade = 1:size(results,2)
            ix = t > results(1, iSaccade) & t < results(2,iSaccade);
            gpyrTrial(kTrial).saccades(ix) = true;
        end
        
        gpyrTrial(kTrial).eyePosAtFrame = [gpyrTrial(kTrial).eyeXPx(irow) gpyrTrial(kTrial).eyeYPx(irow)];
        
        gpyrTrial(kTrial).rngs       = PDS{i}.data{thisTrial}.(stim).rngs;
        gpyrTrial(kTrial).n          = PDS{i}.data{thisTrial}.(stim).n;
        gpyrTrial(kTrial).xpos       = PDS{i}.data{thisTrial}.(stim).xpos;
        gpyrTrial(kTrial).ypos       = PDS{i}.data{thisTrial}.(stim).ypos;
        gpyrTrial(kTrial).gridpos    = PDS{i}.data{thisTrial}.(stim).gridpos;
        gpyrTrial(kTrial).contrast   = repmat(PDS{i}.data{thisTrial}.(stim).n.mypars(1,:), numel(gpyrTrial(kTrial).frameTimes), 1);
        
        gpyrTrial(kTrial).xposEye    = bsxfun(@minus, gpyrTrial(kTrial).xpos, gpyrTrial(kTrial).eyePosAtFrame(:,1));
        gpyrTrial(kTrial).yposEye    = -bsxfun(@minus, gpyrTrial(kTrial).ypos, gpyrTrial(kTrial).eyePosAtFrame(:,2));
        
        
    end
    
    trialNum = kTrial;
end


%% concatenate trials

flipTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), gpyrTrial, 'UniformOutput', false)');
xposRel   = cell2mat(arrayfun(@(x) x.xposEye, gpyrTrial, 'UniformOutput', false)');
yposRel   = cell2mat(arrayfun(@(x) x.yposEye, gpyrTrial, 'UniformOutput', false)');
eyeposX   = cell2mat(arrayfun(@(x) x.eyeXPx, gpyrTrial, 'UniformOutput', false)');
eyeposY   = cell2mat(arrayfun(@(x) x.eyeYPx, gpyrTrial, 'UniformOutput', false)');
contrast  = cell2mat(arrayfun(@(x) x.contrast, gpyrTrial, 'UniformOutput', false)');
saccades  = cell2mat(arrayfun(@(x) x.saccades(:), gpyrTrial, 'UniformOutput', false)');
validTimes = ~(any(isnan(xposRel),2) | any(isnan(yposRel),2)) & ~saccades;
t0 = flipTimes(1);
tEnd = flipTimes(end);



flipTimes = flipTimes(validTimes);
xposRel   = xposRel(validTimes, :);
yposRel   = yposRel(validTimes, :);
contrast  = contrast(validTimes,:);

[flipTimes, id]=sort(flipTimes);

xposRel = xposRel(id,:);
yposRel = yposRel(id,:);
contrast = contrast(id,:);


%% Build design matrix

sc = gpyrTrial(end).n.scale*gpyrTrial(end).n.sc;
lvls = sort(unique(sc), 'descend');

xwin  = [-5 5]*ppd;
ywin  = [-5 5]*ppd;

binSize = 1; %lvls(kLevel)*37.4400;
xax = xwin(1):binSize:xwin(2);
yax = ywin(1):binSize:ywin(2);

[xx,yy]=meshgrid(xax, yax);

xax = xax/ppd;
yax = yax/ppd;

sz = ceil([diff(ywin)/binSize diff(xwin)/binSize]);

binfun = @(x) (x==0) + ceil(x/binSize);

X = [];

gkern = {};
levels = 1;
for kLevel = 1:numel(levels)
    thisLevel = lvls(levels(kLevel));
    levelIx = sc == thisLevel;
    xpos = xposRel(:,levelIx);
    ypos = yposRel(:,levelIx);
    
    c    = contrast(:,levelIx);

% ixGood = (xpos > xwin(1)) & (xpos < xwin(2)-1) & (ypos > ywin(1)) & (ypos < ywin(2)-1);



gkern{kLevel} = exp(-((xx-mean(xwin)).^2 + (yy-mean(ywin)).^2)/(2*thisLevel.^2));

xgrid = binfun(xpos - xwin(1));
ygrid = binfun(ypos - ywin(1));

xgrid(xgrid<1) = nan;
ygrid(ygrid<1) = nan;

xgrid(xgrid>sz(2)) = nan;
ygrid(ygrid>sz(1)) = nan;

gridpos = sub2ind(sz, ygrid, xgrid);

% gridpos(~ixGood)=nan;

ix = ~isnan(gridpos);

[frameNumber,gaussianNumber] = find(ix);

X = [X sparse(frameNumber, gridpos(ix), c(ix), numel(flipTimes), prod(sz))];

end

addpath l1_ls_matlab\


%% sanity check
figure(1); clf
kFrame = randi(numel(flipTimes));
tmp = sum(reshape(X(kFrame,:), [], numel(levels)),2);
imagesc(xax*ppd, yax*ppd, reshape(tmp, sz)); hold on


colormap gray
plot(xposRel(kFrame,levelIx), yposRel(kFrame,levelIx), 'or')
%% measure STA
addBias = 0;
if addBias
    Xd = [X ones(size(X,1), 1)]; %X.^2;
else
    Xd = X;
end
% C = Xd'*Xd;

figure(1); clf
nf = 1;


inds = unitIds([4 6]);
nn = numel(inds);
ax = pdsa.tight_subplot(nn, nf, .11, .1);

for kUnit = 1:nn
%
st = sp.st(sp.clu==(inds(kUnit)));
%
y = histc(st, flipTimes);
bad = diff(flipTimes)> (2*ifi);
y(bad) = 0;



% C = Xd'*Xd;

%


for ts = 1:nf
    set(gcf, 'currentaxes', ax((kUnit-1)*nf + ts))
%     subplot(nUnits,nf,(kUnit-1)*nf + ts)
    ysamp = y;
    if nf==1
        ysamp = smooth(ysamp, 3);
        ysamp(1:end-(7-1)) = ysamp(7:end);
    else
        ysamp(1:end-(ts-1)) = ysamp(ts:end);
    end
    
   
    ysamp = ysamp - mean(ysamp);
%     ysamp = ysamp ./ std(ysamp);
    w =(Xd'*ysamp);
% %     
% %     % ridge regression
% % %     if binS
%     Cprior = speye(size(Xd,2));
%     Cprior = toeplitz([-1 1 zeros(1,size(Xd,2)-3)];
%     
%     w = (C + 10e2*Cprior)\w;

%     w=l1_ls(Xd, ysamp, 8);

%     w = fastASD(X, ysamp, sz, .1);
    if addBias
        w(end)=[];
    end
    I = zeros(sz);
    for kLevel = 1:numel(gkern)
        nb = prod(sz);
        ix = (kLevel-1)*nb + (1:nb);
        I_ = reshape(w(ix), sz);
        I_ = conv2(I_, gkern{kLevel}, 'same');
        I = I + I_;
    end
    I = I - mean(I(:));
    I = sign(I) .* I.^2;
    
    imagesc(xax, yax, I); %, [-1 1]*90);
    axis xy
    colormap(gray)
    drawnow
    grid on
    set(gca, 'XTick', [-5:1:5], 'YTick', [-5:1:5], 'GridColor', 'y')
%     imagesc(xax, yax, I);
%     title(sprintf('frame -%d', ts))
end

end

