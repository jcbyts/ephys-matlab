
addpath(genpath('C:\Users\Jake\Repos\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\Jake\Repos\npy-matlab')) % path to npy-matlab scripts
addpath(genpath('C:\Users\Jake\Repos\spikes')) % path to npy-matlab scripts
addpath(genpath('C:\Users\Jake\Repos\sortingQuality')) % path to npy-matlab scripts
addpath(genpath('C:\Users\Jake\Dropbox\MatlabCode\Repos\RigBuildPhotodiodeTest\analysis-tools-master'))
addpath(genpath('C:\Users\Jake\Dropbox\MatlabCode\Repos\pdstools'))

addpath C:\Users\Jake\Dropbox\MatlabCode\Repos\scalablerf\code
addpath C:\Users\Jake\Dropbox\MatlabCode\Repos\scalablerf\code\code_fastASD
addpath C:\Users\Jake\Dropbox\MatlabCode\Repos\scalablerf\code\tools_optim
addpath C:\Users\Jake\Dropbox\MatlabCode\Repos\ncclabcode\tools_dft\
addpath C:\Users\Jake\Dropbox\MatlabCode\Repos\ncclabcode\tools_kron\

addpath('C:\Users\Jake\Dropbox\MatlabCode\Repos\pds-stimuli')

% add PLDAPS (https://github.com/HukLab/PLDAPS.git)
addpath(genpath('C:\Users\Jake\Dropbox\MatlabCode\Repos\PLDAPS'))

% add PLDAPStools (https://github.com/jonaskn/PLDAPStools.git)
addpath('C:\Users\Jake\Dropbox\MatlabCode\Repos\PLDAPStools\')

% add edfmex for reading eyelink files (https://github.com/HukLab/edfmex)
addpath('C:\Users\Jake\Dropbox\MatlabCode\Repos\edfmex\')

%%
% oepath = 'C:\Data\Ellie_3-sessions_2017-07-24_2017-07-26\';
% oepath = 'C:\Data\Ellie_2017-07-27_14-17-50_Shank2D6\';
oepath = 'C:\Data\Ellie_2017-07-31_15-21-11_Shnkd8';

oepath = uigetdir


%%
[session, ops, info] = io.loadSession(oepath);

info.timestamps = info.timestamps/info.sampleRate;
info.timestamps = info.timestamps*datenum(0000,00,00,00,00,1) + info.dateNum;

PDS = io.getPds(session);

%%
[data, timestamps, elInfo] = io.getEdf(ops, PDS);

sp = io.getSpikesFromKilo(ops, info);

[lfp, lfpTime, lfpInfo] = io.getLFP(ops);

%% plot spike waveforms?
[~, depthIdx] = sort(sp.clusterDepths);

clustId = sp.cids(depthIdx);
cmap = lines;

figure(1); clf
fname = ops.fbinary;
fid = fopen(fname, 'r');
buffer = [32 50];

if exist(ops.chanMap, 'file')
    load(ops.chanMap)
else
    chanMap = 1:ops.Nchan;
end

for kClust = 1:numel(clustId)
    
iix = sp.clu==clustId(kClust);
n = sum(iix);
ss = sp.ss(iix);

for i = ss(1:ceil((n/500)):n)'
       
    fseek(fid, (i-10)*2*ops.Nchan, 'bof');
    data = double(fread(fid, buffer, '*int16'));
    data = data(chanMap,:)';
    data = bsxfun(@minus, data, mean(data([1:10 (buffer(2)-10):buffer(2)],:)));
    
    wf = bsxfun(@plus, data*ops.bitVolts, 5*flipud(sp.yc)');
    plot((1:buffer(2))+kClust*buffer(2), wf, 'Color', cmap(kClust,:)); hold on
    drawnow
end

end

fclose(fid);







%%
% [data, timestamps, elInfo] = io.getEdf(ops, PDS);

ix = any(data(1:2,:)==elInfo.bitDeg(2));
data(1,ix) = nan;
data(2,ix) = nan;
figure(1); clf
plot(timestamps, data(1:2,:)')

ix = ~any(isnan(data));
[saccades] = pdsa.detectSaccades(timestamps(ix), data(1:2,ix), 'verbose', false);
% lfpInfo.timestamps = lfpInfo.timestamps-lfpInfo.timestamps(1);
%% Saccade-Triggered CSD
et = saccades(1,:);

ixLong = [0 diff(et)>.2];

theta = cart2pol(saccades(7,:)-saccades(5,:), saccades(8,:)-saccades(6,:));
theta = theta*180/pi;
theta = wrapTo360(theta);

ixTheta = theta > 0 & theta < 360;
et = et(ixTheta & ixLong);


ev = io.convertTimeToSamples(et, lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));

[sta,~, xax] = pdsa.eventTriggeredAverage(lfp, ev(:), [-200 1e3]);

% X = imresize(sta, 5);
% X = sta;

figure(2); clf
subplot(1,2,1)
ims = 15;
ch0 = -(1:32)*50;
chInd = 1:1:32;
staUp = imresize(sta(:,chInd), ims);
chans = imresize(ch0(chInd), ims);
time  = imresize(xax, ims);
time  = time(1,:);
chans = chans(1,:);
CSD = diff(staUp, [], 2)';

imagesc(time, chans, CSD-mean(CSD(:))); axis xy
colormap jet
hold on
plot(xax, bsxfun(@plus, sta, -(1:32)*50), 'k')

xlim([-150 200])

subplot(1,2,2)
CSD = csd.splineCSD(sta', 'el_pos', (1:32)*.05);
imagesc(xax, ch0, CSD)

%% Try CSD analysis
stim = 'csdFlash';

hasStim = io.findPDScontainingStimModule(PDS, stim);

csdTrial = struct();
trialNum = 0;

for i = find(hasStim(:)')
    
    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 

    stimTrials = find(trialIx);
    
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
    
        kTrial = trialNum + j;
        
        csdTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
        csdTrial(kTrial).start      = csdTrial(kTrial).frameTimes(1);
        csdTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - csdTrial(kTrial).start;
    
        csdTrial(kTrial).on         = PDS{i}.data{thisTrial}.(stim).on;
        csdTrial(kTrial).onset      = csdTrial(kTrial).frameTimes(diff(csdTrial(kTrial).on)==1);
    end
    
    trialNum = kTrial;
    
end

%%
et = [csdTrial.onset];
sacTimes = saccades(1,:);

isClean = false(numel(et), 1);
for iEv = 1:numel(et)
   sdiff = sacTimes - et(iEv);
   if ~any(sdiff > -.1 & sdiff < 0.1)
       isClean(iEv) = true;
   else
       isClean(iEv) = false;
   end
end
sum(isClean)


ev = io.convertTimeToSamples(et(isClean), lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));

%%
[b, a] = butter(1, 50/lfpInfo.sampleRate*2, 'low');

lfp0 = filtfilt(b,a,lfp);

%%
[sta,~, xax] = pdsa.eventTriggeredAverage(lfp0, ev(:), [-100 250]);

% X = imresize(sta, 5);
% X = sta;
%%
figure(1); clf
ch0 = fliplr((1:32)*50);

subplot(1,4,1)
chInd = 1:1:32;
ims = 5;
staUp = imresize(sta(:,chInd), ims);
chans = imresize(ch0(chInd), ims);

time  = imresize(xax, ims);
time  = time(1,:);
chans = chans(1,:);
CSD = diff(staUp, [], 2)';
imagesc(time, chans, CSD);
title('2nd spatial derivative')
hold on
plot(xax, bsxfun(@plus, sta, ch0), 'Color', repmat(.5, 1, 3))


subplot(1,4,2)
CSD = csd.standardCSD(sta', 'el_pos', ch0/1e3); colormap jet
imagesc(xax, ch0, CSD);
hold on
plot(xax, bsxfun(@plus, sta, ch0), 'Color', repmat(.5, 1, 3))
title('Standard')

subplot(1,4,3)
CSD = csd.stepCSD(sta', 'el_pos', ch0/1e3);
imagesc(xax, ch0, CSD);
hold on
plot(xax, bsxfun(@plus, sta, ch0), 'Color', repmat(.5, 1, 3))
title('inverse CSD method')


subplot(1,4,4)
CSD = csd.stepCSD(sta', 'el_pos', ch0/1e3);
imagesc(xax, ch0, CSD);
hold on
plot(xax, bsxfun(@plus, sta, ch0), 'Color', repmat(.5, 1, 3))
title('iCSD spline')


%% spikes
% nUnits = numel(sp.cids);
% for i = 1:nUnits
% %     subplot(nUnits, 1, i)
%     pdsa.plotRaster(sp.st(sp.clu==sp.cids(i)), et, [-.1 .5], 10);
%     pause
% end
% %% spike triggered
% nUnits = numel(sp.cids);
% figure(3); clf
% ax = pdsa.tight_subplot(1, nUnits, .01, .01);
% for i = 1:nUnits
% %     subplot(nUnits, 1, i)
% set(gcf, 'currentaxes', ax(i))
%     st = sp.st(sp.clu==sp.cids(i));
%     ev = io.convertTimeToSamples(st, lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));
%     
%     [sta,~, xax] = pdsa.eventTriggeredAverage(lfp, ev(:), [-100 500]);
%     
%     ims = 15;
%     ch0 = -(1:32)*50;
%     chInd = 1:1:32;
%     staUp = imresize(sta(:,chInd), ims);
%     chans = imresize(ch0(chInd), ims);
%     time  = imresize(xax, ims);
%     time  = time(1,:);
%     chans = chans(1,:);
%     CSD = diff(staUp, [], 2)';
%     
%     imagesc(time, chans, CSD); axis xy
%     colormap jet
%     hold on
%     plot(xax, bsxfun(@plus, sta*5, -(1:32)*50), 'k')
%     
%     xlim([-50 100])
%     drawnow
% end

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
    
%     if any(cellfun(@(x) any(strcmp(fieldnames(x), stim)), PDS{i}.conditions))
%         keyboard
%     end
        

    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 
    
    stimTrials = find(trialIx);
    
    if isempty(stimTrials)
        continue
    end
        
    kxs=PDS{i}.data{stimTrials(1)}.hartley.kxs;
    kys=PDS{i}.data{stimTrials(1)}.hartley.kys;
    
    % check that the stimuli were the same grid the whole time
%     assert(size(unique(cell2mat(cellfun(@(x) x.(stim).kxs, PDS{i}.data(stimTrials), 'UniformOutput', false)'), 'rows'), 1) ==1, ...
%     'kx grid changed!')
% 
%     assert(size(unique(cell2mat(cellfun(@(x) x.(stim).kys, PDS{i}.data(stimTrials), 'UniformOutput', false)'), 'rows'), 1) ==1, ...
%     'ky grid changed!')

    
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


%% example trial
figure(10); clf

kTrial = (numel(hartleyTrial))-2;

stairs(hartleyTrial(kTrial).frameTimes-hartleyTrial(kTrial).start, 10*hartleyTrial(kTrial).on+15, 'k'); hold on
plot(hartleyTrial(kTrial).eyeSampleTime, (hartleyTrial(kTrial).eyeXPx-960)/37, 'k')
plot(hartleyTrial(kTrial).eyeSampleTime, (hartleyTrial(kTrial).eyeYPx-540)/37, 'Color', repmat(.5, 1, 3))
% xlim([10 15])
ylim([-10 40])

set(gcf, 'renderer', 'painters')
% set(ax(k), 'XTick', -1:.1:.5, 'XTickLabel', -1:.1:.5, 'TickDir', 'out')
set(gcf, 'PaperSize', [4 2], 'PaperPosition', [0 0 4 2])
% saveas(gcf, 'hartley_trial', 'epsc')

figure(11); clf
ix = abs((hartleyTrial(kTrial).eyeXPx-960)/37) < 10;
% ix = 2e3:4e3;
plot((hartleyTrial(kTrial).eyeXPx(ix)-960)/37,(hartleyTrial(kTrial).eyeYPx(ix)-540)/37, 'r')
xlim([-15 15])
ylim([-15 18])
set(gcf, 'renderer', 'painters')
% set(ax(k), 'XTick', -1:.1:.5, 'XTickLabel', -1:.1:.5, 'TickDir', 'out')
set(gcf, 'PaperSize', [4 4], 'PaperPosition', [0 0 4 4])
% saveas(gcf, 'hartley_trialEye', 'epsc')

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
nUnits = numel(sp.cids);
ax = pdsa.tight_subplot(nUnits,2,.001, .1, .1);

[~, id] = sort(sp.clusterDepths);
unitIds = sp.cids(id);
for kUnit = 1:nUnits
    
    
    st = sp.st(sp.clu==unitIds(kUnit)) - t0;
    
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
    
    [u,s,v] = svd(full(sta));
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
saveas(gcf, 'hartley_RFs', 'epsc')

%%
frozenTrials = find([hartleyTrial.frozenSequence]);
if any(frozenTrials)
   sequenceStarts = cell2mat(arrayfun(@(x) x.frameTimes(1:x.frozenSequenceLength:end)', hartleyTrial(frozenTrials), 'UniformOutput', false)');
   
   cmap = lines(nUnits);
   figure(11); clf
   
   for kUnit = 1:nUnits
    st = sp.st(sp.clu==unitIds(kUnit));
   
    seqLength = hartleyTrial(frozenTrials(1)).frozenSequenceLength;
   
   [spcnt, bcenters] = pdsa.binSpTimes(st, sequenceStarts, [0 seqLength/120], 1e-3);
   
   [i, j] = find(spcnt);
   
   uIx = sp.cids==unitIds(kUnit);
   
   if sp.uQ(uIx)>20 && sp.cR(uIx) < .2
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

