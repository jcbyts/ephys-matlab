
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

%% make sure all folders are imported
% 'Ellie_2017-07-19_15-49-43_Shnk2', ...
% 'Ellie_2017-07-20_09-41-04_Shank2D2', ...
% 'Ellie_2017-07-24_15-59-49_shank2D3', ...      
% 'Ellie_2017-07-25_13-29-24_Shank2D4', ...      
% 'Ellie_2017-07-26_13-51-39_Shank2D5', ...      
% 'Ellie_2017-07-27_14-17-50_Shank2D6', ...      
% 'Ellie_2017-07-28_15-21-04_Shank2D7', ...    


% 'Ellie_2017-07-31_15-21-11_Shnkd8', ...      
% 'Ellie_2017-08-01_11-27-24_Shank2D9', ...  

session_list = {'Ellie_2017-08-01_11-27-24_Shank2D9', ...  
    'Ellie_2017-08-02_14-17-43_ShankD10', ...      
'Ellie_2017-08-03_15-14-58_Shank2D11'};

for iSession = 1:numel(session_list)
    oepath = fullfile('C:\Data', session_list{iSession});
    fprintf('loading [%s]\n', oepath)
    ops = io.oe2dat(oepath);
end


%% import LFP for each session
for iSession = 1:numel(session_list)
   
    oepath = fullfile('C:\Data', session_list{iSession});
    fprintf('loading [%s]\n', oepath)
    
    [sess, ops, info] = io.loadSession(oepath);

    PDS = io.getPds(sess);

    [data, timestamps, elInfo] = io.getEdf(ops, PDS);

% sp = io.getSpikesFromKilo(ops, info);

    [lfp, lfpTime, lfpInfo] = io.getLFP(ops);
    
    disp('***************************************************************')
    disp('***************************************************************')
    disp('***************************************************************')
end

%% plot CSD for each session
clear CSD
sessionIx = 1:numel(session_list);
figure(1); clf
ax = pdsa.tight_subplot(1,numel(sessionIx), 0.01, 0.1, 0.1);
for i = 1:numel(sessionIx)
    iSession = session_list{sessionIx(i)};
    oepath = fullfile('C:\Data', iSession);
    
    set(gcf, 'currentaxes', ax(i))
    
    CSD(i) = session.computeCSD(oepath);
    [sess, ~, ~] = io.loadSession(oepath); 
    
    imagesc(CSD(i).time, CSD(i).depths, CSD(i).csd, 8*[-1 1]); hold on
    plot(CSD(i).time(1:CSD(i).imRescale:end), bsxfun(@plus, CSD(i).pots, CSD(i).depths(1:CSD(i).imRescale:end)), 'Color', repmat(.5, 1, 3))
    colormap jet
    title(sess.date)
    
    if i > 1
        set(gca, 'YTickLabel', [])
    else
        ylabel('Depth (\mum)')
    end
    
    xlabel('Time (ms)')
    xlim([-50 200])
    
    drawnow
end
figDir = 'C:\Users\Jake\Dropbox\Projects\MarmosetTraining';
pdsa.fixfigure(gcf, 10, [8 4]);
saveas(gcf, fullfile(figDir, 'CSDacrossSessions.pdf'))
%%
cmap = lines;

nSessions = numel(session_list);

% figure(1); clf
% ax = pdsa.tight_subplot(nSessions, 1, .01, .01);
% 
% 
% figure(2); clf
% ax2 = pdsa.tight_subplot(nSessions, 1, .01, .01);
%  figure(iSession+nSessions); clf
% 	ax = pdsa.tight_subplot(1,numel(unitIdx),.001, .001);
plotUnits = [6 6 6 5];
plotUnits = [7 7 6];
unitFileName = sprintf('unit%02.0f', plotUnits(1));

itr = 1;
figure(1); clf
figure(2); clf
for iSession = 1:nSessions
    
    %     set(gcf, 'currentaxes', ax(iSession))
    
    oepath = fullfile('C:\Data', session_list{iSession});
    [sess, ops, info] = io.loadSession(oepath);
    
    PDS = io.getPds(sess);
    
    sp = io.getSpikesFromKilo(ops, info);
    
    goodUnits = find(sp.uQ > 20);
    [~, depthIdx] = sort(sp.clusterDepths(goodUnits));
    unitIdx = sp.cids(goodUnits(depthIdx));
    unitIdx = unitIdx(plotUnits(iSession));
    
    
    
    fname = ops.fbinary;
    fid = fopen(fname, 'r');
    buffer = [32 40];
    
    if exist(ops.chanMap, 'file')
        load(ops.chanMap)
    else
        chanMap = 1:ops.Nchan;
    end
    
   
   
    for kUnit = 1:numel(unitIdx)
        
        figure(1);
        
        iix = sp.clu==unitIdx(kUnit);
        n = sum(iix);
        ss = sp.ss(iix);
        
        for i = ss(1:ceil((n/50)):n)'
            
            fseek(fid, (i-10)*2*ops.Nchan, 'bof');
            data = double(fread(fid, buffer, '*int16'));
            data = data(chanMap,:)';
            data = bsxfun(@minus, data, mean(data([1:10 (buffer(2)-10):buffer(2)],:)));
            
            wf = bsxfun(@plus, data*ops.bitVolts, 10*(sp.yc)');
            plot((1:buffer(2))+itr*buffer(2), wf, 'Color', cmap(itr,:)); hold on
            
        end
%         axis off
        ylim(sp.clusterDepths(sp.cids==unitIdx(kUnit))*10 + [-300 300]*10)
        set(gca, 'YTickLabel' , fliplr(get(gca, 'YTick')/10), 'box', 'off', 'XTick', [])
        drawnow
        
        figure(2);
%         set(gcf, 'currentaxes', ax(kUnit))
        bins = 0:5:1.5e3;
        cnt = histc(diff(double(ss)), bins); %histogram(diff(ss), 0:5:500, 'FaceColor', cmap(kUnit,:), 'EdgeColor', 'none')
        cnt = cnt'/sum(cnt);
        bar(itr*max(bins)*2 + 5 + [-bins(2:end) bins], [cnt(2:end) cnt], 'FaceColor', cmap(itr,:), 'EdgeColor', 'none'); hold on
%         xlim([0 500])
        text(itr*max(bins)*2 + 5, .038, sess.date)
        set(gca, 'YTick', [])
        title(unitIdx(kUnit))
        
        RF = hackyGetHarlteyFunction(PDS, sp, unitIdx(kUnit));
        figure(3)
        subplot(2,nSessions, iSession)
        imagesc(RF.kxs, RF.kys, RF.spatialRF);
        colormap(gray.^2)
        subplot(2,nSessions, iSession+nSessions)
        plot(RF.time, RF.temporalRF)
        axis tight
        
        itr = itr+1;
    end
    title(sess.date)
    ylim([-.1 .8])
    set(gca, 'YAxisLocation', 'right', 'YTick', 0:.2:.8)
    fclose(fid);
end

pdsa.fixfigure(1, 8, [4 4])
set(1, 'renderer', 'painters')
saveas(1, fullfile(figDir, [unitFileName '_Waveforms.pdf']), 'epsc')

pdsa.fixfigure(2, 8, [4 2])
saveas(gcf, fullfile(figDir, [unitFileName '_ISI.pdf']))

pdsa.fixfigure(3, 8, [4 4])
saveas(gcf, fullfile(figDir, [unitFileName '_RF.pdf']))


%% find csd shift amount

baseCsd = 4;

shiftSize = nan(numel(CSD),1);

for kCsd = 1:numel(CSD)
    
    
sz = size(CSD(baseCsd).csd);
csdFilt = zeros(max([sz(1)+sz(1)-1, sz(1), sz(1)]),sz(2));
nTimeBin = sz(2);
maxShift = 4;
% shift by integer channel numbers
nChannels = ceil(sz(1)/CSD(baseCsd).imRescale);

baseIdx = (maxShift+1):(nChannels-maxShift);

shifts = -maxShift:maxShift;
nShifts = numel(shifts);

csd1 = CSD(baseCsd).csd(1:CSD(baseCsd).imRescale:end,:);
csd2 = CSD(kCsd).csd(1:CSD(kCsd).imRescale:end,:);

figure(1); clf

for iShift = 1:nShifts
    
    
    shift = shifts(iShift);
    testIdx = baseIdx + shift;
    
    
    subplot(1,3,1)
    imagesc(csd1(baseIdx,:));
    subplot(1,3,2)
    imagesc(csd2(testIdx,:));
    subplot(1,3,3)
    
    a(iShift) = sum(sum(csd1(baseIdx,:) .* csd2(testIdx,:)));
    plot(shifts(1:iShift), a(1:iShift), '-o'); hold on
	pause
    
end
[~, shiftId] = max(a);
shiftSize(kCsd) = shifts(shiftId);
end

%%
    
    
    
    


for t = 1:nTimeBin
    csdFilt(:,t) = conv(CSD(baseCsd).csd(:,t), CSD(kCsd).csd(:,t), 'full');
    figure(1); clf
    plot(CSD(kCsd).csd(:,t)); hold on
    plot(CSD(baseCsd).csd(:,t));
    plot(csdFilt(:,t))
    pause
end


figure(1); clf

%%

% oepath = 'C:\Data\Ellie_3-sessions_2017-07-24_2017-07-26\';
% oepath = 'C:\Data\Ellie_2017-07-27_14-17-50_Shank2D6\';
oepath = 'C:\Data\Ellie_2017-07-31_15-21-11_Shnkd8';

oepath = uigetdir


%%
[sess, ops, info] = io.loadSession(oepath);

PDS = io.getPds(sess);

[data, timestamps, elInfo] = io.getEdf(ops, PDS);

sp = io.getSpikesFromKilo(ops, info);

[lfp, lfpTime, lfpInfo] = io.getLFP(ops);

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
[sta,~, xax] = pdsa.eventTriggeredAverage(lfp0, ev(:), [-100 600]);

% X = imresize(sta, 5);
% X = sta;
%%
figure(1); clf
% subplot(1,2,1)
ims = 5;
ch0 = -(1:32)*50;
chInd = 1:1:32;
staUp = imresize(sta(:,chInd), ims);
chans = imresize(ch0(chInd), ims);

% gkern = exp( -(1:32).^2/(2*2^2));
% staUp = filter(gkern, 1, staUp');
% staUp = flipud(staUp);
% staUp = filter(gkern, 1, staUp);
% staUp = flipud(staUp)';

time  = imresize(xax, ims);
time  = time(1,:);
chans = chans(1,:);
CSD = diff(staUp, [], 2)';

imagesc(time, chans, CSD-mean(CSD(:))); axis xy
colormap jet
hold on
plot(xax, bsxfun(@plus, sta, -(1:32)*50), 'k')

xlim([-50 200])
% 
% 
% 
% X = 0:0.1:(max(elPos)+.01);
% 
% elPos = (1:32)'*.05;
% pots2 = sta'; %((sta-mean(sta(:))) / std(sta(:)))';
% k = kCSD1d(elPos, pots2, 'X', X, 'h', .1);
% k.estimate;
% 
% subplot(1,2,2)
% imagesc(k.csdEst);
% colormap jet


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

