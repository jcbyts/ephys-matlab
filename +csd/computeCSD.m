function CSD = computeCSD(ops, varargin)


% check if argument 1 is a path or an ops
if ~isstruct(ops)    
    fdir = ops;
    ops = io.loadOps(fdir);
end

ip = inputParser();
ip.addParameter('method', 'standard')
ip.addParameter('window', [-.1 .25])
ip.addParameter('binSize', 1e-3)
ip.addParameter('verbose', true)

ip.parse(varargin{:})

[session, ops, ~] = io.loadSession(ops.root);

PDS = io.getPds(session);

[lfp, timestamps, lfpInfo] = io.getLFP(ops);

% --- find CSD flash trials
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

% --- compute stimulus-triggered average field potentials
et = [csdTrial.onset];
isClean = true(numel(et),1);

% sacTimes = saccades(1,:);
% 
% isClean = false(numel(et), 1);
% for iEv = 1:numel(et)
%    sdiff = sacTimes - et(iEv);
%    if ~any(sdiff > -.1 & sdiff < 0.1)
%        isClean(iEv) = true;
%    else
%        isClean(iEv) = false;
%    end
% end
% sum(isClean)


ev = io.convertTimeToSamples(et(isClean), lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));


[sta,~, xax] = pdsa.eventTriggeredAverage(lfp, ev(:), ip.Results.window/ip.Results.binSize);

% --- compute CSD
ch0 = fliplr((1:32)*50);

chInd = 1:1:32;
ims = 5;
staUp = imresize(sta(:,chInd), ims);
chans = imresize(ch0(chInd), ims);

time  = imresize(xax, ims);
time  = time(1,:);
chans = chans(1,:);
csd = diff(staUp, [], 2)';

CSD.method = ip.Results.method;
CSD.time   = time;
CSD.depths = chans;
CSD.csd    = csd;
CSD.imRescale = ims;
CSD.pots      = sta;
CSD.potsUp    = staUp;
% 
% imagesc(time, chans, csd);
% title('2nd spatial derivative')
% hold on
% plot(xax, bsxfun(@plus, sta, ch0), 'Color', repmat(.5, 1, 3))
