function [OE2PTBfit, PTB2OEfit,PTB2OE, maxreconstructionerror ] = sync2OeClock_trial(PDS, filenameE)
% SYNCOPENEPHYSCLOCK synchronizes the pdlaps PTB clock with open-ephys recording
% Inputs:
%   PDS       - PDS struct (or cell-array of PDS structs)
%   filenameE - string (or cell-array of strings) path to OE events file (*.kwe, *.events)
% Outputs:
%   OE2PTBfit - parameters of fit
%   OE2PTB@function_handle - converts OE time to PTB time
%   PTB2OE@function_handle - convers PTB time to OE time
%   maxreconstructionerror - maximum error in the synchonization
% Example Call:
%   [OE2PTBfit, OE2PTB,PTB2OE, maxreconstructionerror ] = syncOpenEphysClock(PDS, filenameE)

% 2017.08.14     jly     modified Jonas' version
verbose = true;
debug = false;

if ~iscell(filenameE)
    filenameE={filenameE};
end

nOeFiles=numel(filenameE);

bitNumber=[];
timestamps=[];
highlow=[];

for kOeFile=1:nOeFiles
    [~, ~, ext]=fileparts(filenameE{kOeFile});
    switch ext
        case '.kwe'
            tmp_timestamps = hdf5read(filenameE{kOeFile}, '/event_types/TTL/events/time_samples');
            tmp_highlow = hdf5read(filenameE{kOeFile}, '/event_types/TTL/events/user_data/eventID');
            tmp_bitNumber = hdf5read(filenameE{kOeFile}, '/event_types/TTL/events/user_data/event_channels');
            timestamps=[timestamps; tmp_timestamps(:)]; %#ok<*AGROW>
            bitNumber=[bitNumber; tmp_bitNumber(:)];
            highlow=[highlow; tmp_highlow(:)];
        case '.events'
            [tmp_bitNumber, tmp_timestamps, info]=load_open_ephys_data_faster(filenameE{kOeFile});
            tmp_highlow=info.eventId;
            timestamps=[timestamps; tmp_timestamps(:)]; %#ok<*AGROW>
            bitNumber=[bitNumber; tmp_bitNumber(:)];
            highlow=[highlow; tmp_highlow(:)];
    end
end
%%

% --- find when the strobe signal occured
strobeSet   = find(bitNumber==7 & highlow==1); % strobe came on
strobeUnset = find(bitNumber==7 & highlow==0); % strobe came off
strobeUnset = [1; strobeUnset];

% extract strobed words
value=nan(size(strobeSet));
for iStrobe=1:length(strobeSet)

    % look between the last unset and the next unset
    twin = timestamps([strobeUnset(iStrobe) strobeUnset(iStrobe+1)]);

    % find all times
    ts = timestamps >= twin(1) & timestamps < twin(2);
    
    % throw out the strobe bit and only keep valid flipped bits
    ts = ts & highlow==1 & bitNumber ~= 7;

    bn = bitNumber(ts);
    bn = unique(bn); % this is unlikely
       
    value(iStrobe)=sum(2.^bn);  
    
    if debug
        figure(1); clf; plot(timestamps(highlow==1), bitNumber(highlow==1), '.'); hold on; plot(timestamps(highlow==0), bitNumber(highlow==0), '.'); hold on
        xlim(timestamps(strobeSet(iStrobe)) + [-.02 .02])
    
        ylim([0 8])
    
        plot(timestamps(strobeSet(iStrobe))*[1 1], ylim, 'k--')
        plot(twin(1)*[1 1], ylim, 'b--')
        plot(twin(2)*[1 1], ylim, 'b--')
        plot(timestamps(ts), bitNumber(ts), 'o')
        title(iStrobe)
        drawnow
        
        disp(bn')
        disp(value(iStrobe))
    end
     
     
end

if debug
clf
plot(value, '.')
end

%%

sixletsOE=fliplr(conv2(value,eye(6)));
sixletsOE=sixletsOE(6:end,:);

sixletsOEts=fliplr(conv2(timestamps(strobeSet),eye(6)));
sixletsOEts=sixletsOEts(6:end,:);

sixletsOE7=fliplr(conv2(value,eye(7)));
sixletsOE7=sixletsOE7(7:end,:);

% sixletsOE7ts=fliplr(conv2(timestamps(strobeSet),eye(7)));
% sixletsOE7ts=sixletsOE7ts(7:end,:);

sixletsPTB=cell2mat(arrayfun(@(X) X.unique_number, PDS,'uni',0));
sixletsPTB = mod(sixletsPTB, 2^7);
sixletsPTBts = cell2mat(arrayfun(@(X) X.datapixx.unique_number_time(:,1)', PDS,'uniformOutput',false));

% sort PTB trials so they all came in order
[~, ind] = sort(sixletsPTBts(:,1));
sixletsPTBts = sixletsPTBts(ind,:);
sixletsPTB = sixletsPTB(ind,:);

sixletsOEts_ = sixletsOEts;
sixletsPTBts_ = sixletsPTBts;


% if using datetime as the unique wor then the first two columns will be
% the year and month. This lets us reduce the number of possible options
% substantially and try to catch any errors
if all(sixletsPTB(:,1)==mode(sixletsPTB(:,1))) && all(sixletsPTB(:,2)==mode(sixletsPTB(:,2)))
    year = mode(sixletsPTB(:,1));
    month = mode(sixletsPTB(:,2));
    day = mode(sixletsPTB(:,3));

    putativeUniqueNumbers = sixletsOE(:,1) == year & sixletsOE(:,2) == month & sixletsOE(:,3) == day;

    nOEtrials = sum(putativeUniqueNumbers);

    nPTBtrials = size(sixletsPTB,1);
    
    if verbose
        fprintf('Datetime was used as unique word.\n')
        fprintf('Found %d PTB trials and %d strobes on the OE recording\n', nPTBtrials, nOEtrials)
    end
    
    sixletsOE = sixletsOE(putativeUniqueNumbers,:);
    sixletsOEts = sixletsOEts(putativeUniqueNumbers,:);
    
    sixletsOE7 = sixletsOE7(putativeUniqueNumbers,:);
%     sixletsOE7ts = sixletsOE7ts(putativeUniqueNumbers,:);

    % offset year by 2000 (this is future proof until the year 3000)
    sixletsOE(:,1) = sixletsOE(:,1)+2e3;
    sixletsPTB(:,1) = sixletsPTB(:,1)+2e3;
end

[~, hind] = ismember(sixletsPTB, sixletsOE, 'rows');

goodPTBindex = hind~=0;
if verbose
    fprintf('%d strobes match\n', sum(goodPTBindex))
end

ptbTs = sixletsPTBts(goodPTBindex,:);
oeTs  = sixletsOEts(hind(goodPTBindex),:);
ptbTs = ptbTs(:);
oeTs = oeTs(:);

OE2PTBfit=[oeTs(:) ones(numel(oeTs),1)]\ptbTs(:);
% OE2PTB=@(x) x*OE2PTBfit(1) + OE2PTBfit(2);
PTB2OE=@(x) (x - OE2PTBfit(2))/OE2PTBfit(1);

if mean(goodPTBindex) > 0.5
    PTB2OEfit = OE2PTBfit;
    d = oeTs - PTB2OE(ptbTs);

    plot(oeTs, d*1e3, '.')
    ylim([-1 1])
    ylabel('error (ms)')
    xlabel('oe timestamp')
    OE2PTBfit = [];
    maxreconstructionerror = max(abs(d));
    fprintf('%d strobes match\n', sum(goodPTBindex))
    fprintf('exiting\n')
    return
end

    %% over half the strobes are missing
if mean(goodPTBindex) < 0.5 && isempty(intersect(unique(bitNumber), 5))
    % the 6th bit wasn't included. This affected a few datasets where the
    % wire physically came unplugged. Luckily, we used the datetime as our
    % unique words and we can therefor reconstruct where the corrupted
    % strobes are. Because bit 6 was missing, all minutes > 32 or seconds
    % greater than 32 are corrupt. If we only take strobes that occured
    % in the first half of the hour, we should be okay
    
    figure(1); clf
    subplot(1,2,1)
    plot(sixletsOEts, datenum(sixletsOE), '.b'); hold on
    subplot(1,2,2)
    plot(sixletsPTBts, datenum(sixletsPTB), 'b.'); hold on
    
    warning('sync2OeClock: bit 5 missing. Trying to repair parts of the clock');
    
    nOEtrials = size(sixletsOE,1);
    validOEindex = false(nOEtrials, 1);
    
    % find any times the hour changed
    hourChange = find(diff(sixletsOE(:,4)));
    minuteChange = find(diff(sixletsOE(:,5)) < 0 );
    
%     minuteChange(find(diff(minuteChange) == 1)+1)
    
    k = find(hourChange == minuteChange);
    n = numel(k);
    for i = 1:n
        if numel(minuteChange) == k(i)
            validOEindex((minuteChange(k(i))+1):end) = true;
        else
            validOEindex((minuteChange(k(i))+1):minuteChange(k(i)+1)) = true;
        end
    end
    
    if ~any(k==1)
        validOEindex(1:minuteChange(1)) = true;
    end
    
    sixletsOE(~validOEindex,5) = sixletsOE(~validOEindex,5)+32;
%     validOEindex(abs(detrend(datenum(sixletsOE))) > .01) = false;
 
    subplot(1,2,1)
    plot(sixletsOEts, datenum(sixletsOE), '.r'); hold on
    
    % include trial info this time
    [~, hind] = ismember([sixletsPTB(:,1:5) arrayfun(@(x) x.pldaps.iTrial, PDS)],[sixletsOE(:,1:5) sixletsOE7(:,end)], 'rows');

    goodPTBindex = hind~=0;
    if verbose
        fprintf('%d strobes match\n', sum(goodPTBindex))
    end
    
end

%% take matched trials and try to recover the reward timestamps as well
ptbIndex = find(goodPTBindex(:)');
% trialHadReward = ptbIndex(arrayfun(@(x) ~isempty(x.reward.log), PDS(ptbIndex)));
trialHadReward = ptbIndex(arrayfun(@(x) ~isempty(x.behavior.reward.timeReward), PDS(ptbIndex)));

oeRewardTimes = [];
pdsRewardTimes = [];

for iTrial = 1:numel(trialHadReward)
    thisTrial = trialHadReward(iTrial);

%     assert(all(sixletsPTB(thisTrial,:) == sixletsOE(hind(thisTrial),:)), 'trial does not match')

    twin = sixletsOEts(hind(thisTrial), end);
    twin = twin + [0.01 diff(PDS(thisTrial).timing.flipTimes(1,[1 end]))]; % cutoff 10ms from either end

    ts = timestamps(strobeSet);
    ix = ts > twin(1) & ts < twin(2);
    
    if ~any(ix)
        continue
    end
    
    oe_reward_times = ts(ix);
    oe_reward_values = value(ix);

    goodRewardIndex = oe_reward_values == PDS(thisTrial).pldaps.iTrial;
    oe_reward_times  = oe_reward_times(goodRewardIndex);
    pds_reward_times = PDS(thisTrial).behavior.reward.timeReward(1,:);

    if debug
    figure(1); clf; plot(timestamps(highlow==1), bitNumber(highlow==1), '.'); hold on; plot(timestamps(highlow==0), bitNumber(highlow==0), '.'); hold on; xlim(twin)
    plot([1; 1]*ts(:)', ylim'*ones(1,numel(ts)), 'k--')
    end

    if debug
        figure(2); clf; plot(diff(oe_reward_times), '.'); hold on; plot(diff(pds_reward_times), '.')
    end
    
    if numel(oe_reward_times) == numel(pds_reward_times)
        oeRewardTimes = [oeRewardTimes; oe_reward_times(:)];
        pdsRewardTimes = [pdsRewardTimes; pds_reward_times(:)];
    end
end

figure(1); clf
plot(oeRewardTimes, pdsRewardTimes, '.')
%%

ptbTs = sixletsPTBts(goodPTBindex,:);
oeTs  = sixletsOEts(hind(goodPTBindex),:);

ptbTs = [ptbTs(:); pdsRewardTimes(:)];
oeTs  = [oeTs(:); oeRewardTimes(:)];

goodTs = abs(zscore(oeTs - ptbTs)) < 8;

ptbTs = ptbTs(goodTs);
oeTs  = oeTs(goodTs);

pdsFileStarts = find(arrayfun(@(x) x.pldaps.iTrial==1, PDS));
pdsFileStops  = [pdsFileStarts(2:end) + 1; numel(PDS)];
pdsTrialStarts = arrayfun(@(x) x.timing.flipTimes(1), PDS);

splineKnots = pdsTrialStarts(sort([pdsFileStarts; pdsFileStops]));
% splineKnots(1) = splineKnots(1) - 1;
% splineKnots(1) = splineKnots(end) + 1;
PTB2OEfit = splinefit(ptbTs, oeTs, splineKnots, 2);
PTB2OE = @(x) ppval(PTB2OEfit, x);

figure(1); clf
subplot(1,2,1)
plot(ptbTs, oeTs, '.'); hold on
xlabel('PTB clock')
ylabel('OE clock')


xx = linspace(min(ptbTs), max(ptbTs), 100);
plot(xx, PTB2OE(xx), 'r')
plot([1; 1]*splineKnots(:)', ylim'*ones(1,numel(splineKnots)), 'k--')

subplot(1,2,2)
% ptbTs = sixletsPTBts(goodPTBindex,:);
% oeTs  = sixletsOEts(hind(goodPTBindex),:);

d = oeTs - PTB2OE(ptbTs);

plot(oeTs, d*1e3, '.')
ylim([-1 1])
ylabel('error (ms)')
xlabel('oe timestamp')
OE2PTBfit = [];
maxreconstructionerror = max(abs(d));


% 
% %%
% nPTBtrials = size(sixletsPTB,1);
% 
% oeMatch = nan(nPTBtrials,1);
% for iTrial = 1:nPTBtrials
%     
%     ptbSixlet = sixletsPTB(iTrial,:);
%     
%     % check for matches
%     matches = find(all(bsxfun(@eq, sixletsOE, ptbSixlet),2));
%     
%     switch numel(matches)
%         case 1
%             oeMatch(iTrial) = matches;
%         case 2
%             
%             % check if a trial has already been matched -- they should be
%             % sequential because they were sorted above meaning we can
%             % discount this match
%             [intersection, ~, ib] = intersect(oeMatch, matches);
%             if verbose
%                 fprintf('detected 2 matches\n')
%                 if numel(intersection) > 0
%                     fprintf('removing matches that were already counted\n')
%                 else
%                     fprintf('Using the first match\n')
%                 end
%             end
%                 
% %             matches(ib) = [];
% %             oeMatch(iTrial) = matches(1);
%             
%             
%         otherwise
%     end
%             
%     
% end
% 
% goodPTBindex = ~isnan(oeMatch);
% 
% if verbose
%     fprintf('Found %d/%d trials matched\n', sum(goodPTBindex), nPTBtrials)
% end
% 
% oeMatch(~goodPTBindex) = [];
% 
% clf
% plot(datenum(sixletsPTB(goodPTBindex,:)),datenum(sixletsOE(oeMatch,:)), '.')
% 
% %%
% sixletsOEts = sixletsOEts(oeMatch,:); %double([timestamps(strobeSet(oeMatch)) timestamps(strobeSet(oeMatch+1)) timestamps(strobeSet(hind+2)) timestamps(strobeSet(hind+3)) timestamps(strobeSet(hind+4)) timestamps(strobeSet(hind+5))]);
% sixletsPTBts=sixletsPTBts(goodPTBindex,:);
% 
% clf
% plot(sixletsPTBts, sixletsOEts - sixletsPTBts, '.')
% %%
% % sixletsDPts=sixletsDPts(goodPTBindex,:);
% 
% 
% OE2PTBfit=[sixletsOEts(:) ones(numel(sixletsOEts),1)]\sixletsPTBts(:);
% OE2PTB=@(x) x*OE2PTBfit(1) + OE2PTBfit(2);
% PTB2OE=@(x) (x - OE2PTBfit(2))/OE2PTBfit(1);
% 
% % pp2 = splinefit(sixletsPTBts(:), sixletsOEts(:), linspace(min(sixletsPTBts(:)), max(sixletsPTBts(:)), 100));
% PTB2OEfit = splinefit(sixletsPTBts(:), sixletsOEts(:), 10, 2);
% PTB2OE = @(x) ppval(PTB2OEfit, x);
% 
% % get a reconstruction estimate
% maxreconstructionerror = max(abs(PTB2OE(sixletsPTBts(:))-sixletsOEts(:)));
% 
% %%
% if verbose
%     fprintf('Max reconstruction error = %02.3f ms\n', maxreconstructionerror*1e3)
%     figure(1); clf
%     plot(PTB2OE(sixletsPTBts(:)), '.', 'MarkerSize', 10); hold on
%     plot(sixletsOEts(:), '.');
%     xlabel('stobe #')
%     ylabel('OE time')
%     legend('converted PTB', 'OE')
%     
%     sixletsPTBts_ = sixletsPTBts_';
%     sixletsPTBts_ = sixletsPTBts_(:);
%     
% %     sixletsOEts_ = sixletsOEts_';
% %     sixletsOEts_ = sixletsOEts_(:);
%     
%    %% debugging
%     d = abs(filter(ones(6, 1), 1, diff(sixletsOEts_(:,1))));
%     clf; plot(d); ylim([0 1])
%     sixletsOEts_b = sixletsOEts_(d<0.01);
%     d = abs(filter(ones(5, 1), 1, diff(sixletsPTBts_(:))));
%     sixletsPTBts_b = sixletsPTBts_(d<0.01);
% %     
%     t1 = sixletsPTBts_b(:);
%     t2 = sixletsOEts_b(:);
% % %     PTB2OE = @(x) x - b;
% % %     
%     t1 = PTB2OE(t1);
% %     
% % %     f = @(w) sum(abs(Y - X*w));
% % %     what = fmincon(f, OE2PTBfit)
% % 
%     figure(2); clf
%     plot( PTB2OE([t1 t1])', [0 1]'*ones(1,numel(t1)), 'r'); hold on
%     plot( [t2 t2]', [0 1]'*ones(1,numel(t2)), '--b');
%     
%     nPtb = numel(sixletsPTBts_b);
%     OeMatch = nan(nPtb,1);
%     err = nan(nPtb,1);
%     for i = 1:nPtb
%         
%         [err(i), OeMatch(i)] = min((t1(i) - t2).^2);
%         
%     end
%     
%     tol = 1e-6;
%     goodIndex = err < tol;
%     t1 = sixletsPTBts_b(:);
%     t2 = sixletsOEts_b(:);
%     
%     t1 = t1(goodIndex);
%     t2 = t2(OeMatch(goodIndex));
% %     
% %     X = [t2 ones(numel(t1),1)];
% %     Y = t1;
% %     
% %     OE2PTBfit = (X'*X)\(X'*Y(:));
% %     OE2PTB=@(x) x*OE2PTBfit(1) + OE2PTBfit(2);
% %     PTB2OE=@(x) (x - OE2PTBfit(2))/OE2PTBfit(1);
% %     
% %     w0 = OE2PTBfit;
% %     
% %     fmap = @(t1, w) ( (t1 -w(2))./w(1)  + w(5)*sin( w(3)*(t1 - w(2))./w(1) +w(4)));
% %     f = @(w) sum((t2 - fmap(t1, w) ).^2);
% %     w0 = [w0(1) w0(2) 0 0 1]';
% %     f(w0)
% %     
% %     what = fminsearch(f, w0);
% %     f(what)-f(w0)
% %     %%
%     clf
% %     t1 = sixletsPTBts(:); 
% %     t2 = sixletsOEts(:);
%     plot(t1, PTB2OE(t1) - t2, '.'); hold on
%     
%     PTB2OEfit = splinefit(t1, t2, 10);
%     
%     PTB2OE = @(x) ppval(PTB2OEfit, x);
% 
%     
%     %%
% 	bad = find(abs(zscore(PTB2OE(t1) - t2)) > 5);
%     
%     t1(bad) = [];
%     t2(bad) = [];
%     
%     PTB2OEfit = splinefit(t1,t2,10);  % 11 breaks, 10 pieces
%     pp3 = splinefit(t2, t1, 10);
%     PTB2OE = @(x) ppval(PTB2OEfit, x);
%     OE2PTB = @(x) ppval(pp3, x);
%     
%     
%     OE2PTBfit=pp3;
%     PTB2OEfit = PTB2OEfit;
%     % X = [sixletsOEts(:) ones(numel(sixletsOEts),1)];
%     % Y = sixletsPTBts(:);
% % OE2PTBfit = (X'*X)\(X'*Y(:));
% %     OE2PTB=@(x) x*OE2PTBfit(1) + OE2PTBfit(2);
% %     PTB2OE=@(x) (x - OE2PTBfit(2))/OE2PTBfit(1);
% 
%     clf; plot(PTB2OE(t1), t2, '.')
% % get a reconstruction estimate
%     maxreconstructionerror = max(abs(PTB2OE(sixletsPTBts(:))-sixletsOEts(:)));
%     fprintf('Max reconstruction error = %02.3f ms\n', maxreconstructionerror*1e3)
% %     plot(t1, fmap(t1, what)-t2, '-')
% %     plot(t1, sin(550 - 500)/5e3)
%     
%     
%     clf
%     plot(t1, PTB2OE(t1) - t2, '.'); hold on
%     drawnow
% end
