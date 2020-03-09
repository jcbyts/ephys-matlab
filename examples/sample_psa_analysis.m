

%% Load up the data

% first get the path to where the data live
fp = fullfile(getpref('EPHYS', 'SERVER_DATA'), 'pds_backups');

% find all PDS file that start with Ellie
fl = dir(fullfile(fp, 'Ellie*.PDS'));

% search for runGaborTargetSelection (this might've been possible in the
% previous line by using more *
idx = arrayfun(@(x) any(strfind(x.name, 'runGaborTargetSelection')), fl);

fl = fl(idx);

% load a subset of the analyzable files
filesToLoad = numel(fl)-4; %numel(fl)-4:1:numel(fl);
nFiles = numel(filesToLoad);
PDS = cell(nFiles,1);
for iFile = 1:nFiles

    thisFile = filesToLoad(iFile);
    tmp = load(fullfile(fp, fl(thisFile).name), '-mat');
    PDS{iFile} = tmp.PDS;
    % since we don't have open ephys time, just align all time to the start
    % of the first session
    if iFile == 1
        PDS{iFile}.PTB2OE = @(x) x-tmp.PDS.initialParametersMerged.session.experimentStart;
    else
        PDS{iFile}.PTB2OE = @(x) x-PDS{1}.initialParametersMerged.session.experimentStart;
    end
end
    

%% concatenate and convert to psaForage object
D = session.psaForage(PDS);

%% analyze proportion chosen as a function of trial-lag-from-switch-time

goodTrials = ~isnan([D.trial.targChosen])';
nGood = sum(goodTrials);
nTrials = numel(goodTrials);
fprintf('%d / %d completed Trials (%02.0f%% completion rate)\n', nGood, nTrials, nGood/nTrials*100)

targRewarded = nan(nTrials, 1);
targChosen   = nan(nTrials, 1);

for iTrial = 1:nTrials
    
    targChosen(iTrial) = D.trial(iTrial).targChosen;
    
    % TODO: this needs to be more sophisticated if we use probabilistic reward
    % rates -- IMPORTANT!!
    targRewarded(iTrial) = find(D.trial(iTrial).isRewarded);
    
end

% proportion of times she chose target 1 given that target 1 was rewarded
ix = targRewarded == 1;
targ1Rewarded = targChosen(ix & goodTrials)==1;
pTarg1 = mean(targ1Rewarded);
pTarg1SE = std(targ1Rewarded)/sqrt(numel(targ1Rewarded));

fprintf('Chose target 1 %02.2f%% +- %02.2f%% of the time target 1 was rewarded\n', pTarg1, pTarg1SE)

%% find the switches
figure(1); clf
plot(targRewarded); hold on

switchTrials = find(diff(targRewarded)~=0)+1;
plot([1; 1]*switchTrials', ylim'*ones(1, numel(switchTrials)), 'r')
xlabel('Trial #')
ylabel('Targ Rewarded')
legend({'targ rewarded', 'switch time'})


%% loop over switches and calculate lagged probability of choosing a target

nSwitches = numel(switchTrials)-1; % ignore the last switch so we don't excede the number of trials we have
trialLags = -1:7;
laggedChoice = nan(nSwitches, numel(trialLags));
laggedReward = nan(nSwitches, numel(trialLags));

% TODO: we don't want to count trials across sessions. Because we
% concatenated them, there's an error here that may or may not detect a
% switch when a new session started... need to figure out what to do here
for iSwitch = 1:nSwitches

    thisSwitch  = switchTrials(iSwitch);
    laggedIndex = thisSwitch + trialLags;
    
    laggedChoice(iSwitch, :) = targChosen(laggedIndex);
    laggedReward(iSwitch, :) = targRewarded(laggedIndex);
    
    % find if another switch happened
    nextSwitch=find(diff(laggedReward(iSwitch, trialLags >= 0))~=0, 1, 'first');
    if isempty(nextSwitch)
        continue
    else
        laggedChoice(iSwitch,nextSwitch:end) = nan;
        laggedReward(iSwitch,nextSwitch:end) = nan;
    end
    
end

rewardState = laggedReward(:,trialLags==0);

% convert the 1s and 2s to 0 and 1s appropriate while keeping nans
choice1reward1 = abs(laggedChoice(rewardState==1,:)-2);
choice1reward2 = abs(laggedChoice(rewardState==2,:)-2);

figure(2); clf
propCho1R1 = nanmean(choice1reward1);
propCho1R2 = nanmean(choice1reward2);
propCho1R1SE = nanstd(choice1reward1)./sqrt(sum(~isnan(choice1reward1)));
propCho1R2SE = nanstd(choice1reward2)./sqrt(sum(~isnan(choice1reward2)));

errorbar(trialLags, propCho1R1, propCho1R1SE); hold on
errorbar(trialLags, propCho1R2, propCho1R2SE)


legend({'p(C=1|R=1)', 'p(C=1|R=2)'})
ylabel('p(C==1)')
xlabel('Trial from switch')

 
 
% %%
% targRewarded = cell2mat(arrayfun(@(x) x.rewardRate, D.trial(:), 'uni', false));
% targChosen   = arrayfun(@(x) x.targChosen, D.trial(:));
% rewardSwitches = find([D.trial.switchReward]);
% 
% trialLags = -1:8; % check 8 trials past switch
% 
% nSwitches = numel(rewardSwitches);
% nLags = numel(trialLags);
% completed = nan(nSwitches, nLags);
% choice = nan(nSwitches, nLags);
% rewardState = nan(nSwitches, nLags);
% 
% for iSwitch = 1:nSwitches
%     ix = rewardSwitches(iSwitch) + trialLags;
%     iix = 1:numel(ix);
%     % make sure index is valid
%     iix(ix<1 | ix>nTrials) = []; 
%     ix(ix<1 | ix>nTrials) = [];
%     
%     if iSwitch < nSwitches
%         nextSwitch = rewardSwitches(iSwitch + 1);
%         
%         iix(ix > nextSwitch) = [];
%         ix(ix > nextSwitch) = [];
%     end
%     
%     rewardState(iSwitch,iix) = targRewarded(ix,1);
%     completed(iSwitch, iix) = ~isnan(targChosen(ix));
%     choice(iSwitch,iix) = targChosen(ix);
%     
% end
% %%
% 
% figure(1); clf
% subplot(1,2,1)
% completionRate = nanmean(completed);
% completionSE   = nanstd(completed)./sqrt(sum(~isnan(completed)));
% errorbar(trialLags, completionRate, completionSE)
% xlabel('Trial Lag from switch')
% ylabel('Completion Rate')
% 
% subplot(1,2,2)
% choiceRateGivenRight = nanmean(choice(rewardState(:,1)==0,:)==1);
% choiceRateGivenLeft = nanmean(choice(rewardState(:,1)==1,:)==1);
% 
% choiceSEGivenRight = nanstd(choice(rewardState(:,1)==0,:)==1)./sqrt(sum(~isnan(choice(rewardState(:,1)==0,:))));
% choiceSEGivenLeft = nanstd(choice(rewardState(:,1)==1,:)==1)./sqrt(sum(~isnan(choice(rewardState(:,1)==1,:))));
% 
% errorbar(trialLags, choiceRateGivenLeft, choiceSEGivenLeft); hold on
% errorbar(trialLags, choiceRateGivenRight, choiceSEGivenRight);
% 
% figure(2); clf
% 
% plot(targRewarded(:,1));  hold on
% plot(targChosen-1 == targRewarded(:,1), '.');
% plot(targChosen-1 ~= targRewarded(:,1), '.');

%%
ev = nan(D.numTrials,1);
choice = nan(D.numTrials, 1);

for i = 1:D.numTrials
    ev(i) = D.trial(i).choiceTime;
    choice(i) = D.trial(i).targChosen;
end


figure(1); clf
[m, s, bc]=pdsa.eventPsth(sp{1}.st, ev(choice==1), [-1, 1.5], .05);
errorbar(bc, m, s)
hold on

[m, s, bc]=pdsa.eventPsth(sp{1}.st, ev(choice==2), [-1, 1.5], .05);
errorbar(bc, m, s)

%%
