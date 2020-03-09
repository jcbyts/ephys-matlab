function [headstage] = combineHeadstages(headstages)

% --- build channel map info
numHeadstages = numel(headstages);

nHeadstageChannels = zeros(numHeadstages,1);
for iHeadstage = 1:numHeadstages
    nHeadstageChannels(iHeadstage) = sum(~isnan(headstages{iHeadstage}.channelMap));
end

% build new combined headstage
headstage = hardware.headstage.headstage;

chOffset = zeros(numHeadstages,1);
chCount  = zeros(numHeadstages,1);
for iHeadstage = 1:numHeadstages
    headstage.name         = [headstage.name headstages{iHeadstage}.name ','];
    headstage.manufacturer = [headstage.manufacturer headstages{iHeadstage}.manufacturer ','];
    headstage.model        = [headstage.model headstages{iHeadstage}.model ','];
    headstage.connector    = [headstage.connector headstages{iHeadstage}.connector ','];
    headstage.filter       = [headstage.filter; headstages{iHeadstage}.filter];
    headstage.samplingRate = [headstage.samplingRate headstages{iHeadstage}.samplingRate];
    headstage.gains        = [headstage.gains headstages{iHeadstage}.gains];
    tmp = max(headstage.channelMap);
    if isempty(tmp)
        chOffset(iHeadstage) = 0;
    else
        chOffset(iHeadstage) = tmp;
    end
    headstage.channelMap = [headstage.channelMap headstages{iHeadstage}.channelMap + chOffset(iHeadstage)];
    chCount(iHeadstage)  = numel(headstages{iHeadstage}.channelMap);
end

% remove comma at end
headstage.name(end) = [];
headstage.manufacturer(end) = [];
headstage.model(end) = [];
headstage.connector(end) = [];