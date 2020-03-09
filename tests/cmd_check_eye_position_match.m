

%% Load a session
% rerun getExperimentsAnd to refresh the meta data
meta = io.getExperimentsAnd(); % get all experiments meta data

%thisSession = meta(148,:); % again grab the session you are working with
thisSession = meta(end,:);
disp(thisSession)

fprintf('Loading Behavioral data from the server...')
tic
PDS = io.getPds(thisSession);
fprintf(' [%02.2f]\n', toc)


% load Eye position
[data, timestamps, elInfo] = io.getEdf(thisSession, PDS, false);

% --- remove bad samples
ix = any(data(1:2,:) == elInfo.bitDeg(2));
data(1,ix) = nan;
data(2,ix) = nan;

%% plot
figure(1); clf
subplot(2,1,1)
plot(timestamps, data(1,:), 'k'); hold on

subplot(2,1,2)
plot(timestamps, data(2,:), 'k'); hold on


nPds = numel(PDS);
for kPds = 1:nPds
    
    flipTimes = cell2mat(cellfun(@(x) PDS{kPds}.PTB2OE(x.timing.flipTimes(1,:))', PDS{kPds}.data, 'uni', false)');
    
    eyeAtFrame = cell2mat(cellfun(@(x) x.behavior.eyeAtFrame', PDS{kPds}.data, 'uni', false)');
    
    % center
    eyeAtFrame = bsxfun(@minus, eyeAtFrame, PDS{kPds}.initialParametersMerged.display.ctr(1:2));
    
    % flip y
    eyeAtFrame(:,2) = -eyeAtFrame(:,2);
    
    % pixels to degrees
    eyeAtFrame = pds.px2deg(eyeAtFrame', PDS{kPds}.initialParametersMerged.display.viewdist, PDS{kPds}.initialParametersMerged.display.px2w)';
    
    subplot(2,1,1)
    plot(flipTimes, eyeAtFrame(:,1), '.'); hold on
    
    subplot(2,1,2)
    plot(flipTimes, eyeAtFrame(:,2), '.'); hold on
    
end
    