function CsdBasic(lfp, eventTimes, lfpInfo, varargin)
% CsdBasic(lfp, eventTimes, lfpInfo, varargin)
% 'upSampleFactor'
% ip.addParameter('window', [-100 200])
% ip.addParameter('channelDepths', [])
% ip.addParameter('skipChannels', 1)

ip = inputParser();
ip.addParameter('upSampleFactor', 5)
ip.addParameter('window', [-100 200])
ip.addParameter('channelDepths', [])
ip.addParameter('skipChannels', 1)
ip.parse(varargin{:});

eventTimes = eventTimes(:);

% conver times to samples
ev = io.convertTimeToSamples(eventTimes, lfpInfo.sampleRate, lfpInfo.timestamps(:), lfpInfo.fragments(:));

[sta,~, xax] = pdsa.eventTriggeredAverage(lfp, ev(:), ip.Results.window);


ims = ip.Results.upSampleFactor;
if isempty(ip.Results.channelDepths)
    ch0 = -(1:32)*50;
end
chInd = 1:ip.Results.skipChannels:32;
staUp = imresize(sta(:,chInd), ims);
chans = imresize(ch0(chInd), ims);
time  = imresize(xax, ims);
time  = time(1,:);
chans = chans(1,:);
CSD = diff(staUp, [], 2)';

imagesc(time, chans, CSD-mean(CSD(:))); axis xy
colormap jet
hold on
plot(xax, bsxfun(@plus, sta, ch0), 'Color', repmat(.5, 1, 3))

xlim(ip.Results.window)
