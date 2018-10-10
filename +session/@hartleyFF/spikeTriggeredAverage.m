function sta = spikeTriggeredAverage(h, spikeTimes)

%            bins = (1:size(h.design.Xd))*h.display.ifi + h.design.rowTimes(1);
bins = h.design.rowTimes;
% bin spikes at the frame rate
y = histc(spikeTimes, bins);

%             figure(1); clf
%             plot(y)
y(diff(bins)>h.display.ifi) = 0;

y = y(:);

sta.w = (h.design.Xd'*y);
sta.xy = sta.w;
sta.yy = y'*y;
sta.ny = numel(y);
sta.w = (h.design.XX + 10e2 * speye(size(h.design.XX,2)))\sta.w;
sta.fullRF  = full(reshape(sta.w(1:end-1), [h.design.nkTime h.design.nkx*h.design.nky]));

[u,~,v] = svd(sta.fullRF);
u(:,1) = u(:,1) - mean(u(1:5,1));
[~, im] = max(abs(u(:,1)));
sflip = sign(u(im,1));

sta.kxs = h.design.kxs;
sta.kys = h.design.kys;
sta.RF = reshape(sflip*v(:,1), h.design.nkx, h.design.nky);
sta.time = (1:h.design.nkTime)*h.display.ifi;
sta.RFtime = flipud(sflip*u(:,1));

end