function sta = AsdRf(h, spikeTimes)

bins = (1:size(h.design.Xd))*h.display.ifi + h.design.rowTimes(1);
% bin spikes at the frame rate
y = histc(spikeTimes, bins);

y = y(:);

sta.w = fastASD_2D_hyper_dual(h.design.Xd(:,setdiff(1:size(h.design.Xd,2), h.design.biasCol)), y-mean(y), [h.design.nkTime h.design.nkx*h.design.nky], 1);
sta.fullRF  = full(reshape(sta.w, [h.design.nkTime h.design.nkx*h.design.nky]));

[u,~,v] = svd(sta.fullRF);
u(:,1) = u(:,1) - mean(u(1:5,1));
[~, im] = max(abs(u(:,1)));
sflip = sign(u(im,1));

sta.kxs = h.design.kxs;
sta.kys = h.design.kys;
sta.RF = reshape(sflip*v(:,1), h.design.nkx, h.design.nky);
sta.time = (1:h.design.nkTime)*h.display.ifi;
sta.RFtime = sflip*u(:,1);

end