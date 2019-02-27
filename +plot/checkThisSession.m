function checkThisSession(thisSession)

PDS = io.getPds(thisSession);
spikes = io.getSpikes(thisSession);

hart = session.hartleyFF(PDS);
sp = spikes{1};
nUnits = numel(sp.cids);
sx = ceil(sqrt(nUnits));
sy = round(sqrt(nUnits));

if hart.numTrials > 10
    hart.buildDesignMatrix
    
    
    
    figure(1); clf
    
    for iUnit = 1:nUnits
        sta = hart.spikeTriggeredAverage(sp.st(sp.clu==sp.cids(iUnit)));
        
        subplot(sx, sy, iUnit)
        if sp.isiV(iUnit)<.25
            plot(sta.time, sta.RFtime, 'k')
        else
            plot(sta.time, sta.RFtime, 'b')
        end
        drawnow
    end
end

%%
csd = session.csdFlash(PDS);

if csd.numTrials > 5
ops = io.loadOps(thisSession);
figure
csd.computeCsd(ops(1), 'plot', true);
end
