function removeArtifactsManual(thisSession, blockSize)


% --- Load session
[~, ops, info] = io.loadSession(thisSession);

% --- Loop over loading

nSampleTotal = sum(info.fragments);

if nargin < 2
    blockSize = 300;  % 5 minute blocks
end

if isinf(blockSize)
    blockSize = nSampleTotal;
else
    blockSize = info.sampleRate*blockSize;
end

fprintf('Loading recording in 5 minute blocks\n')
nBlocks = ceil(nSampleTotal/blockSize);

wfs  = [];
ssix = [];

thresh = 2e4;

% build filter
[b,a] = butter(3, [(300/30e3*2), (5000/30e3*2)]);

for iBlock = 1:nBlocks
    fprintf('Loading block %d/%d\n', iBlock, nBlocks)
    % blockSize = nSampleTotal;
    
    win = (iBlock - 1) * blockSize + [1 blockSize];
    win(2) = min(win(2), nSampleTotal);
    
    dataRAW = io.loadRaw(ops, win, true);
        
    
    fprintf('\tFiltering\n')
    dataRAW = dataRAW';
    dataHP = filter(b,a,dataRAW);
    dataHP = flipud(dataHP);
    dataHP = filter(b,a,dataHP);
    dataHP = flipud(dataHP);
    
    % common average reference
    dataHP = bsxfun(@minus, dataHP, mean(dataHP,2));
    
    % threshold based on the energy
    En = max(dataHP.^2,[],2); % maximum across channels
    f = figure(2); clf
    f.Position = [500 200 500 500];
    h = histogram(En);
    set(gca, 'yscale', 'log')
    xlabel('Energy')
    ylabel('Count')
    
    
    figure(1); clf
    plot(En)
    axis tight
    xlabel('Sample #')
    ylabel('Energy')
    hold on
    
    
    if iBlock==1
        title('Set high threshold for clipping artifacts')
        drawnow
        
        [x, thresh] = ginput(1);
    end
    
    plot(xlim, thresh*[1 1], 'r--')
	drawnow
    
    % detect threshold crossings
    fprintf('\tThresholding\n')
    ss = find(any(En > thresh,2));
    ss(diff(ss) < 32) = [];
    spWin = -10:32;
    
    ss((ss + min(spWin)) < 1)= []; 
    ss((ss + max(spWin)) > size(dataHP,1)) = [];
    
    fprintf('\tClipping Waveforms\n')
    wf = ephys.spikeSorting.extractWaveforms(dataHP, ss, spWin);
    X = permute(wf, [1 3 2]);
    sz = size(X);
    X = reshape(X, [], sz(3));
    
    wfs = [wfs X];
    ssix = [ssix; (iBlock - 1) * blockSize + ss];
end

% --- Run waveform selector
assert(numel(ssix)<10e4, 'Too many waveforms selected. Run again and set a higher threshold.')

app = WaveformSelector(ssix, wfs);

waitfor(app.UIFigure)

% --- Select waveforms for deletion
removeInd = setdiff(ssix, app.SS);

removeInd(diff(removeInd) < info.sampleRate) = [];

nSelected = numel(removeInd);

fprintf('%d segments set to be evluated\n', nSelected)

isArtifact = false(nSelected,1);
%%

fig = uifigure;
ax  = uiaxes('Parent', fig, 'Units', 'pixels', ...
    'Position', [50 80 400 250]);
buttonIs = uibutton(fig, 'push');
buttonIs.ButtonPushedFcn = @(btn,event) buttonIsPushed(btn,ax);
buttonIs.Position = [450 190 100 22];
buttonIs.Text = 'Yes is artifact';

buttonIsnt = uibutton(fig, 'push');
buttonIsnt.ButtonPushedFcn = @(btn,event) buttonIsntPushed(btn,ax);
buttonIsnt.Position = [450 160 100 22];
buttonIsnt.Text = 'No it isnt';

buttonAll = uibutton(fig, 'push');
buttonAll.ButtonPushedFcn = @(btn,event) buttonAllPushed(btn,ax);
buttonAll.Position = [450 130 100 22];
buttonAll.Text = 'Mark All Yes';

kInd = 1;
drawAxes   
    
        
waitfor(fig)


removeInd = removeInd(isArtifact);

% remove 1 second after each artifact
artInds = unique(bsxfun(@plus, removeInd, -10e3:info.sampleRate));
artInds(artInds < 1 | artInds > sum(info.fragments)) = [];
fname = fullfile(ops.root, 'ephys_info.mat');
einfo = load(fname);

einfo.artifacts = artInds;

fprintf('Saving new artifact indices\n')

save(fname, '-v7.3', '-struct', 'einfo');

fprintf('Done\n')

close all force

    function buttonIsntPushed(btn,ax)
        isArtifact(kInd) = false;
        drawAxes
        checkFinish
    end
    
    function  buttonIsPushed(btn,ax)
        isArtifact(kInd) = true;
        drawAxes
        checkFinish
    end

    function  buttonAllPushed(btn,ax)
        for kInd = 1:nSelected
            isArtifact(kInd) = true;
        end
        kInd = kInd + 1;
        checkFinish
    end

    function drawAxes()
        loadWin = [-info.sampleRate/2 info.sampleRate];
        win_ = loadWin + removeInd(kInd);
        win_ = max(win_, 0);
        win_ = min(win_, sum(info.fragments));
        % --- load on the whole dataset right now
        [dataRAW, ~] = io.loadRaw(ops, win_, true);
    
        dataRAW = dataRAW';
        
        ix = win_(1):win_(2);
        ix = ix - removeInd(kInd);
        ix = ix/30e3;
        plot(ax, ix, dataRAW, 'k')
        title('Is this an artifact?')
        drawnow
    end

    function checkFinish
        kInd = kInd + 1;
        if kInd > nSelected
            close(fig)
        end
        
    end

end





               


% % Create the function for the ButtonPushedFcn callback
% function out = buttonIsPushed(btn,ax)
%         x = linspace(0,2*pi,100);
%         y = sin(x);
%         plot(ax,x,y)
% 
% 
% function buttonIsPushed(btn,ax)
%         x = linspace(0,2*pi,100);
%         y = sin(x);
%         plot(ax,x,y)













% %%
% nCh = size(dataRAW,2);
% dataNEW = dataRAW;
% 
% % TODO: remove from each block -- then resave
% blockStart = ((1:nBlocks) - 1) * blockSize + 1;
% 
% for iSeg = 1:nSegments
%     % clip out a whole second
%     widx = (-16e3:32e3) + removeInd(iSeg); 
%     repix = widx(widx > 0 & widx <= n);
%     for i = 1:nCh
%         dataNEW(repix, i) = nan;
%     end
% end
% 
% for i = 1:nCh
%     fprintf('channel %d\n', i)
%     dataNEW(:,i) = repnan(dataNEW(:,i), 'spline');
% end
%     
% %     figure(1); clf
% %     plotidx = (-100e3:100e3) + ixRemove(iSeg); 
% %     plotIx = plotidx(plotidx > 0 & plotidx <= n);
% %     plot(plotIx, bsxfun(@plus, dataRAW(plotIx,:), (1:nCh)*200), 'k'); hold on
% %     plot(plotIx, bsxfun(@plus, dataNEW(plotIx,:), (1:nCh)*200), 'b');
% %    	drawnow
% % end
% i = 0;
% %%
% i = i + 1;
% figure(1); clf
% n = size(dataRAW,1);
% ix = 1:n;
% times = io.convertSamplesToTime(ix, 32e3, info.timestamps, info.fragments);
% plot(times, dataRAW(:, i)); hold on
% plot(times, dataNEW(:, i));
% plot([1; 1]*ts', ylim'*ones(1, numel(ts)), 'k')
% 
% plot([1; 1]*[nats.trial.start], ylim'*ones(1, numel(nats.trial)), 'r')
% plot([1; 1]*[hart.trial.start], ylim'*ones(1, numel(hart.trial)), 'g')
% plot([1; 1]*[csd.trial.start], ylim'*ones(1, numel(csd.trial)), 'b')
% 
% %%
% PDS = io.getPds(thisSession);
% 
% ts = cell2mat(cellfun(@(x) cellfun(@(y) x.PTB2OE(y.timing.flipTimes(1)), x.data)', PDS, 'uniformOutput', false));
% 
% nats = session.natImgBackground(PDS);
% hart = session.hartleyFF(PDS);
% csd  = session.csdFlash(PDS);
% 
% %%
% 
% nlEn = abs(dataRAW(2:end-1,:).* (dataRAW(2:end-1,:)-dataRAW(1:end-2,:)) .*dataRAW(3:end,:));
% en = dataRAW(2:end-1,:).^2;
% 
% ch = 0;
% %%
% ch = ch + 1;
% clf
% 
% plot(zscore(nlEn(:,ch))+15)
% 
% hold on
% plot(zscore(en(:,ch))+4)
% 
% plot(zscore(dataRAW(:,ch)));
% 
% %%
% clf
% ch = ch+1;
% plot(en(:,ch), nlEn(:,ch), '.')
% hold on
% mu = mean([en(:,ch), nlEn(:,ch)]);
% C = cov([en(:,ch), nlEn(:,ch)]);
% plotellipse(mu, C, 40, 'Linewidth', 4)


