function RF = hackyGetHarlteyFunction(PDS, sp, unitIds)

stim = 'hartley';

hasStim = io.findPDScontainingStimModule(PDS, stim);

hartleyTrial = struct();
trialNum = 0;

for i = find(hasStim(:)')
    
%     if any(cellfun(@(x) any(strcmp(fieldnames(x), stim)), PDS{i}.conditions))
%         keyboard
%     end
        

    trialIx = cellfun(@(x) isfield(x, stim), PDS{i}.data); 
    
    stimTrials = find(trialIx);
    
    if isempty(stimTrials)
        continue
    end
        
    kxs=PDS{i}.data{stimTrials(1)}.hartley.kxs;
    kys=PDS{i}.data{stimTrials(1)}.hartley.kys;
    
    % check that the stimuli were the same grid the whole time
%     assert(size(unique(cell2mat(cellfun(@(x) x.(stim).kxs, PDS{i}.data(stimTrials), 'UniformOutput', false)'), 'rows'), 1) ==1, ...
%     'kx grid changed!')
% 
%     assert(size(unique(cell2mat(cellfun(@(x) x.(stim).kys, PDS{i}.data(stimTrials), 'UniformOutput', false)'), 'rows'), 1) ==1, ...
%     'ky grid changed!')

    
    for j = 1:numel(stimTrials)
        thisTrial = stimTrials(j);
    
        kTrial = trialNum + j;
        
        hartleyTrial(kTrial).frameTimes = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,1:end-1));
        hartleyTrial(kTrial).start      = hartleyTrial(kTrial).frameTimes(1);
        hartleyTrial(kTrial).duration   = PDS{i}.PTB2OE(PDS{i}.data{thisTrial}.timing.flipTimes(1,end-1)) - hartleyTrial(kTrial).start;
    
        if isfield(PDS{i}.conditions{thisTrial}, stim)
            if isfield(PDS{i}.conditions{thisTrial}.(stim), 'setupRNG')
                if strcmp(PDS{i}.conditions{thisTrial}.(stim).setupRNG, 'frozenSequence')
                    hartleyTrial(kTrial).frozenSequence = true;
                    hartleyTrial(kTrial).frozenSequenceLength = PDS{i}.conditions{thisTrial}.(stim).sequenceLength;
                else
                    hartleyTrial(kTrial).frozenSequence = false;
                    hartleyTrial(kTrial).frozenSequenceLength = nan;
                end
            
            else
                hartleyTrial(kTrial).frozenSequence = false;
                hartleyTrial(kTrial).frozenSequenceLength = nan;
            end
            
        else
                hartleyTrial(kTrial).frozenSequence = false;
                hartleyTrial(kTrial).frozenSequenceLength = nan;
        end
        
        hartleyTrial(kTrial).kx         = PDS{i}.data{thisTrial}.(stim).kx;
        hartleyTrial(kTrial).ky         = PDS{i}.data{thisTrial}.(stim).ky;
        hartleyTrial(kTrial).on         = PDS{i}.data{thisTrial}.(stim).on;
        
        eyepos = io.getEyePosition(PDS{i}, thisTrial);
        hartleyTrial(kTrial).eyeSampleTime = eyepos(:,1);
        hartleyTrial(kTrial).eyeXPx        = eyepos(:,2);
        hartleyTrial(kTrial).eyeYPx        = eyepos(:,3);
        hartleyTrial(kTrial).pupilArea     = eyepos(:,4);
    end
    
    trialNum = kTrial;
    
end

%%
flipTimes = cell2mat(arrayfun(@(x) x.frameTimes, hartleyTrial, 'UniformOutput', false))';
kx        = cell2mat(arrayfun(@(x) x.kx', hartleyTrial, 'UniformOutput', false))';
ky        = cell2mat(arrayfun(@(x) x.ky', hartleyTrial, 'UniformOutput', false))';
on        = ~(isnan(kx) | isnan(ky));

kxs = unique(kx(on));
kys = unique(ky(on));

x=arrayfun(@(x) find(x==kxs), kx(on));
y=arrayfun(@(x) find(x==kys), ky(on));

ind=sub2ind([numel(kys) numel(kxs)], y, x);

%% Build design matrix

t0 = flipTimes(1);
binsize = PDS{1}.initialParametersMerged.display.ifi;
binfun = @(t) (t==0) + ceil(t/binsize);

stimOnset = binfun(flipTimes(on)-t0);
% stimulus
X = sparse(stimOnset, ind, ones(size(ind,1),1), max(stimOnset), max(ind));

ntk=30;
Xd=rfmap.makeStimRowsSparse(X, ntk);
Xd = [Xd ones(size(X,1),1)];

%% analyze RFs binned spikes
% figure(10); clf

if ~exist('unitIds', 'var')
    [~, id] = sort(sp.clusterDepths);
    unitIds = sp.cids(id);
end
nUnits = numel(unitIds);


% ax = pdsa.tight_subplot(nUnits,2,.001, .1, .1);


clear RF
for kUnit = 1:nUnits
    
    
    st = sp.st(sp.clu==unitIds(kUnit)) - t0;
    
    ss = binfun(st(st>0 & st < flipTimes(end)-t0));
    ss(ss > size(X,1)) = [];
    
    y = sparse(ss, ones(numel(ss), 1), ones(numel(ss), 1), size(X,1), 1);
    
    
    % STA
    % sta =
    ttsta = (Xd'*Xd + 10e2*speye(size(Xd,2)))\(Xd'*y);
    RF2{kUnit} =  ttsta(1:end-1);
    sta =reshape(RF2{kUnit}, ntk, []);
 
%     RF{kUnit} = fastASD(Xd(:,1:end-1), y-mean(y), [ntk max(ind)], .1);
%     sta =reshape(RF{kUnit}, ntk, []);
%     
    
%     set(gcf, 'currentaxes', ax((kUnit-1)*2 + 1))
%     sta =reshape(sta(1:end-1), ntk, []);
%     
    
    [u,s,v] = svd(full(sta));
    u(:,1) = u(:,1) - mean(u(1:5,1));
    [~, im] = max(abs(u(:,1)));
    sflip = sign(u(im,1));
%     sflip = sign(sum(v(:,1)));
    
    
    spatialRF{kUnit} = reshape(sflip*v(:,1), numel(kxs), numel(kys));
%     imagesc(kxs, kys, spatialRF{kUnit})
%     colormap(gray.^2)
    % subplot(2,32,kUnit +32)
%     axis off
    
%     set(gcf, 'currentaxes', ax((kUnit-1)*2 + 2))
    
%     plot((1:ntk)*binsize,sflip*u(:,1), 'k-')
    
    RF(kUnit).time = -(1:ntk)*binsize;
    RF(kUnit).temporalRF = flipud(sflip*u(:,1));
    RF(kUnit).kxs = kxs;
    RF(kUnit).kys = kys;
    RF(kUnit).spatialRF = spatialRF{kUnit};
    
%     axis off
%     axis tight
%     drawnow
end

% set(gcf, 'renderer', 'painters')
% % set(ax(k), 'XTick', -1:.1:.5, 'XTickLabel', -1:.1:.5, 'TickDir', 'out')
% set(gcf, 'PaperSize', [1 5], 'PaperPosition', [0 0 1 5])
% saveas(gcf, 'hartley_RFs', 'epsc')