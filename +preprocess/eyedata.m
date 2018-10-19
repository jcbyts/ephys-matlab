function [hObj, sac, blnk, purs] = eyedata(timestamps, eyeDeg, elInfo)
% [hObj, saccades, blinks, pursuit] = eyedata(timestamps, eyeDeg, elInfo)
verbose = false;
%--- curate data
ix = any(eyeDeg(1:2,:)==elInfo.bitDeg(2));
eyeDeg(1,ix) = nan;
eyeDeg(2,ix) = nan;

nans = find(any(isnan(eyeDeg),1));
dn = diff(nans);
segment_starts = nans(dn > 1) + 1;
segment_stops   = segment_starts + dn(dn > 1) - 2;

if verbose
    figure(1); clf
    plot(timestamps, eyeDeg(1,:)); hold on
    plot(timestamps(segment_starts), eyeDeg(1,segment_starts), '.')
    plot(timestamps(segment_stops), eyeDeg(1,segment_stops), '.')
end

%%
numSegments = numel(segment_starts);
fprintf('Found %d segments in the eye position recording\n', numSegments)
fprintf('Looping over now\n')

blinks   = repmat(struct('tstart', [], 'tend', [], 'idx', [], 'delta', []), numSegments, 1);
saccades = repmat(struct('tstart', [], 'tend', [], 'idx', [], 'delta', []), numSegments, 1);
pursuit  = repmat(struct('tstart', [], 'tend', [], 'idx', [], 'delta', []), numSegments, 1);

clear hObj
for iSegment = 1:numSegments
    ix = segment_starts(iSegment):segment_stops(iSegment);
    
    hObj(iSegment) = session.eyepos(timestamps(ix), eyeDeg(1,ix)', eyeDeg(2,ix)', sqrt(eyeDeg(3,ix)/pi)', sqrt(eyeDeg(3,ix)/pi)');
    
    nsamples = nan;
    try
        fprintf('Resampling eye position at 1kHz\n')
        tic;
        hObj(iSegment) = hObj(iSegment).resample('debug', false, 'method', 'spline');
        fprintf('Done [%02.2f]\n', toc)
        
        nsamples = sum(arrayfun(@(x) numel(x.tsample), hObj(1:iSegment-1)));
        
        fprintf('Segment %d: Detecting Blinks\n', iSegment); tic;
        [blinks(iSegment).tstart,blinks(iSegment).tend,blinks(iSegment).idx,blinks(iSegment).delta]         = hObj(iSegment).findBlinks('debug', false);
        fprintf('Done [%02.2f]\n', toc)
        
        fprintf('Segment %d: Detecting Saccades\n', iSegment); tic;
        [saccades(iSegment).tstart,saccades(iSegment).tend,saccades(iSegment).idx,saccades(iSegment).delta] = hObj(iSegment).findSaccades('debug', false);
        fprintf('Done [%02.2f]\n', toc)
        
        fprintf('Segment %d: Detecting Pursuit\n', iSegment); tic;
        
        [pursuit(iSegment).tstart,pursuit(iSegment).tend,pursuit(iSegment).idx,pursuit(iSegment).delta]     = hObj(iSegment).findPursuit('debug', false);
        fprintf('Done [%02.2f]\n', toc)
        
        blinks(iSegment).idx = blinks(iSegment).idx + nsamples;
        saccades(iSegment).idx = saccades(iSegment).idx + nsamples;
        pursuit(iSegment).idx = pursuit(iSegment).idx + nsamples;
        
        if verbose
            figure(1); clf
            subplot(2,1,1)
            plot(timestamps(ix), eyeDeg(1,ix), '--'); hold on
            plot(hObj(iSegment).tsample, hObj(iSegment).x, 'k');
            ylabel('x position')
            
            xsacs = [saccades(iSegment).tstart,saccades(iSegment).tend,saccades(iSegment).tend,saccades(iSegment).tstart];
            
            xblinks = [blinks(iSegment).tstart,blinks(iSegment).tend,blinks(iSegment).tend,blinks(iSegment).tstart];
            
            xpurs   = [pursuit(iSegment).tstart,pursuit(iSegment).tend,pursuit(iSegment).tend,pursuit(iSegment).tstart];
            yy = kron(get(gca,'YLim'),[1,1]);
            
            % plot saccades
            lh = [];
            lstr = {};
            nSacs = size(xsacs,1);
            for iSac = 1:nSacs
                fhs = fill(xsacs(iSac,:)',yy(:),zeros(1,3));
                set(fhs,'FaceAlpha',0.1,'LineStyle','none');
            end
            if nSacs >=1
                lh = [lh fhs];
                lstr = [lstr {'saccades'}];
            end
            
            % plot blinks
            nBlinks = size(xblinks,1);
            for iBlink = 1:nBlinks
                fhb = fill(xblinks(iBlink,:)',yy(:),[1 0 0]);
                set(fhb,'FaceAlpha',0.1,'LineStyle','none');
            end
            if nBlinks >=1
                lh = [lh fhb];
                lstr = [lstr {'blinks'}];
            end
            
            % plot pursuit
            nPursuit = size(xpurs,1);
            for iPursuit = 1:nPursuit
                fhp = fill(xpurs(iPursuit,:)',yy(:),[0 1 0]);
                set(fhp,'FaceAlpha',0.1,'LineStyle','none');
            end
            
            if nPursuit >= 1
                lh = [lh fhp];
                lstr = [lstr {'pursuit'}];
            end
            
            legend(lh, lstr, 'Location', 'Best')
            
            subplot(2,1,2)
            plot(timestamps(ix), eyeDeg(2,ix), '.'); hold on
            plot(hObj(iSegment).tsample, hObj(iSegment).y, 'k');
            yy = kron(get(gca,'YLim'),[1,1]);
            % plot saccades
            lh = [];
            lstr = {};
            nSacs = size(xsacs,1);
            for iSac = 1:nSacs
                fhs = fill(xsacs(iSac,:)',yy(:),zeros(1,3));
                set(fhs,'FaceAlpha',0.1,'LineStyle','none');
            end
            if nSacs >=1
                lh = [lh fhs];
                lstr = [lstr {'saccades'}];
            end
            
            % plot blinks
            nBlinks = size(xblinks,1);
            for iBlink = 1:nBlinks
                fhb = fill(xblinks(iBlink,:)',yy(:),[1 0 0]);
                set(fhb,'FaceAlpha',0.1,'LineStyle','none');
            end
            if nBlinks >=1
                lh = [lh fhb];
                lstr = [lstr {'blinks'}];
            end
            
            % plot pursuit
            nPursuit = size(xpurs,1);
            for iPursuit = 1:nPursuit
                fhp = fill(xpurs(iPursuit,:)',yy(:),[0 1 0]);
                set(fhp,'FaceAlpha',0.1,'LineStyle','none');
            end
            
            if nPursuit >= 1
                lh = [lh fhp];
                lstr = [lstr {'pursuit'}];
            end
            
            drawnow
        end
    end % try
end

%% concatentate


tsample = cell2mat(arrayfun(@(x) x.tsample, hObj(:), 'uni', 0));
x = cell2mat(arrayfun(@(x) x.x, hObj(:), 'uni', 0));
y = cell2mat(arrayfun(@(x) x.y, hObj(:), 'uni', 0));

iix = (abs(x) > 20 | abs(y) > 20);
x(iix) = nan;
y(iix) = nan;

pwdth = cell2mat(arrayfun(@(x) x.pwdth, hObj(:), 'uni', 0));
phght = cell2mat(arrayfun(@(x) x.phght, hObj(:), 'uni', 0));

plot(tsample, x)
%%
ix = arrayfun(@(x) ~isempty(x.idx), saccades);
sac_idx = cell2mat(arrayfun(@(x) x.idx, saccades(ix), 'uni', 0));


blink_start = cell2mat(arrayfun(@(x) x.tstart(:), blinks, 'uni', 0));
blink_end = cell2mat(arrayfun(@(x) x.tend(:), blinks, 'uni', 0));
sac_start = cell2mat(arrayfun(@(x) x.tstart(:), saccades, 'uni', 0));
sac_end = cell2mat(arrayfun(@(x) x.tend(:), saccades, 'uni', 0));
pursuit_start = cell2mat(arrayfun(@(x) x.tstart(:), pursuit, 'uni', 0));
pursuit_end   = cell2mat(arrayfun(@(x) x.tend(:), pursuit, 'uni', 0));

% convert times to indices (again)
nSaccades = numel(sac_start);
sac = struct('tstart', sac_start,...
    'tend', sac_end, ...
    'duration', sac_end-sac_start,...
    'size', nan(nSaccades, 1), ...
    'startXpos', nan(nSaccades, 1), ...
    'startYpos', nan(nSaccades, 1), ...
    'endXpos', nan(nSaccades, 1), ...
    'endYpos', nan(nSaccades, 1), ...
    'startIndex', nan(nSaccades, 1), ...
    'endIndex', nan(nSaccades, 1), ...
    'dx', nan(nSaccades, 1), ...
    'dy', nan(nSaccades, 1), ...
    'vel', nan(nSaccades, 1));

sac.startIndex = arrayfun(@(x) find((tsample - x)>=0, 1), sac_start);
sac.endIndex   = arrayfun(@(x) find((tsample - x)>=0, 1), sac_end);
sac.startXpos = x(sac.startIndex);
sac.startYpos = y(sac.startIndex);
sac.endXpos   = x(sac.endIndex);
sac.endYpos   = y(sac.endIndex);
sac.dx = sac.endXpos - sac.startXpos;
sac.dy = sac.endYpos - sac.startYpos;
sac.size = sqrt(sac.dx.^2 + sac.dy.^2);
sac.vel = sac.size./sac.duration;

nPursuit = numel(pursuit_start);
purs = struct('tstart', pursuit_start,...
    'tend', pursuit_end, ...
    'duration', pursuit_end-pursuit_start,...
    'startXpos', nan(nPursuit, 1), ...
    'startYpos', nan(nPursuit, 1), ...
    'endXpos', nan(nPursuit, 1), ...
    'endYpos', nan(nPursuit, 1), ...
    'startIndex', nan(nPursuit, 1), ...
    'endIndex', nan(nPursuit, 1), ...
    'dx', nan(nPursuit, 1), ...
    'dy', nan(nPursuit, 1));

purs.startIndex = arrayfun(@(x) find((tsample - x)>=0, 1), pursuit_start);
purs.endIndex   = arrayfun(@(x) find((tsample - x)>=0, 1), pursuit_end);
purs.startXpos = x(purs.startIndex);
purs.startYpos = y(purs.startIndex);
purs.endXpos   = x(purs.endIndex);
purs.endYpos   = y(purs.endIndex);
purs.dx = purs.endXpos - purs.startXpos;
purs.dy = purs.endYpos - purs.startYpos;
purs.size = sqrt(purs.dx.^2 + purs.dy.^2);

blnk = struct('tstart', blink_start,...
    'tend', blink_end, ...
    'duration', blink_end-blink_start,...
    'startIndex', arrayfun(@(x) find((tsample - x)>=0, 1), blink_start), ...
    'endIndex', arrayfun(@(x) find((tsample - x)>=0, 1), blink_end));

hObj = session.eyepos(tsample, x, y, pwdth, phght);
