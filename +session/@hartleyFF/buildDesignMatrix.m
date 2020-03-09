function buildDesignMatrix(h, varargin)

ip = inputParser();
ip.addParameter('trialIdx', 1:h.numTrials)
ip.addParameter('nTimeLags', 30)
ip.parse(varargin{:})


flipTimes = cell2mat(arrayfun(@(x) x.frameTimes(:), h.trial(ip.Results.trialIdx), 'UniformOutput', false))';
kx        = cell2mat(arrayfun(@(x) x.kx(:), h.trial(ip.Results.trialIdx), 'UniformOutput', false))';
ky        = cell2mat(arrayfun(@(x) x.ky(:), h.trial(ip.Results.trialIdx), 'UniformOutput', false))';
on        = ~(isnan(kx) | isnan(ky));

% --- find total stimulus space
kxs = unique(kx(on));
kys = unique(ky(on));

nTotalFrames = sum(arrayfun(@(x) numel(x.frameTimes), h.trial));

h.design.trialIdx = ip.Results.trialIdx;
h.design.nkTime   = ip.Results.nTimeLags;
h.design.nkx      = numel(kxs);
h.design.nky      = numel(kys);
h.design.kxs      = kxs;
h.design.kys      = kys;


iFrame = 0;

sz = [h.design.nkx h.design.nky];
X = zeros(nTotalFrames, prod(sz));



nTrials = numel(h.trial);
for kTrial = 1:nTrials
    fprintf('%d / %d \n', kTrial, nTrials)
    frameIdx = iFrame + (1:numel(h.trial(kTrial).frameTimes));
    
    nFrames = numel(frameIdx);
    
    xtmp = zeros(nFrames, prod(sz));
    
    stimOn = find(h.trial(kTrial).on(1:nFrames));
    [kxj, ~] = find( bsxfun(@eq, h.trial(kTrial).kx(stimOn), h.design.kxs(:)')');
    tmp = bsxfun(@eq, h.trial(kTrial).ky(stimOn), h.design.kys(:)');
    [kyj, ~] = find( tmp' );
    
    kxy = sub2ind(sz, kxj, kyj);
    
    ind = sub2ind([nFrames prod(sz)], stimOn, kxy);
    
    xtmp(ind) = 1;
    
    X(frameIdx,:) = xtmp;
    
    
    iFrame = iFrame + nFrames;
    
end

tic
fprintf('Unwrapping convolution for STA, regression... \t')
Xd  = rfmap.makeStimRowsDense(X, h.design.nkTime);
fprintf('[%02.0fs]\n', toc)

h.design.biasCol = size(Xd,2)+1;

h.design.rowTimes = flipTimes(:);

h.design.Xd  = [Xd ones(size(X,1),1)];

% sample covariance matrix (for regression)
h.design.XX = h.design.Xd'*h.design.Xd;

fprintf('Done\n')
fprintf('%02.0f total seconds of stimulus\n', size(h.design.Xd,1)*h.display.ifi)
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%             x=arrayfun(@(x) find(x==kxs), kx(on));
%             y=arrayfun(@(x) find(x==kys), ky(on));
%
%             ind=sub2ind([numel(kys) numel(kxs)], y, x);
%
%
%             fprintf('Building the Design Matrix for %d trials\n', numel(ip.Results.trialIdx))
%
%             t0 = flipTimes(1);
%             binsize = h.display.ifi;
%
%             temporalBinningFunction = @(t) (t==0) + ceil(t/binsize);
%
%             stimOnset = temporalBinningFunction(flipTimes(on)-t0);
%
%             X = sparse(stimOnset, ind, ones(size(ind,1),1), max(stimOnset), max(ind));
%
%             ntk = ip.Results.nTimeLags;
%             Xd  = rfmap.makeStimRowsSparse(X, ntk);
%             h.design.biasCol = size(Xd,2)+1;
%
%             h.design.rowTimes = flipTimes(:);
%             h.design.binfun   = temporalBinningFunction;
%
%             h.design.Xd  = [Xd ones(size(X,1),1)];
%
%             h.design.XX = h.design.Xd'*h.design.Xd;
end