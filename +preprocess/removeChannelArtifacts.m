function [data,bad] = removeChannelArtifacts(data, sizeThresh, chanThresh, win)

% verbose = true;


% z = abs(mean(data));
% crossings =  z > sizeThresh;
% crossings = crossings(:);

% if verbose
%     figure; clf
%     plot(z); hold on
%     t = nan(size(crossings)); t(crossings) = 1;
%     plot(t*max(z), '.')
% end
sz = size(data);
if sz(2) > sz(1)
    trpose = true;
    data = data';
else
    trpose = false;
end
crossings = sum(abs(data)>sizeThresh,2)>=chanThresh;
bad=unique(bsxfun(@plus,find(crossings), -win:win));
bad(bad<1)=[];
bad(bad>max(sz))=[];

data(bad,:)=0;

if trpose
    data = data';
end