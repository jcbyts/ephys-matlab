function [data,bad] = artifRemovAcrossChannel(data, sizeThresh, chanThresh, win)


bad=unique(bsxfun(@plus, find(sum(abs(data)>sizeThresh,2)>chanThresh), -win:win));
bad(bad<1)=[];
bad(bad>size(data,1))=[];

data(bad,:)=0;