function [OE2PTBfit, OE2PTB,PTB2OE, maxreconstructionerror ] = sync2OeClock(PDS, filenameE)
% SYNCOPENEPHYSCLOCK synchronizes the pdlaps PTB clock with open-ephys recording
% Inputs:
%   PDS       - PDS struct (or cell-array of PDS structs)
%   filenameE - string (or cell-array of strings) path to OE events file (*.kwe, *.events)
% Outputs:
%   OE2PTBfit - parameters of fit
%   OE2PTB@function_handle - converts OE time to PTB time
%   PTB2OE@function_handle - convers PTB time to OE time
%   maxreconstructionerror - maximum error in the synchonization
% Example Call:
%   [OE2PTBfit, OE2PTB,PTB2OE, maxreconstructionerror ] = syncOpenEphysClock(PDS, filenameE)

% 2017.08.14     jly     modified Jonas' version

% handling of multiple files
if ~iscell(PDS)
    PDS={PDS};
end

if ~iscell(filenameE)
    filenameE={filenameE};
end

nPdsFiles=numel(PDS);
nOeFiles=numel(filenameE);

bitNumber=[];
timestamps=[];
highlow=[];

for kOeFile=1:nOeFiles
    [~, ~, ext]=fileparts(filenameE{kOeFile});
    switch ext
        case '.kwe'
            tmp_timestamps = hdf5read(filenameE{kOeFile}, '/event_types/TTL/events/time_samples');
            tmp_highlow = hdf5read(filenameE{kOeFile}, '/event_types/TTL/events/user_data/eventID');
            tmp_bitNumber = hdf5read(filenameE{kOeFile}, '/event_types/TTL/events/user_data/event_channels');
            timestamps=[timestamps; tmp_timestamps(:)]; %#ok<*AGROW>
            bitNumber=[bitNumber; tmp_bitNumber(:)];
            highlow=[highlow; tmp_highlow(:)];
        case '.events'
            [tmp_bitNumber, tmp_timestamps, info]=load_open_ephys_data_faster(filenameE{kOeFile});
            tmp_highlow=info.eventId;
            timestamps=[timestamps; tmp_timestamps(:)]; %#ok<*AGROW>
            bitNumber=[bitNumber; tmp_bitNumber(:)];
            highlow=[highlow; tmp_highlow(:)];
    end
end

strobeSet=find(bitNumber==7 & highlow==1);
strobeUnset=find(bitNumber==7 & highlow==0);
strobeUnset=[1; strobeUnset];
% extract strobe values
value=nan(size(strobeSet));
for iStrobe=1:length(strobeSet)
     ts=timestamps <= timestamps(strobeSet(iStrobe)) & timestamps >= timestamps(strobeUnset(iStrobe)) & bitNumber~=7;
     value(iStrobe)=sum(2.^bitNumber(ts) .* highlow(ts));    
end

sixletsOE=fliplr(conv2(value,eye(6)));
sixletsOE=sixletsOE(6:end,:);
sixletsPTB=[];
sixletsPTBts=[];
sixletsDPts=[];

for kPtbFile=1:nPdsFiles
    tmp_sixletsPTB=cellfun(@(X) X.unique_number, PDS{kPtbFile}.data,'uniformOutput',false);
    tmp_sixletsPTB=mod(vertcat(tmp_sixletsPTB{:}),2^7);
	tmp_sixletsPTBts=cellfun(@(X) X.datapixx.unique_number_time(:,1), PDS{kPtbFile}.data,'uniformOutput',false);
    tmp_sixletsPTBts=[tmp_sixletsPTBts{:}]';
    tmp_sixletsDPts=cellfun(@(X) X.datapixx.unique_number_time(:,2), PDS{kPtbFile}.data,'uniformOutput',false);
    tmp_sixletsDPts=[tmp_sixletsDPts{:}]';
    sixletsPTB=[sixletsPTB; tmp_sixletsPTB];
	sixletsPTBts=[sixletsPTBts; tmp_sixletsPTBts];
    sixletsDPts=[sixletsDPts; tmp_sixletsDPts];
end


[~, hind]=ismember(datenum(sixletsPTB),datenum(sixletsOE));
if any(hind==0)
    warning('sync2OeClock: some strobes were missed')
end
goodPTBindex=hind~=0;
hind(hind==0)=[];
sixletsOEts=double([timestamps(strobeSet(hind)) timestamps(strobeSet(hind+1)) timestamps(strobeSet(hind+2)) timestamps(strobeSet(hind+3)) timestamps(strobeSet(hind+4)) timestamps(strobeSet(hind+5))]);
sixletsPTBts=sixletsPTBts(goodPTBindex,:);
sixletsDPts=sixletsDPts(goodPTBindex,:);


OE2PTBfit=[sixletsOEts(:) ones(numel(sixletsOEts),1)]\sixletsPTBts(:);
OE2PTB=@(x) x*OE2PTBfit(1) + OE2PTBfit(2);
PTB2OE=@(x) (x - OE2PTBfit(2))/OE2PTBfit(1);

OE2DPfit=[sixletsOEts(:) ones(numel(sixletsOEts),1)]\sixletsDPts(:);
OE2DP=@(x) x*OE2DPfit(1) + OE2DPfit(2);
% DP2OE=@(x) (x - OE2DPfit(2))/OE2DPfit(1);


% get a reconstruction estimate
% mean(abs(((sixletsDPts(:)-OE2DP(sixletsOEts(:))))))
maxreconstructionerror = max(((sixletsDPts(:)-OE2DP(sixletsOEts(:)))));
