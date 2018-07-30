function thisSession = runSingleChannelSpikeSortThreshold(ops)

if istable(ops)
    thisSession = ops;
    ops = io.loadOps(thisSession);
else
    thisSession = [];
end

info = load(fullfile(ops.root, 'ephys_info.mat'));
data = io.loadRaw(ops, [], true); % load data in mV

% % highpass filter
% %[b,a] = butter(3, 300/30e3*2, 'high');
% bandpass filter between 300 and 6000 Hz
[b,a] = butter(3, [(300/ops.fs*2), (6000/ops.fs*2)]);

size(data)

% detect artifacts
%[data, bad] = preprocess.removeChannelArtifacts(data, 1000, 1, 30, false);
%data(bad,:) = 0;

sp = struct('ss', [], 'clu', []);

nChannels = size(data,1);
Fs = info.sampleRate; % sampling rate
clustOffset = 0;
refms = 1.0;
ref_period = floor(Fs*refms/1000);  % 2 ms
threshold = -4;  % SD of data

for iCh = 1:nChannels
    
  fprintf('Sorting channel %d\n',iCh);
  
  zdata = data(iCh,:)';
  zdata = filter(b,a,zdata);
  zdata = flipud(zdata);
  zdata = filter(b,a,zdata);
  zdata = flipud(zdata);
  SD_data = std(zdata);
   

  fthresh = threshold * SD_data;
  got_thresh = 0;  % cycle command line to pick threshold

  while (got_thresh == 0)

    %********* plot the total data stream
    figure(1); clf
    %xx = 1:30000;
    %plot(data(xx,iCh)); hold on;
    plot(zdata); hold on;
    axis tight;
    ylabel('mV');
    xlabel('samples');
    V = axis;
    Vax = abs( fthresh * 6);  % zoom in near thresh
    axis([V(1) V(2) -Vax Vax]);
    plot([V(1),V(2)],[fthresh,fthresh],'r--');
    
    %******** compute the spike threshold crossing and mean rate
    ss = [];
    s = find( (zdata(2:end ) < fthresh) & ...
               (zdata(1:(end - 1)) >= fthresh) );
    if (~isempty(s))
      s = s(s > 10 & s < size(zdata, 1) - 25);    % remove spikes close to boundaries
      %**** impose refractory
      z = find( ( s(2:end) - s(1:(end-1))) > ref_period);
      if (~isempty(z))
        ss = s(z);
      end
    end
    %*******
    if (~isempty(ss))
      clu = ones(size(ss));  % all same for now, just one MU unit
    else
      clu = [];
    end
    %*******************************************

    mrate = size(ss,1)/(size(zdata,1)/Fs);
    title(sprintf('Mean Rate %5.2f Hz Chan(%d)',mrate,iCh));
    %***************
    
    fprintf('Current threshold: %5.2f std\n',threshold);
    nval = input('Return to accept, or enter new value: ');
    if isempty(nval)
        got_thresh = 1;
        %****** finalize spike structure
       sp.ss = [sp.ss; ss(:)];
       sp.clu = [sp.clu; clu(:)+clustOffset];
    
       clustOffset = max(sp.clu);
       %***************
    else
        threshold = nval;
        fthresh = threshold * SD_data;
    end    
  end  % end of while loop
end
disp('Finished sorting');


%% save spikes
sp.st = io.convertSamplesToTime(sp.ss, Fs, info.timestamps(:), info.fragments(:));
sp.cids = unique(sp.clu)';

load(ops.chanMap)
sp.yc = ycoords;
sp.ycoords = ycoords;
sp.xc = xcoords;
sp.xcoords = xcoords;
n = numel(sp.cids);
sp.cgs = zeros(n,1);
sp.cgs2= zeros(n,1);
sp.cR  = zeros(n,1);
sp.clusterAmps   = zeros(1,n);
sp.clusterDepths = zeros(1,n);
sp.firingRates   = zeros(1,n);
sp.isiV = zeros(1,n);


save(fullfile(ops.root, 'sp-threshold.mat'), '-v7.3', '-struct', 'sp')

if istable(thisSession)
    newThisSession = thisSession;
    
    ssList = thisSession.SpikeSorting{1};
    if ~isnan(ssList)
        ssList = regexp(ssList, ',', 'split');
        ssList = union(ssList, {'threshold'});
        ssList = sprintf('%s,', ssList{:});
        newThisSession.SpikeSorting{1} = ssList(1:end-1);
    else
        ssList = {'threshold'};
        ssList = sprintf('%s,', ssList{:});
        newThisSession.SpikeSorting{1} = ssList(1:end-1);
    end
    
    io.writeMeta(newThisSession, 2)
    thisSession = newThisSession;
end

