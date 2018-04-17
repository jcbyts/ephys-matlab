function sp = getSpikesFromKilo(ops, info)
% getSpikesFromKilo gets the spike sorting output from Kilosort
% requires the spikeTools from cortexlab github

if isa(ops, 'table')
    assert(any(strcmp(ops.SpikeSorting, 'Kilo')), 'Kilosort has not been run on this session')
    ops = io.loadOps(ops);
    sp = [];
    for i = 1:numel(ops)
        try
            sp = [sp get_spikes_kilo_helper(ops(i))];
        end
    end
    
else
    if ~exist('info', 'var')
        info = load(fullfile(ops.root, 'ephys_info.mat'));
    end
    sp = get_spikes_kilo_helper(ops, info);
end


function sp = get_spikes_kilo_helper(ops, info)

if nargin < 2
    info = load(fullfile(ops.root, 'ephys_info.mat'));
end

load(ops.chanMap)
folderNames = {ops.root};

Fs = ops.fs;
%%
yc = ycoords(connected); xc = xcoords(connected);
yc=yc(:);
xc=xc(:);

for f = 1:numel(folderNames)    

% clu is a length nSpikes vector with the cluster identities of every
% spike
clu = readNPY(fullfile(folderNames{f},  'spike_clusters.npy'));

% ss is a length nSpikes vector with the spike time of every spike (in
% samples)
ss = readNPY(fullfile(folderNames{f},  'spike_times.npy'));

% convert to times in seconds
% st = double(ss)/Fs;
st = io.convertSamplesToTime(double(ss), Fs, info.timestamps(:), info.fragments(:));

% spikeTemplates is like clu, except with the template numbers rather than
% cluster numbers. Each spike was extracted by one particular template
% (identified here), but when templates were merged in the manual sorting,
% the spikes of both take on a new cluster identity in clu. So
% spikeTemplates reflects the original output of the algorithm; clu is the
% result of manual sorting. 
spikeTemplates = readNPY(fullfile(folderNames{f},  'spike_templates.npy')); % note: zero-indexed

% tempScalingAmps is a length nSpikes vector with the "amplitudes":
% each spike is extracted by a particular template after scaling the
% template by some factor - that factor is the amplitude. The actual
% amplitude of the spike is this amplitude multiplied by the size of the
% template itself - we compute these later. 
tempScalingAmps = readNPY(fullfile(folderNames{f},  'amplitudes.npy'));

% cids is a length nClusters vector specifying the cluster IDs that are used
% cgs is a length nClusters vector specifying the "group" of each cluster:
% 0 = noise; 1 = MUA; 2 = Good; 3 = Unsorted
% Noise spikes should be excluded from everything; we do this in a moment. 
% Both MUA and Unsorted reflect real spikes (in my judgment) but ones that
% couldn't be isolated to a particular neuron. They could include in
% analyses of population activity, but might include spikes from multiple
% neurons (or only partial spikes of a single neuron. Good clusters are
% ones I judged to be well-isolated based on a combination of subjective
% criteria: how clean the refractory period appeared, how large the spike
% amplitudes were, how unique the waveform shapes were given the
% surrounding context. These things can be quantified using the code at
% (https://github.com/cortex-lab/sortingQuality), but I have not done so
% when making the "Good" versus "MUA" decision. I think that due to the
% data quality and density of sites on these probes, the isolation quality
% is in general extremely good. However you should certainly interpret
% "single neurons" with caution and very carefully examine any effects you find
% that might be corrupted by poor isolation. 
%
% Note this function is included in the folder with this data (and also the
% "spikes" repository)

[cgs2, uQ, cR, isiV] = sqKilosort.computeAllMeasures(folderNames{f});


% if ~exist('cgs', 'var')
%     cgs = cgs2;
% end

% temps are the actual template waveforms. It is nTemplates x nTimePoints x
% nChannels (in this case 1536 x 82 x 374). These should be basically
% identical to the mean waveforms of each template
temps = readNPY(fullfile(folderNames{f}, 'templates.npy'));

% The templates are whitened; we will use this to unwhiten them into raw
% data space for more accurate measurement of spike amplitudes; you would
% also want to do the same for spike widths. 
winv = readNPY(fullfile(folderNames{f},  'whitening_mat_inv.npy'));

% compute some more things about spikes and templates; see function for
% documentation
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW] = ...
    templatePositionsAmplitudes(temps, winv, yc, spikeTemplates, tempScalingAmps);

% convert to uV according to the gain and properties of the probe
% 0.6 is the range of voltages acquired (-0.6 to +0.6)
% 512 is the bit range (-512 to +512, 10bits)
% 500 is the gain factor I recorded with 
% 1e6 converts from volts to uV
% spikeAmps = spikeAmps*0.6/512/500*1e6;
spikeAmps = spikeAmps * info.bitVolts;
fname = fullfile(folderNames{f},  'cluster_groups.csv');
if exist(fname, 'file')
    [cids, cgs] = readClusterGroupsCSV(fname);
    
else
    cids = unique(clu);
    cgs = cgs2;
end
    
    spikeCounts = zeros(1, numel(cids));
    for iUnit = 1:numel(cids)
        spikeCounts(iUnit) = sum(clu==cids(iUnit));
    end
    
    
    
%     % find and discard spikes corresponding to noise clusters
%     noiseClusters = cids(cgs==0);
%     
%     st = st(~ismember(clu, noiseClusters));
%     ss = ss(~ismember(clu, noiseClusters));
%     spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
%     tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));
%     clu = clu(~ismember(clu, noiseClusters));
%     cgs = cgs(~ismember(cids, noiseClusters));
%     cids = cids(~ismember(cids, noiseClusters));
%     savefname = 'sp-Kilo.mat';
    
notNoiseClusters = cids(cgs~=0 & cids ~=0 & spikeCounts~=0);
bigix   = ismember(clu, notNoiseClusters);
smallix = ismember(cids, notNoiseClusters);
st = st(bigix);
ss = ss(bigix);
spikeTemplates   = spikeTemplates(bigix);
tempScalingAmps  = tempScalingAmps(bigix);
spikeAmps        = spikeAmps(bigix);
spikeDepths      = spikeDepths(bigix);

clu = clu(bigix);

cgs  = cgs(smallix);
cids = cids(smallix);
uQ   = uQ(smallix);
cR   = cR(smallix);
isiV = isiV(smallix);

savefname = 'sp-Kilo.mat';
    
% else
%     savefname = 'sp-Kilo.mat';
%     cids = unique(clu);
% end


sp(f).name = folderNames{f};
sp(f).clu = clu;
sp(f).ss = ss;
sp(f).st = st;
sp(f).spikeTemplates = spikeTemplates;
sp(f).tempScalingAmps = tempScalingAmps;
sp(f).cgs = cgs;
sp(f).cids = cids;
sp(f).yc = yc;
sp(f).xc = xc;
sp(f).ycoords = ycoords;
sp(f).xcoords = xcoords;
sp(f).temps = temps;
sp(f).spikeAmps = spikeAmps;
sp(f).templateYpos = templateYpos;
sp(f).tempAmps = tempAmps;
sp(f).spikeDepths = spikeDepths;
sp(f).tempsUnW = tempsUnW;

sp(f).cgs2 = cgs2;
sp(f).uQ = uQ;
sp(f).cR = cR;
sp(f).isiV = isiV;

end


%% depths and amplitudes of clusters (as the mean depth and amplitude of all of their constituent spikes)
% get firing rates here also
recordingDuration = sp(1).st(end)-sp(1).st(1);

for s = 1:numel(folderNames)
    clu = sp(s).clu;
    sd = sp(s).spikeDepths;
    sa = sp(s).spikeAmps;    
    
    % using a super-tricky algorithm for this - when you make a sparse
    % array, the values of any duplicate indices are added. So this is the
    % fastest way I know to make the sum of the entries of sd for each of
    % the unique entries of clu
%     spikeCounts = zeros(numel(cids),1);
%     for iUnit = 1:numel(cids)
%         spikeCounts(iUnit) = sum(clu==cids(iUnit));
%     end
    [cids, spikeCounts] = countUnique(clu);    
    q = full(sparse(double(clu+1), ones(size(clu)), sd));
    q = q(cids+1);
    clusterDepths = q./spikeCounts; % had sums, so dividing by spike counts gives the mean depth of each cluster
    
    q = full(sparse(double(clu+1), ones(size(clu)), sa));
    q = q(cids+1);
    clusterAmps = q./spikeCounts;
    
    sp(s).clusterDepths = clusterDepths';
    sp(s).clusterAmps = clusterAmps';
    sp(s).firingRates = spikeCounts'./recordingDuration;
    
    sptmp = sp(s);
    save(fullfile(folderNames{s}, savefname), '-v7.3', '-struct', 'sptmp');
    
end

