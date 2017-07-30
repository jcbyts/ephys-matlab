
oepath = 'C:\Data\Ellie_3-sessions_2017-07-24_2017-07-26\';

[session, ops, info] = io.loadSession(oepath);

PDS = io.getPds(session);

sp = io.getSpikesFromKilo(ops, info);


%%

% [data, timestamps, bad] = io.loadAndPreprocess(ops);

inds = 1500e3 + (1:1200e3);

Fs = info.sampleRate;

[data, timestamps] = io.loadRaw(ops, inds);


%%


dproc = preprocess.highpass(data, Fs, 300, 7500);

[dproc, bad]  = preprocess.removeChannelArtifacts(dproc, 10, 10, 50);

data(:,bad) = 0;

figure(1); clf
plot(bsxfun(@plus, data', (1:ops.Nchan)*200))

%%

dproc2 = preprocess.removeLineNoise(dproc, Fs, [60 180], 1);




close all
%%
data = io.loadRaw(ops, inds);


