meta = io.getExperimentsAnd();
thisSession = meta(19,:);
[sessionInfo, ops, info] = io.loadSession(thisSession);

%%

dataRAW = io.loadRaw(ops(1), [], true);
%%
figure(1); clf
plot(dataRAW(1,:))
%%

dataRAW = io.loadRaw(ops(1), 50e5 + [1 100e5], true);

%%
figure(3); clf
plot(bsxfun(@plus, dataRAW', (1:32)*400), 'k')


%%

figure(2); clf
x = dataRAW(14,:);
fs = info.sampleRate;

subplot(3,1,1)
plot((1:numel(x))/fs, x, 'k');
ylim([-1e3 1e3])

subplot(3,1,2:3)


nwin = 32e3;
wind = hanning(nwin);
nlap = ceil(.9*nwin);
freq = 0:2:96;

[S, f, t] = spectrogram(x, nwin, nlap, freq, fs, 'yaxis');

imagesc(t, f, log(abs(S))); axis xy
xlabel('Time')
ylabel('Frequency')
colormap jet