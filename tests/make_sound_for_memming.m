
thisSession = meta(20,:);
ops = io.loadOps(thisSession)

%%
ops.fslow = 3e3;
ops.fshigh = 250;

raw = io.loadRaw(ops, 120*30e3 + [1 30e3*10], true, false);
% raw = raw / max(abs(raw(:)));

if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
    [b1, a1] = butter(5, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
else
    [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
end


raw = raw';

datr = filter(b1, a1, raw);
datr = flipud(datr);
datr = filter(b1, a1, datr);
datr = flipud(datr);

datr = datr / max(abs(datr(:)));
ch = 0;
%%
ch = ch + 1;
figure(1); clf
plot(datr(:,ch));
title(ch)

%%
% ch = 8
sound(datr(:,ch), 30e3)

audiowrite(sprintf('channel%02.0f.wav', ch), datr(:,ch), 30e3);
