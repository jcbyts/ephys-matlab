function [y, pss] = removeLineNoise(dat,sr,freqs,freqrange,pss,showoutput)
% Removes line noise from a signal
% [y, pss] = removeLineNoise(dat,sr,freqs,freqrange,pss,showoutput)
% Input:
%   data  [n x 1]    vector containing continuous data
%   sr    [1 x 1] 	 the sampling rate of the signal
%   freqs [1 x m]    a vector of frequencies around which to remove line
%                    noise, for example [60,180]
%   freqrange (optional) a scalar specifying the range of frequencies
%                        to use for fitting a parametric function for
%                        line noise (default: 2Hz)
%   pss (optional)       parameters of filtering
%   showOutput (opt)     if true, show plots of fits (default: false)
%   Algo:
%       1) abs(fft(dat)).
%       2) Extract the vector corresponding to [freq(1)-freqrange, freq(1)+freqrange].
%       3) Fit a function thefun = @(p,x) p(1)*exp(-p(2)*abs(x-p(4)).^(p(5)))+p(3);
%          to this range of freqs.
%       4) Invert the function to obtain a delining filter and apply it.
% copied from p. mineault (xcorr.net)

if nargin < 6
    showoutput = true;
    if nargin < 5
        pss = nan(5,numel(freqs));
        if nargin < 4
            freqrange = 2;
        end
    end
end

sz = size(dat);
if sz(2) > sz(1)
    dat = dat';
    sz = size(dat);
    trpose = true;
else
    trpose = false;
end

y = nan(sz);

for i = 1:sz(2)
    [y(:,i), pss] = removeLineNoiseChannel(dat(:,i), sr, freqs, freqrange, pss, showoutput);
end

if trpose
    y = y';
end

end

function [y, pss] = removeLineNoiseChannel(dat,sr,freqs,freqrange,pss,showoutput)

ae = [];
if mod(length(dat),2) == 1
    ae = dat(end);
    dat = dat(1:end-1);
end
fftdat = fft(double(dat));
a = abs(fftdat);

%Now remove line noise from datlo to obtain y
%
%Eliminate line noise at target frequencies
thefilt = ones(size(a));

winlen = round(length(dat)/sr*freqrange);

%Fit a curve to this chunk of frequencies
opts = optimset('Display','Off','Jacobian','on','Algorithm','levenberg-marquardt');


n = 1;
for tgtr = freqs
    
    peak = tgtr/sr*length(dat);
    
    rg = round(((peak-winlen):(peak+winlen)))';
    datrg = a(rg);
    
    x = (-winlen:winlen)'/winlen*freqrange;
    
    
    %Only adjust a few parameters at a time
    %convergence is better this way
    %everything but the exponent
    %Find the peak
    datrgsm = conv(datrg,ones(21,1),'same');
    [~,peakloc] = max(datrgsm);
    
    %Set the initial parameters
    x0 = [max(datrg)-median(datrg),1/.2^2,median(datrg),(peakloc-1-winlen)/winlen*freqrange,1]';
    
    if ~any(isnan(pss(:,n)))
        x0([2,4,5]) = pss([2,4,5],n);
    end
    
    
    [ps] = lsqcurvefit(@(x,y) thefun([x;x0(5)],y,[1;1;1;1;0]),x0(1:4),x,datrg,[],[],opts);
    xd = ps(4);
    
    %Everything but the center
    [ps] = lsqcurvefit(@(x,y) thefun([x(1:3);xd;x(4)],y,[1;1;1;0;1]),[ps(1:3);x0(5)],x,datrg,[],[],opts);
    
    %Everything
    [ps] = lsqcurvefit(@(x,y) thefun(x,y,[1;1;1;1;1]),[ps(1:3);xd;ps(4)],x,datrg,[],[],opts);
    
    pss(:,n) = ps;
    
    %Good, now adjust the filter in this range accordingly
    thefilt(rg) = ps(3)./thefun(ps,x);
    b = thefilt(rg);
    thefilt(end-rg+2) = b;
    
    if showoutput
        subplot(length(freqs),1,n);
        plot(x+tgtr,datrg,x+tgtr,thefun(ps,x));
        title(sprintf('%3.1f Hz',tgtr));
        drawnow;
        
        [ps(2),ps(4),ps(5)]
    end
    n=n+1;
end

% y is datlo with line noise removed
a = fftdat.*thefilt;
y = [real(ifft(a));ae];

end

function [y,J] = thefun(p,x,mask)
E = abs(x-p(4)).^(p(5));
M = exp(-p(2)*E);
y = p(1)*M+p(3);
if nargout > 1
    J = [ M,...
        -p(1)*E.*M,...
        ones(size(x)),...
        p(1)*p(2)*p(5)*sign(x-p(4)).*abs(x-p(4)).^(p(5)-1).*M,...
        -p(1)*p(2)*p(5)*E.*log(abs(x-p(4))+1e-6).*M];
    J = J(:,mask==1);
end
end