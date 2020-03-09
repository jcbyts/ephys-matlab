
function plot_direction_select(PDS,sp,zchan)

debug = 0;

%********
gravedigger = 1;  
Prelude = 0.3;  % plot this many secs before onset
Postlude = 0.9; % plot this many secs after onset
Awin = 0.050;
Bwin = 0.550;
%*********************
Smoothing = 15; % in ms
PsthTime.Zero = floor(1000*Prelude);
PsthTime.Awin = floor(1000*(Prelude+Awin));
PsthTime.Bwin = floor(1000*(Prelude+Bwin));
%*********

for z = 1:size(zchan,2)

 StimDirs = [];
 StimSpks = [];
 TrialCnt = 1;

 chan = zchan(z);
 
 %*****
 SessionNum = size(PDS,1);
 for k = 1:SessionNum
   TrialNum = size(PDS{k}.data,2);
   for tk = 1:(TrialNum-1)
       
      %screen flip times in ephys time coordinates
      if (isfield(PDS{k}.data{tk},'DotMotionMapping'))
        
        if (gravedigger == 1)
           tims = PDS{k}.PTB2OE(PDS{k}.data{tk}.timing.flipTimes(1,1:(end-1)));
        else
           tims = PDS{k}.data{tk}.tims;
        end
        dirs = PDS{k}.data{tk}.DotMotionMapping(1).direction;  % motion directions
        onstim = PDS{k}.data{tk}.DotMotionMapping(1).on;  % on stimulus, 1 or 0
        
        spikes = find( sp{1}.clu == chan)';
        
        sptims = sp{1}.st(spikes);   % same coordinates as PTB2OE
        % spcids = sp{1}.cids(spikes);
        spclus = sp{1}.clu(spikes);
      
        %******** show raw data
        if (debug == 1)
          figure(1); hold off;
          H = plot(tims,dirs,'k-');  
          set(H,'Linewidth',2);
          axis tight; hold on;
          H = plot(tims, (100 * onstim),'r-');
          set(H,'Linewidth',2);
      
          colo = 'brgycmkkkkkkkkkkkk';
          % take all clusters, multi-unit activity
          zz = find( (sptims >= tims(1)) & (sptims < tims(end)) );
          if (~isempty(zz))
            for i = zz
              plot([sptims(i),sptims(i)],[0,150],'b-');
            end
          end
          xlabel('Time (secs)');
          V = axis;
          axis([V(1) V(2) 0 360]);
          %***********
          disp(sprintf('Trial %d Session %d',tk,k));
          input('stop');
          %************
        end
        %*******************************
        
        %****** identify onsets ********
        zz = find( onstim(2:end) - onstim(1:(end-1)) == 1);
        onk = zz+1;
        for zk = 1:size(onk,1)
            %********
            ton = tims(onk(zk));
            tstart = (ton - Prelude);
            tend = (ton + Postlude);
            zt = find( (sptims >= tstart) & (sptims < tend) );
            %******
            if ~isempty(zt)
                spt = sptims(zt) - ton;
            else
                spt = [];
            end
            %*************
            StimDirs(TrialCnt) = dirs(onk(zk));
            StimSpks{TrialCnt} = spt;
            TrialCnt = TrialCnt + 1;
            %***********
        end    
        %*************
      end
      
   end
end
%**************************
TrialN = size(StimDirs,2);
T = 1 + floor((Prelude+Postlude)/0.001);
StimRast = zeros(TrialN,T);
for k = 1:TrialN
    spt = StimSpks{k};
    ispt = floor( (spt + Prelude)/0.001 );
    for i = 1:size(ispt,1)
        if (ispt(i) > 0) & (ispt(i) <= T)
           StimRast(k,ispt(i)) = 1;
        end
    end
end
%****

%****** Now build a spike raster from data
udirs = unique( StimDirs );
uspec = mean( abs( udirs(2:end) - udirs(1:(end-1))));
NU = size(udirs,2);
ucounts = zeros(size(udirs));
for cc = 1:size(StimDirs,2)
    idir = find( StimDirs(cc) == udirs );    
    ucounts(idir) = ucounts(idir) + 1;
end
maxcnt = max(ucounts);
%******
spcnts = cell(1,NU);
%*****
icounts = zeros(size(udirs));
xd = [];
yd = [];
%****
for cc = 1:size(StimDirs,2)
    idir = find( StimDirs(cc) == udirs );    
    val = udirs(idir) - (0.45*uspec) + (0.90*uspec*(icounts(idir)/maxcnt));
    icounts(idir) = icounts(idir) + 1;
    spt = StimSpks{cc};
    xd = [xd ; spt];
    yd = [yd ; val*ones(size(spt))];
    %********
    zz = find( (spt >= Awin) & (spt < Bwin) );
    tcnt = (size(zz,1)/(Bwin-Awin));
    spcnts{idir} = [spcnts{idir} tcnt];
    %********
end

%***** compute mean via direction
for di = 1:NU
    uu(di) = mean( spcnts{di} );
    su(di) = std( spcnts{di} ) / sqrt( size(spcnts{di},2) );
end
%***************************
figure(1); 
subplot('position',[0.1 0.4 0.4 0.5]); hold off;
H = plot(xd,yd,'k.'); hold on;
set(H,'Markersize',2);
plot([Awin,Awin],[-60,360],'b-');
plot([Bwin,Bwin],[-60,360],'b-');
axis([-Prelude Postlude -60 360]);
ylabel('Direction (degs)');
title(sprintf('Channel %d',chan));
%*****
subplot('position',[0.1 0.1 0.4 0.25]); hold off;
%************
plot_psth(StimRast,Smoothing,PsthTime);
xlabel('Time (ms)');
ylabel('Mean Rate (hz)');
%********
subplot('position',[0.6 0.6 0.35 0.35]); hold off;
errorbar(udirs,uu,su,'b'); hold on;
axis tight;
V = axis;
axis([V(1) V(2) V(3) V(4)*1.5]);
xlabel('Direction');
ylabel('Mean Rate (hz)');
%*******
hu = uu+su;
lu = uu-su;
subplot('position',[0.6 0.1 0.35 0.35]); hold off;
polar([(udirs*2*pi/360),0],[hu,hu(1)],'b-'); hold on;
H = polar([(udirs*2*pi/360),0],[uu,uu(1)],'b-');
set(H,'Linewidth',2);
polar([(udirs*2*pi/360),0],[lu,lu(1)],'b-'); 
%*****************************************
input('Continue?');

end

return;


function plot_psth(SpikeData,Smoothing,PsthTime)
         
   %*********************
   uu = mean(SpikeData);
   uu = uu * 1000;  % compute mean
   uu = gauss_smooth( uu, Smoothing); % smooth it
   %******* compute the Jacknife
   smoothsub = [];
   N = size(SpikeData,1);
   T = size(SpikeData,2);
   for i = 1:N
    excludeset = [1:(i-1),(i+2):N];
    psth = mean(SpikeData(excludeset,:));
    psth = psth*1000;
    smooth_data = gauss_smooth(psth,Smoothing);  
    smoothsub = [smoothsub; smooth_data];
   end
   u_smooth = mean(smoothsub);
   sem_smooth = std(smoothsub);
   sem_smooth = sem_smooth * sqrt((N-1));   % variance is multiplied by N-1 Jacknife
   %***********************************          
   xx = 1:T;
   uu = mean(SpikeData);
   H = plot(xx,u_smooth,'b-'); hold on;
   set(H,'Linewidth',2);
   plot(xx,(u_smooth+(2*sem_smooth)),'b:');
   plot(xx,(u_smooth-(2*sem_smooth)),'b:');
   axis tight;
   V = axis;
   plot([PsthTime.Zero,PsthTime.Zero],[V(3),V(4)],'k-');
   plot([PsthTime.Awin,PsthTime.Awin],[V(3),V(4)],'b-');
   plot([PsthTime.Bwin,PsthTime.Bwin],[V(3),V(4)],'b-');
   %************
   xlabel('Time(ms)'); 
   ylabel('Rate(hz)');
   
return;


%***** from gauss_smooth_KSW (first pass we made on Sunwoo's code)   
function smo = gauss_smooth(psth,Gsig)

% Make the number of samples depending on the gaussian window size
gaussian_filter_size = 2*Gsig-1; % if Gsig = 10, 19 samples total
                                 % 9 left & 9 right from the mean

% Make smoothing kernel using gaussian filter
for i = 1:gaussian_filter_size
    gauss  = exp(-(((i-Gsig).^2)./(2*Gsig^2)));
    gauss_filter(i,:) = gauss;
end

% Normalize the gaussian filter
gauss_smooth = gauss_filter/sum(gauss_filter);

psth_size    = length(psth);
filter_size  = length(gauss_smooth);
filter_cent = floor((filter_size+1)/2);

for i=1:psth_size   % size_smooth
    
    % Always 0 for the initial value (only sum from product of two vectors)
    smo(i) = 0;
    nomo(i) = 0;
    
    % Apply filter to data
    for j = 1:filter_size
         diff = (j-filter_cent);   % this coordinate, goes - 3*sig up to 3*sig
         samp = (i+diff);
         if ( (samp >= 1) && (samp <= psth_size) )
             smo(i) = smo(i) + (psth(samp) * gauss_smooth(j));
             nomo(i) = nomo(i) + gauss_smooth(j);
         end       
    end
    %********
    if (nomo(i) > 0)
        smo(i) = smo(i) / nomo(i);
    end
    %***********
end

return;






