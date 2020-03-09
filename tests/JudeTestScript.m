
SessionNum = size(PDS,1);


for k = 1:SessionNum

   TrialNum = size(PDS{k}.data,2);
   
   for tk = 1:(TrialNum-1)
       
      %screen flip times in ephys time coordinates
      if (isfield(PDS{k}.data{tk},'DotMotionMapping'))
       
         
          
        %tims = PDS{k}.PTB2OE(PDS{k}.data{tk}.timing.flipTimes(1,1:(end-1)));
        %PDS{k}.data{tk}.tims = tims;
        tims = PDS{k}.data{tk}.tims;
        dirs = PDS{k}.data{tk}.DotMotionMapping(1).direction;  % motion directions
        onstim = PDS{k}.data{tk}.DotMotionMapping(1).on;  % on stimulus, 1 or 0
        sptims = sp{1}.st;   % same coordinates as PTB2OE
        spcids = sp{1}.cids;
        spclus = sp{1}.clu;
      
        %******** show raw data
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
      
        disp(sprintf('Trial %d Session %d',tk,k));
        input('stop');
      end
      
   end
end
