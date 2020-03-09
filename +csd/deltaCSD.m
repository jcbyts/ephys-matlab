function CSD = deltaCSD(pot, varargin)
% calculate current source density using delta methods

import csd.*

ip = inputParser();
ip.addParameter('el_pos', (1:size(pot,1))'*.1) % mm
ip.addParameter('ex_cond', 0.3)% [S/m]
ip.addParameter('top_cond', 0.3) % [S/m]
ip.addParameter('diam', 0.5) % mm
ip.addParameter('Vaknin', false)
ip.addParameter('filter1', 0.5)
ip.addParameter('filter2', .23)
ip.addParameter('dt', 1)
ip.addParemeter('verbose', false)

ip.parse();

pot0 = pot;

% electrical parameters:
cond = ip.Results.ex_cond;
if cond<=0
    errordlg('ex. cond. has to be a positive number');
    return
end

% filter parameters:
b0 = ip.Results.filter1;
b1 = ip.Results.filter2;
if b0+2*b1 == 0 && b1~=0
    errordlg('Singularity: b0+2*b1 cannot equal zero.');
    return
end


% size, potential (m1 has to equal number of electrode contacts)
[m1,m2] = size(pot);

% electrode parameters:
el_pos = ip.Results.el_pos;
el_pos = el_pos * 1e-3; % mm -> m
el_pos_plot = el_pos(2:length(el_pos)-1); % if not Vaknin electrodes
N = length(el_pos);
h = mean(diff(el_pos));
if m1~=N
    errordlg(['Number of electrode contacts has to equal number of rows in potential matrix. Currently there are ',...
        num2str(N),' electrodes contacts, while the potential matrix has ',num2str(m1),' rows.']) 
    return
end

diam     = ip.Results.diam;
cond_top = ip.Results.top_cond;
if cond_top~=cond && (el_pos~=abs(el_pos) || length(el_pos)~=length(nonzeros(el_pos)))
    errordlg('Electrode contact positions must be positive when top cond. is different from ex. cond.')
    return;
end

% --- Compute delta CSD


% compute delta iCSD:
CSD = F_delta(el_pos,diam,cond,cond_top)^-1*pot;

if b1~=0 %filter iCSD
  [n1,n2]=size(CSD);            
  CSD_add(1,:) = zeros(1,n2);   %add top and buttom row with zeros
  CSD_add(n1+2,:)=zeros(1,n2);
  CSD_add(2:n1+1,:)=CSD;        %CSD_add has n1+2 rows
  CSD = S_general(n1+2,b0,b1)*CSD_add; % CSD has n1 rows
end;

% plot CSD
if ip.Results.verbose
    plot_CSD(CSD,el_pos,ip.Results.dt,1,0) %length(el_pos) must equal rows of CSD! 
end