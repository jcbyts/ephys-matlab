function new_CSD_matrix = stepCSD(pot, varargin)
% calculate inverse current source density

import csd.*

ip = inputParser();
ip.addParameter('el_pos', (1:size(pot,1))'*.1) % mm
ip.addParameter('ex_cond', 0.3)% [S/m]
ip.addParameter('top_cond', 0.3) % [S/m]
ip.addParameter('diam', 0.5) % mm
ip.addParameter('gauss_sigma', 0.1)
ip.addParameter('dt', 1)
ip.addParameter('verbose', false)

ip.parse();

% electrical parameters:
cond = ip.Results.ex_cond;
if cond<=0
    errordlg('ex. cond. has to be a positive number');
    return
end

% filter parameters:
gauss_sigma = ip.Results.gauss_sigma*1e-3;
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
if gauss_sigma<0
    errordlg('The gaussian filter width cannot be negative.')
    return
end;


% size, potential (m1 has to equal number of electrode contacts)
[m1,~] = size(pot);

% electrode parameters:
el_pos = ip.Results.el_pos;
el_pos = el_pos * 1e-3; % mm -> m
N = length(el_pos);
if m1~=N
    errordlg(['Number of electrode contacts has to equal number of rows in potential matrix. Currently there are ',...
        num2str(N),' electrodes contacts, while the potential matrix has ',num2str(m1),' rows.']) 
    return
end

diam     = ip.Results.diam*1e-3;
cond_top = ip.Results.top_cond;
if cond_top~=cond && (el_pos~=abs(el_pos) || length(el_pos)~=length(nonzeros(el_pos)))
    errordlg('Electrode contact positions must be positive when top cond. is different from ex. cond.')
    return;
end


% compute step iCSD:
CSD = F_const(el_pos,diam,cond,cond_top)^-1*pot;

%% make CSD continous (~200 points):
el_pos = el_pos(:)';

le = length(el_pos);
h = el_pos(2)-el_pos(1);
first_z = el_pos(1)-h/2; %plot starts at z1-h/2;
mfactor = ceil(200/le);
minizs = 0:h/mfactor:(mfactor-1 )*h/mfactor;
for i=1:size(CSD,1) % all rows
    zs((1:mfactor)+mfactor*(i-1)) = first_z+(i-1)*h+minizs;
    new_CSD_matrix((1:mfactor)+mfactor*(i-1),:) = repmat(CSD(i,:),mfactor,1);
end


if gauss_sigma~=0 %filter iCSD
  [zs,new_CSD_matrix]=gaussian_filtering(zs,new_CSD_matrix,gauss_sigma,filter_range);
end

% plot CSD
if ip.Results.verbose
    plot_CSD2(new_CSD_matrix,zs,ip.Results.dt,1,0) %length(zs) must equal rows of CSD! 
end