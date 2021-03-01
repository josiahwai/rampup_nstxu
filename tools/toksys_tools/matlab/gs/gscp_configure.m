%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gscp_configure
%
%  PURPOSE: Configure gscp with information about the circular plasma model
%
%  INPUTS: config, structure with TokSys tokamak info (tok_data_struct)
%
%  OUTPUTS: Parameters for analysis of circular equilibrium and its response
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander ON 2015-02-27
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(config,'smallest_aminor')
  smallest_aminor = config.smallest_aminor;
end

if isfield(config,'open_field_line_volume')
  open_field_line_volume = config.open_field_line_volume;
end

if isfield(config,'open_field_line_rmaxis')
  open_field_line_rmaxis = config.open_field_line_rmaxis;
end


ds = [];

% Parameters that define the circle
ds.r0 = 'Major radius of circle center, by definition';
r0 = 1;
ds.z0 = 'Height of circle center, by definition';
z0 = 0;
ds.a0 = 'Minor radius of circle, by definition';
a0 = 0;

% Derived circle parameters
ds.rmin = 'Minimum major radius of boundary';
rmin = r0-a0;
ds.rmax = 'Maximum major radius of boundary';
rmax = r0+a0;
ds.zmin = 'Minimum height of boundary';
zmin = z0-a0;
ds.zmax = 'Maximum height of boundary';
zmax = z0+a0;

% The horizontal moving coordinates
ds.nh = 'number of flux points from rmin to rmax';
nrh = 2*nr-1; % number of flux points from rmin to rmax
ds.rh = 'nh major radii from rmin to rmax';
rh = linspace(1,2,nrh);
ds.drh = 'distance between points in vector rh';
drh = 1/(nrh-1);

% The p-vector contains the nh fluxes (unit Wb)
% AND circle definition: r0, z0, a0
ds.np = 'length of p-vector';
np = nrh+3;
ds.p = 'p contains [psih r0 z0 a0]';
p = zeros(1,np);
ds.dp = 'A change of p';
dp = zeros(np,1);
ds.dperr = 'Errors in p';
dperr = zeros(np,1);

% Flux variables continued (h = horizontal line through midplane)
ds.psih = 'Total flux at points rh [Wb]';
psih = zeros(1,nrh);
ds.psihapp = 'Flux from external conductors at points rh [Wb]';
psihapp = zeros(1,np);
ds.psihapp_r = 'Derivative of psihapp w.r.t. rh';
psihapp_r = zeros(1,np);
ds.psihapp_z = 'Derivative of psihapp w.r.t. z0';
psihapp_z = zeros(1,np);
ds.dpsihappdp = 'Derivative of psihapp w.r.t. p';
dpsihappdp = zeros(nrh,np);
ds.dpsihappdis = 'Derivative of psihapp w.r.t. conductor currents (is)';
dpsihappdis = zeros(nrh,nc+nv);
psihbar = zeros(1,nrh);
psihpla = zeros(1,nrh);
dpsihpladp = zeros(nrh,np);
dpsihplads = zeros(np,ns);
psiherr = zeros(1,nrh);

psierr = zeros(1,nrh);

rconti = zeros(1,nr);
rconto = zeros(1,nr);
drcontidp = zeros(nr,np);
drcontodp = zeros(nr,np);

Rmajor = zeros(1,nr);
dRmajordp = zeros(nr,np);
rminor = zeros(1,nr);
drminordp = zeros(nr,np);

L = zeros(1,nr);
dLdp = zeros(nr,np);
A = zeros(1,nr);
dAdp = zeros(nr,np);
V = zeros(1,nr);
dVdp = zeros(nr,np);


I = zeros(nr,1);
dIdp = zeros(nr,np);
dIds = zeros(nr,ns);

W = zeros(nr,1);
dWdp = zeros(nr,np);
dWds = zeros(nr,ns);

T = zeros(1,nr);
dTdp = zeros(nr,np);
dTds = zeros(nr,ns);

rhot = zeros(1,nr);
drhotdp = zeros(nr,np);
drhotds = zeros(nr,ns);

qpsi = zeros(1,nr);
dqpsidp = zeros(nr,np);
dqpsids = zeros(nr,ns);

% Scalar quantities
cpasma = 0;
dcpasmads = zeros(1,ns);
dcpasmadp = zeros(1,np);
integralRjphidA = 0;
dintegralRjphidAds = zeros(1,ns);
dintegralRjphidAdp = zeros(1,np);
rcur = 0;
drcurds = zeros(1,ns);
drcurdp = zeros(1,np);
zcur = 0;
dzcurdp = zeros(1,np);
ipjdA = 0;
dipjdAds = zeros(1,ns);
dipjdAdp = zeros(1,np);
torflux = 0;
Cl = 0;
dCldp = zeros(1,np);
bp2flx = 0;
dbp2flxds = zeros(1,ns);
dbp2flxdp = zeros(1,np);
betap = 0;
dbetapds = zeros(1,ns);
dbetapdp = zeros(1,np);
cbp2v = zeros(1,nr);
dcbp2vdp = zeros(nr,np);
Bp2V3 = 0;
dBp2V3dp = zeros(1,np);
li3 = 0;
dli3ds = zeros(1,ns);
dli3dp = zeros(1,np);
Bp2V = 0;
dBp2Vdp = zeros(1,np);
li = 1;
dlids = zeros(1,ns);
dlidp = zeros(1,np);
drbdefdp = zeros(1,np);
dzbdefdp = zeros(1,np);
nbdef = 0;
dnbdefdp = zeros(1,np);
dpsibrydp = [1 zeros(1,np-1)];
dpsimagdp = zeros(1,np);

% Plasma response to external fluxes and internal profiles
cpM = zeros(np);
cpMi = zeros(np);
dpdx = zeros(np,nx);

dcpasmadx = zeros(1,nx);
dlidx = zeros(1,nx);
dbetapdx = zeros(1,nx);

% For constraining solution to a cpasma, li, betap combination
dcpasmadx3 = zeros(1,nc+nv+4);
dlidx3 = zeros(1,nc+nv+4);
dbetapdx3 = zeros(1,nc+nv+4);
dx1dx3 = zeros(nc+nv+4);
dx3dx1 = zeros(nc+nv+4);
dxdx3 = zeros(nx,nc+nv+4);
dxdx1 = zeros(nx,nc+nv+4);

descriptions = ds;

% Miscellaneous variables 
onep = ones(1,np);
log10e = 1/log(10);
na = nr-1;
vx = cos((1:na-1)*twopi/(na-1));
vy = sin((1:na-1)*twopi/(na-1));
vx2 = vx.^2;
vy2 = vy.^2;
dgdp = zeros(1,np);
dhdp = zeros(1,np);
dxZdp = zeros(1,np);
i4 = nrh+(3:6); % Indices to rmax, zmax, rmin, zmin
