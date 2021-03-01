function [d,e,r,b,p] = gsupdate(inp1,inp2,inp3,inp4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  To initialize: gsupdate(x, e, c) or gsupdate(x, e, c, 1)
%  To update: d = gsupdate(x)       or gsupdate(x, e, c, 0)
%  [d,e,r,b,p] = gsupdate(....) always updates to exactly x
%  d = gsupdate(....) only updates to a close enough x if one exists in memory
%  gsupdate('help') gives a long help
%
%  PURPOSE: Update d,e,r,b,p to state x or d to a near enough x
%
%  INPUTS: x, new state vector
%          e, equilibrium (initialize with gsinit.m)
%          c, configuration created by gsconfig.m          
%
%  OUTPUTS: d, dynamics
%           e, equilibrium
%           r, response
%           b, boundary
%           p, profiles
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  WRITTEN BY:  Anders Welander ON 2016-11-17
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% c = configuration, e = equilibrium, r = response, d = dynamics

% Remember d and some of e,r for nn equilibria, ii points to latest, c_now=c.now
persistent ii nn dd ee rr c_now

% Remember configuration
persistent c

% Number of calls to gsupdate and number of complete response calculations
persistent nupdate ncalc

if isempty(nupdate)
  nupdate = 0;
end
if isempty(ncalc)
  ncalc = 0;
end

if nargin == 1 & ischar(inp1)
  if strcmp(upper(inp1),'HELP')
    unix('more gsupdate_help.txt');
  elseif strcmp(inp1,'ii')
    d = ii;
  elseif strcmp(inp1,'nn')
    d = nn;
  elseif strcmp(inp1,'dd')
    d = dd;
  elseif strcmp(inp1,'ee')
    d = ee;
  elseif strcmp(inp1,'rr')
    d = rr;
  elseif strcmp(inp1,'c')
    d = c;
  end
  return
end
if nargin == 2 & ischar(inp1)
  if strcmp(inp1,'ii')
    ii = inp2;
  elseif strcmp(inp1,'nn')
    nn = inp2;
  elseif strcmp(inp1,'dd')
    dd = inp2;
  elseif strcmp(inp1,'ee')
    ee = inp2;
  elseif strcmp(inp1,'rr')
    rr = inp2;
  elseif strcmp(inp1,'c')
    c = inp2;
  end
  return
end

nupdate = nupdate+1;

% For cubic Hermite splines
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;

% Vacuum permeability
mu0 = 4e-7*pi;

% Configuration data (created by gsconfig, c.info describes the fields of c)
if nargin > 2
  c = inp3;
end
i16 = c.i16;
nz = c.nz;
nr = c.nr;
rg = c.rg;
zg = c.zg;
dr = c.dr;
dz = c.dz;
ngg = c.ngg;
nkp = c.nkp;
nkf = c.nkf;
nic = c.nic;
niv = c.niv;
nih = c.nih;
nsf = c.nsf;
nsp = c.nsp;
nx = c.nx;
ny = c.ny;
nu = c.nu;
ix = c.ix;
iy = c.iy;
iE = c.iE;
mcc = c.mcc;
mcv = c.mcv;
mvv = c.mvv;
mpc = c.mpc;
mpv = c.mpv;
mgg = c.mgg;
psikf = c.psikf;
psikp = c.psikp;
nb = c.nb; % number of points on a contour
np = c.np; % number of contours
nn = c.nn;
nbdtest = c.nbdtest;
nfla = c.nfla;
nflc = c.nflc;
nfltest = c.nfltest;
dpsibarlimit = c.dpsibarlimit;
evolve_option = c.evolve_option;
Pcc = c.Pcc;
Pvc = c.Pvc;
lstarx = c.lstarx;
Rhat = c.Rhat;
Vhat = c.Vhat;

if nargin > 3
  initialization = logical(inp4);
else
  initialization = nargin > 1;
end

% Allocate memory
if isempty(c_now) | c_now ~= c.now
  
  c_now = c.now;
  
  ii = 1;

  initialization = true;

  fl = struct('test',inf(nfltest,nx),'teste',zeros(nfltest,1),'x',zeros(nx,1));

  bd = struct('n',0, 'r',zeros(1,nbdtest), 'z',zeros(1,nbdtest), ...
    'w',zeros(nbdtest,16), 'ii',ones(nbdtest,16), ...
    'psi',zeros(1,nbdtest), 'lim',zeros(1,nbdtest), ...
    'test',inf(nbdtest,nx), 'teste',inf(nbdtest,1), 'x',zeros(nx,1));

  dd = repmat(struct('A', zeros(nx,nx), 'B', zeros(nx,nu), ...
    'C', zeros(ny,nx), 'D', zeros(ny,nu), 'E', zeros(ny,1), ...
    'xdot0', zeros(nx,1), 'y0', zeros(ny,1), 'gamma', 0, ...
    'fl', fl, 'bd', bd, 'de', 1, 'fmax', 1),1,nn);

  ee = repmat(struct('x',zeros(nx,1), 'psizr',zeros(nz,nr), ...
    'rmaxis',nan, 'zmaxis',nan, ...
    'psimag',nan, 'psibry',nan, ...
    'cpasma',0, 'li',nan, 'betap',0, ...
    'alcfs',nan, 'zcur',nan, 'rb',nan(nb,1), 'zb',nan(nb,1) ),1,nn);

  rr = repmat(struct('dpsizrdx',zeros(ngg,nx), 'dpsizrde',zeros(ngg,1), ...
    'drmaxisdx',zeros(1,nx), 'dzmaxisdx',zeros(1,nx), ...
    'dpsimagdx',zeros(1,nx), 'dpsibrydx',zeros(1,nx), ...
    'dcpasmadx',zeros(1,nx), 'dlidx',zeros(1,nx), 'dbetapdx',zeros(1,nx), ...
    'dzcurdx', zeros(1,nx), 'drbdx',zeros(nb,nx), 'dzbdx',zeros(nb,nx)),1,nn);

end

% Here ii points to *latest* output of dd, ee, rr (or 1 if none existed)

if initialization
  iE(:) = 0;
end

% Get the state vector, x
if nargin > 0
  x = inp1;
else
  x = ee(ii).x; % Allow calls with no arguments
end

% The latest dynamics structure, that was returned in the previous call
d = dd(ii);

% Calculate output vector y using the d from the previous call
fl = d.fl;
bd = d.bd;
fltest = d.fl.test*(x-d.fl.x); % validity of linear response
bdtest = d.bd.test*(x-d.bd.x); % test boundary-defining point
flmax1 = max(abs(fltest)); % < 1 for linear approx
bdmax1 = max(bdtest); % > 1 means bdef point jumps
f = 1-max(flmax1,bdmax1); % fraction of d.E to keep
f = max(0,min(d.fmax,f)); % f is in range 0 to d.fmax
yold = d.y0 + d.C*x + d.E*f; % y based on old d

% Look for suitable previously calculated d to use as new d
fs = -inf(1,nn);
for i = 1:nn
  d = dd(i);
  fl = d.fl;
  bd = d.bd;
  fltest = d.fl.test*(x-d.fl.x); % validity of linear response
  bdtest = d.bd.test*(x-d.bd.x); % test boundary-defining point
  flmax = min(inf,max(abs(fltest))); % < 1 for linear approx
  bdmax = min(inf,max(bdtest)); % > 1 means bdef point jumps
  fs(i) = 1-max(flmax,bdmax); % fraction of old y (in d.E) to keep
  if isinf(fs(i))
    break
  end
end
[fnew, ii] = max(fs);

% Now ii points to *best* choice of dd, ee, rr for updating to state x

% The BEST match among previous responses
d = dd(ii);

% Very important that next evaluation in the Simulink subsystem 'linear model'
% gives f > 0, otherwise it will not trigger the update system anymore since
% that only happens when f crosses back to < 0 from having been > 0 at least once.
% A test showed that the x used to calculate the trigger in Simulink is used 
% again the first time new d.A, etc values come in. Therefore the f calculated 
% here will match that next value in 'linear model'. Hence it is sufficient that 
% fnew > 0 to guarantee that f > 0 at least once also in the Simulink 'linear model'.
% Even so, require fnew > 0.2 to avoid creating large E values in code below.
if nargout == 1 & (~initialization & fnew > 0.2 | initialization & fnew == 1)
  % Use the archived d but modify d.E so that: d.y0 + d.C*x + d.E*fnew = yold
  d.E(iE) = (yold(iE) - d.y0(iE) - d.C(iE,:)*x)/fnew;
  d.E(~iE) = 0;
  d.fmax = fnew;
  % Also update dd(ii) with same changes as d, will be used next time to calculate yold
  dd(ii).E = d.E;
  dd(ii).fmax = fnew;
  % Update the special output nupdate
  if isfield(iy,'nupdate')
    d.y0(iy.nupdate) = nupdate;
  end
  % Fast return using archived old calculation of d
  return
end

% Full analysis begins
ncalc = ncalc+1;

% Update psizr to new x
if initialization
  psizr = inp2.psizr;
  rmaxis = nan;
  zmaxis = nan;
else
  dx = x - ee(ii).x;
  psizr = ee(ii).psizr + reshape(rr(ii).dpsizrdx*dx-rr(ii).dpsizrde*d.de,nz,nr);
  rmaxis = ee(ii).rmaxis + rr(ii).drmaxisdx*dx;
  zmaxis = ee(ii).zmaxis + rr(ii).dzmaxisdx*dx;
end

% Full analysis entails:
% step 1, decide pcurrt with this new psizr and x
% step 2, decide dpsizrdxe
% step 3, decide A,B,C,D from differential equations on the form V = R*x+L*xdot

% Extract the different parts of the state data
ic = x(ix.ic);
iv = x(ix.iv);
ih = x(ix.ih);
sf = x(ix.sf);
sp = x(ix.sp);

% Analyze boundary with new psizr
try
  b = gsboundary(c,psizr,rmaxis,zmaxis);
  if ~b.plasma
  save bvac c x psizr rmaxis zmaxis
  end
catch
  save btrouble c x psizr rmaxis zmaxis
  wait('error in gsboundary')
end

% Structure with boundary-defining points
bd.n = b.bd.n;
bd.r = rg(1) + dr*(b.bd.r-1);
bd.z = zg(1) + dz*(b.bd.z-1);
bd.w = b.bd.w;
bd.ii = b.bd.ii;
bd.psi = b.bd.psi;
bd.lim = b.bd.lim;
lim = bd.lim(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary data (b.info describes the fields of b)

n = b.n;

% Cell data
i49 = b.resp.di7x7';
Ac = b.Acell;
Acw = b.resp.Acell.w;
Acb = b.resp.Acell.b;
RAc = b.RAcell;
RAcw = b.resp.RAcell.w;
RAcb = b.resp.RAcell.b;
ARc = b.ARcell;
ARcw = b.resp.ARcell.w;
ARcb = b.resp.ARcell.b;
yc = b.psizrc;
ycw = b.resp.psizrc.w;
ycb = b.resp.psizrc.b;
ycr = b.psizrcr;
ycrw = b.resp.psizrcr.w;
ycrb = b.resp.psizrcr.b;
ycz = b.psizrcz;
yczw = b.resp.psizrcz.w;
yczb = b.resp.psizrcz.b;
rgc = b.rgc;
rgcw = b.resp.rgc.w;
rgcb = b.resp.rgc.b;
zgc = b.zgc;
zgcw = b.resp.zgc.w;
zgcb = b.resp.zgc.b;
igg = b.igg;
egg = b.egg;

% Shape parameters
R8 = b.R8;
Z8 = b.Z8;
rsurf = b.rsurf;
zsurf = b.zsurf;
aminor = b.aminor;
bminor = b.bminor;
elong = b.elong;
triu = b.triu;
tril = b.tril;
squo = b.squo;
squi = b.squi;
sqli = b.sqli;
sqlo = b.sqlo;
drsep = b.drsep;
Cl = b.Cl;
Vtot = b.Vtot;
Atot = b.Atot;
Ltot = b.Ltot;
rsurfp = b.resp.rsurfp;
zsurfp = b.resp.zsurfp;
aminorp = b.resp.aminorp;
bminorp = b.resp.bminorp;
elongp = b.resp.elongp;
triup = b.resp.triup;
trilp = b.resp.trilp;
squop = b.resp.squop;
squip = b.resp.squip;
sqlip = b.resp.sqlip;
sqlop = b.resp.sqlop;
drsepp = b.resp.drsepp;
Clp = b.resp.Clp;
Vtotp = b.resp.Vtotp;
Atotp = b.resp.Atotp;
Ltotp = b.resp.Ltotp;

% Axis and bdef
rmaxis = b.rmaxis;
zmaxis = b.zmaxis;
kmaxis = b.kmaxis;
psimag = b.psimag;
iia = b.resp.iia;
drmaxisdpsi = b.resp.drmaxisdpsi;
dzmaxisdpsi = b.resp.dzmaxisdpsi;
wa = b.resp.wa;
rbdef = b.rbdef;
zbdef = b.zbdef;
psibry = b.psibry;
iib = b.resp.iib;
drbdefdpsi = dr*b.bd.drdpsi(1,:);
dzbdefdpsi = dz*b.bd.dzdpsi(1,:);
wb = b.resp.wb;
wbrr = b.resp.wbrr;
wbrz = b.resp.wbrz;
wbzz = b.resp.wbzz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plasma is true if a closed contour encloses a sufficient area
plasma = b.plasma & b.Atot/c.Ag > 4;

% Assignments that don't depend on plasma
dihdx = [zeros(nih,nic+niv) eye(nih,nih) zeros(nih,nsf+nsp)];
if ~isfield(c,'jhalo')
  mph = zeros(ngg,nih);
  mch = zeros(nic,nih);
  mvh = zeros(niv,nih);
  mhh = zeros(nih,nih);
  psizr_halo = zeros(nz,nr);
  dpsizrdxe = [mpc mpv mph zeros(ngg,nsf+nsp+1)]; % Values will change if plasma
  psizr_app = reshape(mpc*ic + mpv*iv, nz, nr); % externally applied flux
  % plasma influence added later to psic,psiv,psih,dpsicdx,dpsivdx,dpsihdx
  psic = mcc *ic + mcv *iv;
  psiv = mcv'*ic + mvv *iv;
  psih = mch'*ic + mvh'*iv;
  dpsicdxe = [mcc  mcv  zeros(nic,nsf+nsp+1)];
  dpsivdxe = [mcv' mvv  zeros(niv,nsf+nsp+1)];
  dpsihdxe = zeros(0,nx+1);
else % halo = 1 % under development
  Ahalo = c.Acell;             % halo area within cells
  ihalo = Ahalo.*c.jhalo;      % halo current distribution
  ihalo = ihalo/sum(ihalo(:)); % halo current within cells per ih
  mph = mgg*ihalo(:);  % mutual inductance: grid to scalar halo current
  mch = mpc'*ihalo(:); % mutual inductance: coils to scalar halo current
  mvh = mpv'*ihalo(:); % mutual inductance: vessel to scalar halo current
  mhh = ihalo(:)'*mph; % mutual inductance: scalar halo to scalar halo
  psizr_halo = reshape(mph*ih,nz,nr);
  dpsizrdxe = [mpc mpv mph zeros(ngg,nsf+nsp+1)]; % Values will change if plasma
  psizr_app = reshape(mpc*ic + mpv*iv, nz, nr); % externally applied flux
  % plasma influence added later to psic,psiv,psih,dpsicdxe,dpsivdxe,dpsihdxe
  psic = mcc *ic + mcv *iv + mch*ih;
  psiv = mcv'*ic + mvv *iv + mvh*ih;
  psih = mch'*ic + mvh'*iv + mhh*ih;
  dpsicdxe = [mcc  mcv  mch zeros(nic,nsf+nsp+1)];
  dpsivdxe = [mcv' mvv  mvh zeros(niv,nsf+nsp+1)];
  dpsihdxe = [mch' mvh' mhh zeros(nih,nsf+nsp+1)];
end
y = zeros(ny,1); % Will hold values with new equilibrium = d.y0 + d.C*x
xdot0 = zeros(nx,1);
y0 = zeros(ny,1);
C = zeros(ny,nx);
D = zeros(ny,nu);
E = zeros(ny,1);
fl.test = zeros(nfltest,nx);
fl.test(1,ix.ih) = 1e-5;
fl.x = x;
bd.test = zeros(nbdtest,nx);
bd.x = x;

% equilibrium structure, e, with vacuum solution
nbbbs = nb+1;
a.x = 'State vector = [ic;iv;ih;sf;sp]';
e.x = x;
a.ci = 'coil circuit currents = x(ix,ic)';
e.ci = ic;
a.vi = 'vessel circuit currents = x(ix,iv)';
e.vi = iv;
e.x = x;
a.rg = 'R for grid coordinates';
e.rg = rg;
a.zg = 'Z for grid coordinates';
e.zg = zg;
a.psizr = 'Total flux on grid';
e.psizr = psizr;
a.psizr_app = 'Flux on grid from external conductors';
e.psizr_app = psizr_app;
a.psizr_pla = 'Flux on grid from plasma current';
e.psizr_pla = zeros(nz,nr);
a.psizr_halo = 'Flux on grid from halo current';
e.psizr_halo = zeros(nz,nr);
a.psizr_err = 'Numerical error in flux on grid';
e.psizr_err = zeros(nz,nr);
a.psimag = 'Flux at the magnetic axis';
e.psimag = psimag;
a.psibry = 'Flux on the boundary';
e.psibry = psibry;
a.psipla = 'Plasma flux = Integral[j*psi]/Ip';
e.psipla = nan;
a.rbbbs = 'R for (rather evenly spaced) boundary coordinates';
e.rbbbs = nan(nbbbs,1);
a.zbbbs = 'Z for (rather evenly spaced) boundary coordinates';
e.zbbbs = nan(nbbbs,1);
a.nbbbs = 'number of boundary coordinates';
e.nbbbs = nbbbs;
a.rmaxis = 'R for position of magnetic axis';
e.rmaxis = rmaxis;
a.zmaxis = 'Z for position of magnetic axis';
e.zmaxis = zmaxis;
a.rbdef = 'R for point that defines the boundary';
e.rbdef = rbdef;
a.zbdef = 'Z for point that defines the boundary';
e.zbdef = zbdef;
a.psibar = 'Normalized poloidal flux from axis to boundary';
e.psibar = linspace(0,1,nr);
a.rhot = 'sqrt(normalized toroidal flux) from axis to boundary';
e.rhot = nan(1,nr);
a.pres = 'Pressure from axis to boundary';
e.pres = zeros(1,nr)+c.pb+c.vp*sp;
a.fpol = 'fpol = 2e-7*(enclosed poloidal current) from axis to boundary';
e.fpol = zeros(1,nr)+c.fb+c.vf*sf;
a.pprime = 'pprime from axis to boundary';
e.pprime = zeros(1,nr);
a.ffprim = 'ffprim from axis to boundary';
e.ffprim = zeros(1,nr);
a.fprime = 'fprime from axis to boundary';
e.fprime = zeros(1,nr);
a.jpar = 'surface-averaged parallel current density from axis to boundary';
e.jpar = zeros(1,nr);
a.q = 'Safety factor, q from axis to boundary';
e.q = zeros(1,nr);
a.psic = 'Flux at coil circuits';
e.psic = psic;
a.psiv = 'Flux at vessel circuits';
e.psiv = psiv;
a.R8 = 'R for points at Rmax,~45,Zmax,~135,Rmin,~225,Zmin,~315';
e.R8 = R8;
a.Z8 = 'Z for points at Rmax,~45,Zmax,~135,Rmin,~225,Zmin,~315';
e.Z8 = Z8;
a.rsurf = 'R for geometric center';
e.rsurf = rsurf;
a.zsurf = 'Z for geometric center';
e.zsurf = zsurf;
a.aminor = 'Half the radial width of the plasma';
e.aminor = aminor;
a.bminor = 'Half the vertical height of the plasma';
e.bminor = bminor;
a.elong = 'Elongation = bminor/aminor';
e.elong = elong;
a.triu = 'Upper triangularity';
e.triu = triu;
a.tril = 'Lower triangularity';
e.tril = tril;
a.squo = 'Upper outer squareness';
e.squo = squo;
a.squi = 'Upper inner squareness';
e.squi = squi;
a.sqli = 'Lower inner squareness';
e.sqli = sqli;
a.sqlo = 'Lower outer squareness';
e.sqlo = sqlo;
a.drsep = 'flux difference to competing bdef/outboard midplane gradient';
e.drsep = drsep;
a.Cl = 'Contour length = length of the boundary';
e.Cl = Cl;
a.Vtot = 'Integral(2*pi*R*dA) = Total volume';
e.Vtot = Vtot;
a.Atot = 'Integral(dA) = Total area of cross section';
e.Atot = Atot;
a.Ltot = 'Integral(dA/R)';
e.Ltot = Ltot;
a.btsurf = 'Vacuum toroidal field at rsurf = fpol(end)/rsurf';
e.btsurf = nan;
a.Wth = 'Total thermal energy';
e.Wth = 0;
a.betap = 'Poloidal beta';
e.betap = 0;
a.betan = 'Normalized beta';
e.betan = 0;
a.pcurrt = 'Toroidal current within grid cells (grid point at the center)';
e.pcurrt = zeros(nz,nr);
a.rcur = 'R for current centroid';
e.rcur = nan;
a.zcur = 'Z for current centroid';
e.zcur = nan;
a.cpasma = 'Total toroidal plasma current';
e.cpasma = 0;
a.li = 'Normalized plasma inductance';
e.li = nan;
a.y2 = 'Second moment of current distribution';
e.y2 = nan;
a.y2n = 'Normalized second moment of current distibution';
e.y2n = nan;
e.info = a;

% Put current outside the last closed flux surface if its minor radius is too small
if b.aminor < c.amin
  if b.aminor == 0
    % Search for place to create a plasma with b.aminor = c.amin
  end
end

if plasma

  % flag grid cells containing plasma
  iplasma = igg | egg;
  
  % Flux difference between boundary and axis
  dba = psibry-psimag;

  % Normalized flux at grid points
  psibarzr = (psizr-psimag)/dba;

  % Normalized flux at current centers within cells, ync
  ync = (yc-psimag)/dba;

  % Find which spline knot to use for each point on the grid
  ikf = ones(nz,nr);
  for i = 2:nkf
    ikf(ync > psikf(i)) = i; % knots on grid for fpol
  end
  ikp = ones(nz,nr);
  for i = 2:nkp
    ikp(ync > psikp(i)) = i; % knots on grid for pres
  end

  % Coefficients for all spline knots
  fb = c.fb + c.vf*sf; % (R*Btvac)^2/2
  f0 = c.mf0*sf;
  f1 = c.mf1*sf;
  f2 = c.mf2*sf;
  f3 = c.mf3*sf;
  pb = c.pb + c.vp*sp; % pressure at boundary, typically 0
  p0 = c.mp0*sp;
  p1 = c.mp1*sp;
  p2 = c.mp2*sp;
  p3 = c.mp3*sp;

  % pres, fpol^2/2 depend on g
  g = dba/(2*pi);

  % Pressure at cell current centers, rgc, zgc
  preszr = g*(p0(ikp)+p1(ikp).*ync+p2(ikp).*ync.^2+p3(ikp).*ync.^3) + pb;
  preszr(~iplasma) = 0;

  % Thermal energy within cells
  Wzr = 3*pi*RAc.*preszr;

  % Volume integrals of Bpol^2 for cells
  B2Vc = ARc.*(ycr.^2 + ycz.^2)/(2*pi);

  % fpol^2/2 at cell current centers, rgc, zgc
  hfp2zr = mu0*g*(f0(ikf)+f1(ikf).*ync+f2(ikf).*ync.^2+f3(ikf).*ync.^3) + fb^2/2;

  % pcurrt = Integral[(R*pprime+ffprim/mu0/R)*dAcell]
  % Ignore R,Z dependency in pprime, ffprim, use values at rgc, zgc
  ppzr = p1(ikp) + 2*p2(ikp).*ync + 3*p3(ikp).*ync.^2; 
  fmzr = f1(ikf) + 2*f2(ikf).*ync + 3*f3(ikf).*ync.^2;
  pcurrt = RAc.*ppzr + ARc.*fmzr;

  % psizr_pla = Integral(m(R,Z)*i(R,Z))*dA
  % Ignore R, Z dependencies for m within cells
  psizr_pla = reshape(mgg*pcurrt(:),nz,nr);

  % Convergence error, psizr_err
  psizr_err = psizr-psizr_app-psizr_pla-psizr_halo;

  % Adding plasma contribution to fluxes
  psic = psic + mpc'*pcurrt(:); % Coils
  psiv = psiv + mpv'*pcurrt(:); % Vessel
  psih = psih + mph'*pcurrt(:); % Halo

  % Scalar parameters
  cpasma = sum(pcurrt(:));  % Total plasma current (inherited name from EFIT)
  btsurf = fb/rsurf;        % Vacuum toroidal B at rsurf
  Wth = sum(Wzr(:));        % Total thermal energy
  B2V = sum(B2Vc(:));       % Volume-integral of poloidal B
  bp2flx = (mu0*cpasma/Cl)^2;
  psipla = yc(:)'*pcurrt(:)/cpasma;         % Plasma flux
  betap = 4/3*mu0*(Wth/Vtot-1.5*pb)/bp2flx; % Poloidal beta
  ipnorm = abs(cpasma/1e6/aminor/btsurf);
  betat = 400/3*Wth/Vtot/btsurf^2*mu0;
  betan = betat/ipnorm;
  li = B2V/Vtot/bp2flx; % Volume-averaged Bp^2 / bp2flx
  rcur = sum(rgc(:).*pcurrt(:))/cpasma; % R of current centroid
  zcur = sum(zgc(:).*pcurrt(:))/cpasma; % Z of current centroid
  
  % Weighting function for second moment of current distribution, y2
  f2zr = ((rgc-rcur)/(2*rcur)+1).^2.*(rgc-rcur).^2 - ((rgc-rcur)/rcur+1).^2.*(zgc-zcur).^2;
  % Fusion Engineering and Design 43 (1998) 37-58 (equation 50 on page 50)
  
  % The second moment of the current distribution (equation 49 on page 50)
  y2 = f2zr(:)'*pcurrt(:)/cpasma;
  
  % The normalized second moment of the current distribution (equation 52 on page 50)
  y2n = y2*Atot/(f2zr(:)'*Ac(:));


  %%%%%%%%%%%%
  % Profiles %
  %%%%%%%%%%%%
  
  % Normalized flux for profiles
  s = linspace(0,1,np)';
  
  % Create nbbbs boundary points Rb, Zb
  cs = [1:nbbbs]'*n/nbbbs; % Floating index for b.R0, etc, for making rbbbs,zbbbs
  Cc = [0; cumsum(b.dC)]; % Total boundary length from point 1
  % Find floating indices that make boundary lengths between points roughly equal
  j = 2;
  for i = 1:nbbbs
    ct = Cc(n)*(i-1)/(nbbbs-1); % desired boundary length from point 1 to point cs(i)
    while Cc(j) < ct & j < n % second test is redundant, j never exceeds n anyway
      j = j+1;
    end
    % Now Cc(j) >= ct
    if Cc(j) == Cc(j-1)
      cs(i) = j;
    else
      cs(i) = j-1+(ct-Cc(j-1))/(Cc(j)-Cc(j-1));
    end
  end
  ics = floor(cs);
  xcs = cs-ics;
  Rb = b.R0(ics) + b.R1(ics).*xcs + b.R2(ics).*xcs.^2 + b.R3(ics).*xcs.^3;
  Rb(nbbbs) = Rb(1);
  Zb = b.Z0(ics) + b.Z1(ics).*xcs + b.Z2(ics).*xcs.^2 + b.Z3(ics).*xcs.^3;
  Zb(nbbbs) = Zb(1);
  
  % Index spline knots for s
  jkf = ones(np,1);
  for i = 2:nkf
    jkf(s > psikf(i)) = i; % knots on grid for fpol
  end
  jkp = ones(np,1);
  for i = 2:nkp
    jkp(s > psikp(i)) = i; % knots on grid for pres
  end

  % fpol^2/2 profile
  hf2 = mu0*g*[f0(jkf)+f1(jkf).*s+f2(jkf).*s.^2+f3(jkf).*s.^3]' + fb^2/2;

  % fpol profile
  fpol = sign(fb)*sqrt(2*hf2);

  % ffprim profile
  ffprim = mu0*[f1(jkf) + 2*f2(jkf).*s + 3*f3(jkf).*s.^2]';

  % fprime profile
  fprime = ffprim./fpol;
  
  % pressure profile
  pres = g*[p0(jkp)+p1(jkp).*s+p2(jkp).*s.^2+p3(jkp).*s.^3]' + pb;
  pres(end) = pb;

  % pprime profile
  pprime = [p1(jkp) + 2*p2(jkp).*s + 3*p3(jkp).*s.^2]';
  
  % Update equilibrium with plasma solution
  e.psizr_halo = psizr_halo;
  e.psizr_pla = psizr_pla;
  e.psizr_err = psizr_err;
  e.psipla = psipla;
  e.rbbbs = Rb;
  e.zbbbs = Zb;
  e.pres = pres;
  e.fpol = fpol;
  e.pprime = pprime;
  e.ffprim = ffprim;
  e.fprime = fprime;
  e.psic = psic;
  e.psiv = psiv;
  e.btsurf = btsurf;
  e.Wth = Wth;
  e.psipla = psipla;
  e.betap = betap;
  e.betan = betan;
  e.pcurrt = pcurrt;
  e.cpasma = cpasma;
  e.li = li;
  e.rcur = rcur;
  e.zcur = zcur;  
  e.y2 = y2;  
  e.y2n = y2n;  
  
  % Calculate profiles
  try
  p = gsprofiles(e,b);
  catch
    save ptrouble c e x psizr
  end
  
  rb = p.R0(:,end);
  zb = p.Z0(:,end);
  
  % Update equilibrium with more profiles
  e.rhot = p.rhot;
  e.jpar = p.jpar;
  e.q = p.q;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculation of plasma response, dpsizrdxe %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % psp = d(preszr)/d(sp)
  psp = g                       *c.mp0(ikp,:) + ...
        g*ync(:)   *ones(1,nsp).*c.mp1(ikp,:) + ...
        g*ync(:).^2*ones(1,nsp).*c.mp2(ikp,:) + ...
        g*ync(:).^3*ones(1,nsp).*c.mp3(ikp,:) + ...
	ones(ngg,1)*c.vp;

  % igsf = d(pcurrt)/d(sf)
  igsf = ARc(:)           *ones(1,nsf).*c.mf1(ikf,:) + ...
       2*ARc(:).*ync(:)   *ones(1,nsf).*c.mf2(ikf,:) + ...
       3*ARc(:).*ync(:).^2*ones(1,nsf).*c.mf3(ikf,:);

  % igsp = d(pcurrt)/d(sp)
  igsp = RAc(:)           *ones(1,nsp).*c.mp1(ikp,:) + ...
       2*RAc(:).*ync(:)   *ones(1,nsp).*c.mp2(ikp,:) + ...
       3*RAc(:).*ync(:).^2*ones(1,nsp).*c.mp3(ikp,:);

  % mutuals between grid and sf, sp
  mpf = mgg*igsf; % d(psizr)/d(sf)
  mpp = mgg*igsp; % d(psizr)/d(sp)

  % Response of psizr_pla to psizr (and psimag, psibry separately)
  preszrw = zeros(ngg,49);
  preszra = zeros(ngg,1);
  preszrb = zeros(ngg,1);
  Wzrw = zeros(ngg,49);
  Wzra = zeros(ngg,1);
  Wzrb = zeros(ngg,1);
  B2Vcw = zeros(ngg,49);
  B2Vcb = zeros(ngg,1);
  pcurrtw = zeros(ngg,49);
  pcurrta = zeros(ngg,1);
  pcurrtb = zeros(ngg,1);
  jg = zeros(ngg);   % d(dpsizr_pla)/dpsizr
  ja = zeros(ngg,1); % d(dpsizr_pla)/dpsimag
  jb = zeros(ngg,1); % d(dpsizr_pla)/dpsibry
  je = psizr_err(:);
  for i = 1:nz
    for j = 1:nr
      k = i+nz*(j-1);
      ok = k+i49>1 & k+i49<=ngg; % okay indices
      if iplasma(k)
      
	% Spline knots for grid cell k
	h = ikp(k);
	m = ikf(k);

	% Derivatives of pres/g, pprime, ffprim/mu0 w.r.t. ync
	pd = p1(h) + 2*p2(h)*ync(k) + 3*p3(h)*ync(k)^2;
        pp = 2*p2(h) + 6*p3(h)*ync(k);
        fm = 2*f2(m) + 6*f3(m)*ync(k);

	% yncw = d(ync(k))/d(psizr(k+i49))
	yncw = ycw(k,:)/dba;
	
	% How current ig in grid cell k changes when cell coverage changes because psizr changes
        igw = ppzr(k)*RAcw(k,:) + fmzr(k)*ARcw(k,:); % part of d(ig)/d(psizr)

	% Add how current ig changes when current density changes due to psizr change
        igw = igw + RAc(k)*pp*yncw + ARc(k)*fm*yncw; % d(ig)/d(psizr)

	% Response of current within cell k to change in psizr(k+i49)
	pcurrtw(k,:) = igw;

	% Response of plasma flux on the whole grid to changing ig due to psizr change
	jg(:,k+i49(ok)) = jg(:,k+i49(ok)) + mgg(:,k)*igw(ok);
	
        % Response of pressure at rgc(k), zgg(k) to change in psizr(k+i49)
	preszrw(k,:) = g*pd*yncw;
	
	% Response of thermal energy within cell k to change in psizr(k+i49)
        Wzrw(k,:) = 3*pi*RAc(k)*preszrw(k,:) + Wzr(k)/RAc(k)*RAcw(k,:);
        
	% Response of cell integral of Bpol^2 to change in psizr(k+i49)
	B2Vcw(k,:) = ARc(k)/pi*(ycr(k)*ycrw(k,:)+ycz(k)*yczw(k,:)) + B2Vc(k)/ARc(k)*ARcw(k,:);
	
	% ynca = d(ync(k))/d(psimag)
	ynca = (ync(k)-1)/dba;
	
	% How current ig changes when current density changes due to psimag change
	iga = RAc(k)*pp*ynca + ARc(k)*fm*ynca;

	% Response of current within cell k to change in psimag
	pcurrta(k) = iga;

        % Response of plasma flux on the whole grid to changing ig due to psimag change
	ja = ja + mgg(:,k)*iga;
		
        % Response of pressure at rgc(k), zgg(k) to change in psimag
        preszra(k) = g*pd*ynca - preszr(k)/dba;
        
	% Response of thermal energy within cell k to change in psimag
	Wzra(k) = 3*pi*RAc(k)*preszra(k);
	
	% yncb = d(ync(k))/d(psibry)
	yncb = (ycb(k)-ync(k))/dba;
	
	% How current ig changes when cell coverage changes due to psibry change	
	igb = RAc(k)*pp*yncb + ARc(k)*fm*yncb;

	% Add how current ig changes when current density changes due to psibry change
        igb = igb + ppzr(k)*RAcb(k) + fmzr(k)*ARcb(k);

	% Response of current within cell k to change in psibry
	pcurrtb(k) = igb;

        % Response of plasma flux on the whole grid to changing ig due to psibry change
	jb = jb + mgg(:,k)*igb;	
		
        % Response of pressure at rgc(k), zgg(k) to change in psibry
        preszrb(k) = g*pd*yncb + preszr(k)/dba;
        
	% Response of thermal energy within cell k to change in psibry
        Wzrb(k) = 3*pi*RAc(k)*preszrb(k) + Wzr(k)/RAc(k)*RAcb(k);
        
	% Response of cell integral of Bpol^2 to change in psibry
        B2Vcb(k) = ARc(k)/pi*(ycr(k)*ycrb(k)+ycz(k)*yczb(k)) + B2Vc(k)/ARc(k)*ARcb(k);

      end
    end
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Solve the matrix (eye(ngg)-jg) %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  % Solve for change of psizr_pla when x, e, b, a change
  % x = [ic;iv;sf;sp], e = fraction of flux error, b = psibry, a = psimag
  % dpsi = dpsipla+dpsiext = jg*dpsi+dpsiext => (E-jg)*dpsi = dpsiext
  dpsizrdxeba = (eye(ngg)-jg)\[mpc mpv mph mpf mpp je jb ja];
  
  % Alternative way to calculate dpsizrdxeba for use in Simulink matlab function blocks
  % dpsizrdxeba = axbop((eye(ngg)-jg),[mpc mpv mpf mpp je jb ja]);

  % Remove last column of dpsizrdxeba to obtain response to x, e, b
  va = wa*dpsizrdxeba(iia,:); % da = va*[dxeb; da]
  dadxeb = va(1:end-1)/(1-va(end)); % da*(1-va(end)) = va*dxeb
  dpsizrdxeb = dpsizrdxeba(:,1:end-1)+dpsizrdxeba(:,end)*dadxeb;

  % Remove last column of xeb to obtain response to x, e
  vb = wb*dpsizrdxeb(iib,:); % db = vb*[dxe; db]
  dbdxe = vb(1:end-1)/(1-vb(end)); % db*(1-vb(end)) = vb*dxe
  dpsizrdxe = dpsizrdxeb(:,1:end-1) + dpsizrdxeb(:,end)*dbdxe;
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Derive responses from dpsizrdxe %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  % Response of psimag, psibry
  dpsimagdxe = wa*dpsizrdxe(iia,:);
  dpsibrydxe = wb*dpsizrdxe(iib,:);

  % Response of psibarzr
  dpsibarzrdxe = (dpsizrdxe-(1-psibarzr(:))*dpsimagdxe-psibarzr(:)*dpsibrydxe)/dba;

  % Response of preszr, Wzr, B2V, pcurrt
  dpreszrdxe = zeros(ngg,nx+1);
  dWzrdxe = zeros(ngg,nx+1);
  dB2Vcdxe = zeros(ngg,nx+1);
  dpcurrtdxe = zeros(ngg,nx+1);
  drgcdxe = zeros(ngg,nx+1);
  dzgcdxe = zeros(ngg,nx+1);

  for i = 1:nz
    for j = 1:nr
      k = i+nz*(j-1);
      ok = k+i49>1 & k+i49<=ngg; % okay indices
      if egg(k)      
        drgcdxe(k,:) = rgcw(k,ok)*dpsizrdxe(k+i49(ok),:)+rgcb(k)*dpsibrydxe;
        dzgcdxe(k,:) = zgcw(k,ok)*dpsizrdxe(k+i49(ok),:)+zgcb(k)*dpsibrydxe;
      end
      if iplasma(k)
        dpreszrdxe(k,:) = preszrw(k,ok)*dpsizrdxe(k+i49(ok),:) + ...
	  preszra(k)*dpsimagdxe + preszrb(k)*dpsibrydxe;
	dpreszrdxe(k,ix.sp) = dpreszrdxe(k,ix.sp) + psp(k,:);
        dWzrdxe(k,:) = Wzrw(k,ok)*dpsizrdxe(k+i49(ok),:) + ...
	  Wzra(k)*dpsimagdxe + Wzrb(k)*dpsibrydxe;
	dWzrdxe(k,ix.sp) = dWzrdxe(k,ix.sp) + 3*pi*RAc(k)*psp(k,:);
        dB2Vcdxe(k,:) = B2Vcw(k,ok)*dpsizrdxe(k+i49(ok),:) + B2Vcb(k)*dpsibrydxe;
	dpcurrtdxe(k,:) = pcurrtw(k,ok)*dpsizrdxe(k+i49(ok),:) + ...
	  pcurrta(k)*dpsimagdxe + pcurrtb(k)*dpsibrydxe;
        dpcurrtdxe(k,ix.sp) = dpcurrtdxe(k,ix.sp) + igsp(k,:);
        dpcurrtdxe(k,ix.sf) = dpcurrtdxe(k,ix.sf) + igsf(k,:);
      end
    end
  end

  % Adding plasma influence on response of conductor circuit fluxes
  dpsicdxe = dpsicdxe + mpc'*dpcurrtdxe;
  dpsivdxe = dpsivdxe + mpv'*dpcurrtdxe;

  % Response of scalar parameters
  drmaxisdxe = drmaxisdpsi*dpsizrdxe(iia,:);
  dzmaxisdxe = dzmaxisdpsi*dpsizrdxe(iia,:);
  drbdefdxe = drbdefdpsi*dpsizrdxe(iib,:);
  dzbdefdxe = dzbdefdpsi*dpsizrdxe(iib,:);
  dcpasmadxe = sum(dpcurrtdxe);
  drsurfdxe  = rsurfp*dpsizrdxe;
  dzsurfdxe  = zsurfp*dpsizrdxe;
  daminordxe = aminorp*dpsizrdxe;
  dbminordxe = bminorp*dpsizrdxe;
  delongdxe = elongp*dpsizrdxe;
  dtriudxe = triup*dpsizrdxe;
  dtrildxe = trilp*dpsizrdxe;
  dsquodxe = squop*dpsizrdxe;
  dsquidxe = squip*dpsizrdxe;
  dsqlidxe = sqlip*dpsizrdxe;
  dsqlodxe = sqlop*dpsizrdxe;
  dCldxe = Clp*dpsizrdxe;
  dVtotdxe = Vtotp*dpsizrdxe;
  dAtotdxe = Atotp*dpsizrdxe;
  dLtotdxe = Ltotp*dpsizrdxe;
  dbtsurfdxe = -btsurf/rsurf*drsurfdxe;
  dbtsurfdxe(ix.sf) = dbtsurfdxe(ix.sf) + c.vf/rsurf;
  dWthdxe = sum(dWzrdxe);
  dbp2flxdxe = 2*bp2flx*(dcpasmadxe/cpasma-dCldxe/Cl);
  dpsipladxe = ycb(:)'*pcurrt(:)/cpasma*dpsibrydxe + yc(:)'*dpcurrtdxe/cpasma - psipla/cpasma*dcpasmadxe;
  for i = 1:nz
    for j = 1:nr
      k = i+nz*(j-1);
      ok = k+i49>1 & k+i49<=ngg; % okay indices
      if egg(k)
        dpsipladxe = dpsipladxe + pcurrt(k)/cpasma*ycw(k,ok)*dpsizrdxe(k+i49(ok),:);
      elseif igg(k)
        dpsipladxe = dpsipladxe + pcurrt(k)/cpasma*dpsizrdxe(k,:);
      end
    end
  end
  %betap = 4/3*mu0*(Wth/Vtot-1.5*pb)/bp2flx;
  dbetapdxe = 4/3*mu0/bp2flx*(dWthdxe/Vtot-Wth/Vtot^2*dVtotdxe) - betap/bp2flx*dbp2flxdxe;
  dbetapdxe(ix.sp) = dbetapdxe(ix.sp)-2*mu0/bp2flx*c.vp;
  dipnormdxe = ipnorm*(dcpasmadxe/cpasma - daminordxe/aminor - dbtsurfdxe/btsurf);
  dbetatdxe = betat*(dWthdxe/Wth - dVtotdxe/Vtot - 2*dbtsurfdxe/btsurf);
  dbetandxe = betan*(dbetatdxe/betat - dipnormdxe/ipnorm);
  dB2Vdxe = sum(dB2Vcdxe);
  dlidxe = li*(dB2Vdxe/B2V - dVtotdxe/Vtot - dbp2flxdxe/bp2flx);
  drcurdxe = sum(rgc(:)/cpasma*ones(1,nx+1).*dpcurrtdxe + ...
    pcurrt(:)/cpasma*ones(1,nx+1).*drgcdxe) - rcur/cpasma*dcpasmadxe;
  dzcurdxe = sum(zgc(:)/cpasma*ones(1,nx+1).*dpcurrtdxe + ...
    pcurrt(:)/cpasma*ones(1,nx+1).*dzgcdxe) - zcur/cpasma*dcpasmadxe;
  
  % Response of weighting function for second moment of current distribution
  df2zrdrgc = rgc.*(rgc.^2-rcur^2-2*(zgc-zcur).^2)/rcur^2;  
  df2zrdzgc = -((rgc-rcur)/rcur+1).^2.*(zgc-zcur)*2;  
  df2zrdrcur = (rcur^4-rgc.^4+4*zcur^2*rgc.^2-8*zcur*rgc.^2.*zgc+4*rgc.^2.*zgc.^2)/(2*rcur^3);
  df2zrdzcur = ((rgc-rcur)/rcur+1).^2.*(zgc-zcur)*2;
  df2zrdxe = df2zrdrgc(:)*ones(1,nx+1).*drgcdxe + df2zrdzgc(:)*ones(1,nx+1).*dzgcdxe + ...
    df2zrdrcur(:)*drcurdxe + df2zrdzcur(:)*dzcurdxe;
  
  % Response of second moment of the current distribution
  dy2dxe = pcurrt(:)'*df2zrdxe/cpasma + f2zr(:)'*dpcurrtdxe/cpasma - y2/cpasma*dcpasmadxe;
  
  % Response of normalized second moment of the current distribution
  % y2n = y2*Atot/(f2zr(:)'*Ac(:));
  dum = f2zr(:)'*Ac(:);
  ddumdxe = Ac(:)'*df2zrdxe + f2zr(:)'*Acb*dpsibrydxe;
  for i = 4:nz-3
    for j = 4:nr-3
      k = i+nz*(j-1);
      if egg(k)
	ddumdxe = ddumdxe + f2zr(k)*Acw(k,:)*dpsizrdxe(k+i49,:);
      end
    end
  end
  dy2ndxe = y2n/y2*dy2dxe + y2n/Atot*dAtotdxe - y2n/dum*ddumdxe;
      
  dgdxe = (dpsibrydxe-dpsimagdxe)/(2*pi);
  
  % Profile responses

  dhf2dxe = mu0*[f0(jkp)+f1(jkp).*s+f2(jkp).*s.^2+f3(jkp).*s.^3]*dgdxe;
  dhf2dxe(:,ix.sf) = dhf2dxe(:,ix.sf) + fb*ones(nr,1)*c.vf + mu0*g*(...
    c.mf0(jkf,:) + ...
    c.mf1(jkf,:).*[s*ones(1,nsf)] + ...
    c.mf2(jkf,:).*[s.^2*ones(1,nsf)] + ...
    c.mf3(jkf,:).*[s.^3*ones(1,nsf)]);  
  dfpoldxe = sign(fb)./sqrt(2*hf2')*ones(1,nx+1).*dhf2dxe;

  dffprimdxe = zeros(np,nx+1);
  dffprimdxe(:,ix.sf) = mu0*[c.mf1(jkf,:)+c.mf2(jkf,:).*[2*s*ones(1,nsf)]+...
    c.mf3(jkf,:).*[3*s.^2*ones(1,nsf)]];
  
  dfprimedxe = dffprimdxe./[e.fpol'*ones(1,nx+1)]-[e.fprime'./e.fpol'*ones(1,nx+1)].*dfpoldxe;

  dhf2ds = mu0*g*[f1(jkf)+2*f2(jkf).*s+3*f3(jkf).*s.^2]';
  dfpolds = sign(fb)./sqrt(2*hf2).*dhf2ds;
  dffprimds = mu0*[2*f2(jkf) + 6*f3(jkf).*s]';
  
  dffds = ones(nbbbs-1,1)*dffprimds;
  dfprimeds = dffprimds./e.fpol-e.fprime./e.fpol.*dfpolds;
  dfpds = ones(nbbbs-1,1)*dfprimeds;

  dpresdxe = [p0(jkp)+p1(jkp).*s+p2(jkp).*s.^2+p3(jkp).*s.^3]*dgdxe;
  dpresdxe(:,ix.sp) = dpresdxe(:,ix.sp) + ones(nr,1)*c.vp + g*(...
    c.mp0(jkp,:) + ...
    c.mp1(jkp,:).*[s*ones(1,nsp)] + ...
    c.mp2(jkp,:).*[s.^2*ones(1,nsp)] + ...
    c.mp3(jkp,:).*[s.^3*ones(1,nsp)]);  
  
  dpprimedxe = zeros(np,nx+1);
  dpprimedxe(:,ix.sp) = c.mp1(jkp,:)+c.mp2(jkp,:).*[2*s*ones(1,nsp)]+...
    c.mp3(jkp,:).*[3*s.^2*ones(1,nsp)];

  dpprimeds = [2*p2(jkp) + 6*p3(jkp).*s]';
  dppds = ones(nbbbs-1,1)*dpprimeds;
  
  dpsicontdxe = zeros((nbbbs-1)*nr,nx+1);
  k = 0;
  for j = 1:nx
    dpsicontdxe(:,j) = sum(p.resp.W.*dpsizrdxe(k+p.resp.Iw),2);
    k = k+ngg;
  end
  psibarcont = ones(nbbbs-1,1)*s';
%  psibarcont2 = reshape(sum(p.resp.W.*psibarzr(p.resp.Iw),2),nbbbs-1,nr);
  dpsibarcontdxe = (dpsicontdxe-(1-psibarcont(:))*dpsimagdxe-psibarcont(:)*dpsibrydxe)/dba;

  jdSxe = (p.jdSfp(:).*dfpds(:)+p.jdSpp(:).*dppds(:)+p.jdSff(:).*dffds(:))*ones(1,nx+1).*dpsibarcontdxe;
  
  % response of Integral(j*dS)/S along stationary paths
  djparFdxe = p.resp.jFps*dpsizrdxe; % changes in Yr, Yz changes fprime contribution
  k = 0;
  for i = 1:nr
    djparFdxe(i,:) = djparFdxe(i,:) + (...
      sum(jdSxe(k+[1:nbbbs-1],:)) + ... % changes in psibar changes all three
      p.jSfp(i)*dfprimedxe(i,:) + ...
      p.jSpp(i)*dpprimedxe(i,:) + ...
      p.jSff(i)*dffprimdxe(i,:))/p.S(i);
    k = k+nbbbs-1;
  end
  
  % Calculate drbdxe, dzbdxe
  Yg2 = p.Yr(:,end).^2+p.Yz(:,end).^2;
  gb = zeros(nb-1,nx+1);
  k = nb*(np-1);
  for i = 2:nb
    gb(i-1,:) =  sign(dba)/Yg2(i)*...
      (p.resp.W(k+i,:)*dpsizrdxe(p.resp.Iw(k+i,:),:)-dpsibrydxe);
  end
  drbdxe = [drbdefdxe; p.Yr(2:nb,end)*ones(1,nx+1).*gb];
  dzbdxe = [dzbdefdxe; p.Yz(2:nb,end)*ones(1,nx+1).*gb];
  
  
  % Extract d*dx and d*de from d*dxe

  dpsibarzrdx = dpsibarzrdxe(:,1:nx);
  dpsibarzrde = dpsibarzrdxe(:,nx+1);
  
  dpsibrydx = dpsibrydxe(:,1:nx);
  dpsibryde = dpsibrydxe(:,nx+1);
  
  dpsimagdx = dpsimagdxe(:,1:nx);
  dpsimagde = dpsimagdxe(:,nx+1);
  
  dpcurrtdx = dpcurrtdxe(:,1:nx);
  dpcurrtde = dpcurrtdxe(:,nx+1);  
  
  dcpasmadx = dcpasmadxe(:,1:nx);
  dcpasmade = dcpasmadxe(:,nx+1);  
  
  dpprimedx = dpprimedxe(:,1:nx);
  dpprimede = dpprimedxe(:,nx+1);  
  
  dffprimdx = dffprimdxe(:,1:nx);
  dffprimde = dffprimdxe(:,nx+1);  
  
  dfprimedx = dfprimedxe(:,1:nx);
  dfprimede = dfprimedxe(:,nx+1);  
  
  dpreszrdx = dpreszrdxe(:,1:nx);
  dpreszrde = dpreszrdxe(:,nx+1);  
  
  dpsizrdx = dpsizrdxe(:,1:nx);
  dpsizrde = dpsizrdxe(:,nx+1);
  
  djparFdx = djparFdxe(:,1:nx);
  djparFde = djparFdxe(:,nx+1);  
  
  dpresdx = dpresdxe(:,1:nx);
  dpresde = dpresdxe(:,nx+1);
  
  dfpoldx = dfpoldxe(:,1:nx);
  dfpolde = dfpoldxe(:,nx+1);
  
  dpsicdx = dpsicdxe(:,1:nx);
  dpsicde = dpsicdxe(:,nx+1);
  
  dpsivdx = dpsivdxe(:,1:nx);
  dpsivde = dpsivdxe(:,nx+1);
  
  dpsihdx = dpsihdxe(:,1:nx);
  dpsihde = dpsihdxe(:,nx+1);
  
  dB2Vcdx = dB2Vcdxe(:,1:nx);
  dB2Vcde = dB2Vcdxe(:,nx+1);  
  
  dWzrdx = dWzrdxe(:,1:nx);
  dWzrde = dWzrdxe(:,nx+1);  
  
  dlidx = dlidxe(:,1:nx);
  dlide = dlidxe(:,nx+1);
  
  dy2dx = dy2dxe(:,1:nx);
  dy2de = dy2dxe(:,nx+1);
  
  dy2ndx = dy2ndxe(:,1:nx);
  dy2nde = dy2ndxe(:,nx+1);

  % Shape
  
  drsurfdx = drsurfdxe(:,1:nx);
  drsurfde = drsurfdxe(:,nx+1);
  
  dzsurfdx = dzsurfdxe(:,1:nx);
  dzsurfde = dzsurfdxe(:,nx+1);
  
  daminordx = daminordxe(:,1:nx);
  daminorde = daminordxe(:,nx+1);
  
  dbminordx = dbminordxe(:,1:nx);
  dbminorde = dbminordxe(:,nx+1);
  
  delongdx = delongdxe(:,1:nx);
  delongde = delongdxe(:,nx+1);
  
  dtriudx = dtriudxe(:,1:nx);
  dtriude = dtriudxe(:,nx+1);
  
  dtrildx = dtrildxe(:,1:nx);
  dtrilde = dtrildxe(:,nx+1);
  
  dsquodx = dsquodxe(:,1:nx);
  dsquode = dsquodxe(:,nx+1);
  
  dsquidx = dsquidxe(:,1:nx);
  dsquide = dsquidxe(:,nx+1);
  
  dsqlidx = dsqlidxe(:,1:nx);
  dsqlide = dsqlidxe(:,nx+1);
  
  dsqlodx = dsqlodxe(:,1:nx);
  dsqlode = dsqlodxe(:,nx+1);
   
  dCldx = dCldxe(:,1:nx);
  dClde = dCldxe(:,nx+1);
 
  dVtotdx = dVtotdxe(:,1:nx);
  dVtotde = dVtotdxe(:,nx+1);
  
  dAtotdx = dAtotdxe(:,1:nx);
  dAtotde = dAtotdxe(:,nx+1);
  
  dLtotdx = dLtotdxe(:,1:nx);
  dLtotde = dLtotdxe(:,nx+1);
  
  dbtsurfdx = dbtsurfdxe(:,1:nx);
  dbtsurfde = dbtsurfdxe(:,nx+1);

  dWthdx = dWthdxe(:,1:nx);
  dWthde = dWthdxe(:,nx+1);

  dbp2flxdx = dbp2flxdxe(:,1:nx);
  dbp2flxde = dbp2flxdxe(:,nx+1);

  dpsipladx = dpsipladxe(:,1:nx);
  dpsiplade = dpsipladxe(:,nx+1);

  dbetapdx = dbetapdxe(:,1:nx);
  dbetapde = dbetapdxe(:,nx+1);

  dipnormdx = dipnormdxe(:,1:nx);
  dipnormde = dipnormdxe(:,nx+1);

  dbetatdx = dbetatdxe(:,1:nx);
  dbetatde = dbetatdxe(:,nx+1);

  dbetandx = dbetandxe(:,1:nx);
  dbetande = dbetandxe(:,nx+1);

  dB2Vdx = dB2Vdxe(:,1:nx);
  dB2Vde = dB2Vdxe(:,nx+1);

  drcurdx = drcurdxe(:,1:nx);
  drcurde = drcurdxe(:,nx+1);

  dzcurdx = dzcurdxe(:,1:nx);
  dzcurde = dzcurdxe(:,nx+1);

  drmaxisdx = drmaxisdxe(:,1:nx);
  drmaxisde = drmaxisdxe(:,nx+1);

  dzmaxisdx = dzmaxisdxe(:,1:nx);
  dzmaxisde = dzmaxisdxe(:,nx+1);

  drbdx = drbdxe(:,1:nx);
  drbde = drbdxe(:,nx+1);

  dzbdx = dzbdxe(:,1:nx);
  dzbde = dzbdxe(:,nx+1);

  drbdefdx = drbdefdxe(:,1:nx);
  drbdefde = drbdefdxe(:,nx+1);

  dzbdefdx = dzbdefdxe(:,1:nx);
  dzbdefde = dzbdefdxe(:,nx+1);

  %%%%%%%%%%%%
  % Dynamics %
  %%%%%%%%%%%%  
  
  % PONDERINGS ABOUT INERTIA (seems nothing can be done to include it)
  % A 1-parameter model for the effect of vertical inertia:
  % plasma mass, mp = nm*Vtot, where nm is mass density
  % Approximate center of mass with center of current, zcur
  % Vertical force accelerating the plasma, Fz = mp*zcurdotdot
  % And Fz = Integral[jphi*d(psi_app)/dz]
  % For a GS equilibrium Fz = 0, so only the perturbation produces a force
  % Assume a rigid vertical shift (dzm) of jphi produces the force
  % This means that psizr_pla also has the same vertical shift
  % Fz = Integral[jphi_z*dzm*psiapp_z] = mp*dzcurdx*xdotdot
  % dzm = mp*dzcurdx*xdotdot/Integral[jphi_z*psiapp_z]
  % dzm = fzi*xdotdot, where fzi = mp*dzcurdx/Integral[jphi_z*psiapp_z]  
  % Without intertia the total flux change is, dpsizr = dpsizrdx*dx
  % With intertia, dpsizr = dpsizrdx*dx + psizrpla_z*fzi*xdotdot
  
  % mssx = P'*Lext*P; % External inductance in circuit with coil (see gsconfig)
  % lstarx = [mssx zeros(nic+niv,nih+nsf+nsp); zeros(nih+nsf+nsp,nx)];
  lstar = lstarx+[dpsicdx; dpsivdx; dpsihdx; zeros(nsp+nsf,nic+niv+nih) eye(nsp+nsf)];
  if evolve_option == 0 % Mimic scalar Ip evolution done with linear models    
    lstari = inv(lstar);
    A = -lstari*Rhat;
    B = lstari*Vhat;
  elseif evolve_option == 1 % Mimic scalar Ip evolution done with linear models    
    lstar(nic+niv+nih+(1:3),:) = [dpsipladx; dlidx; dbetapdx];
    lstari = inv(lstar);
    A = -lstari*Rhat;
    B = lstari*Vhat;
  elseif evolve_option == 2
    % psis is surface-averaged flux = poloidal flux + toroidal flux/q
    dpsisdx = p.resp.psisps*dpsizrdx + p.resp.psisfp*dfpoldx;
    % Pj projects equations for psis on nsf equations
    lai = inv(lstar([ix.ic ix.iv ix.ih ix.sp],[ix.ic ix.iv ix.ih ix.sp]));
    dxdsf([ix.ic ix.iv ix.ih ix.sp ix.sf],1:nsf) = ...
      [-lai*lstar([ix.ic ix.iv ix.ih ix.sp],ix.sf); eye(nsf)];
    Pj = pinv(dpsisdx*dxdsf);
    lstar(ix.sf,:) = Pj*dpsisdx;
    lstari = inv(lstar);
    A = -lstari*Rhat;
    B = lstari*[eye(nic) zeros(nic,nih+np+nsp); ... % Coil power supplies
      zeros(niv,nu); ... % Vessel has zero applied voltage
      zeros(nih,nic) ones(nih) zeros(nih,np+nsp); ... % halo voltage
      zeros(nsf,nic+nih) Pj zeros(nsf,nsp); ... % voltage profile
      zeros(nsp,nic+nih+np) eye(nsp)];
  end
  [vecs,vals] = eigsort(A);
  gamma = vals(1);
  drcurdv = drcurdx*vecs(:,1);
  dzcurdv = dzcurdx*vecs(:,1);

  % Note that A may become complex by this
  % if gamma>100 & all(imag(vecs(:,1))==0)
  %	A = A - A*vecs(:,1)*vecs(:,1)';
  %	stabilized = true;
  % end
  % A good algorithm for artifical stability should be developed

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % fltest decides when to calculate new response %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  tt = 1+(0:nfla-1)'*(n-1)/nfla; % floating index into R0, Z0
  it = floor(tt); % index into R0, Z0
  xt = tt-it;
  rbfl = b.R0(it) + b.R1(it).*xt + b.R2(it).*xt.^2 + b.R3(it).*xt.^3;
  zbfl = b.Z0(it) + b.Z1(it).*xt + b.Z2(it).*xt.^2 + b.Z3(it).*xt.^3;
  rfl = zeros(nfla*nflc,1);
  zfl = zeros(nfla*nflc,1);
  for i = 0:nflc-1
    rfl(i*nfla+(1:nfla)) = ((nflc-i)*rbfl+i*rmaxis)/nflc;
    zfl(i*nfla+(1:nfla)) = ((nflc-i)*zbfl+i*zmaxis)/nflc;
  end
  rfl = (rfl-rg(1))/dr+1;
  zfl = (zfl-zg(1))/dz+1;
  ifl = floor(zfl);
  jfl = floor(rfl);
  fl.teste = zeros(nfla*nflc,1); % fl.test for removing flux error
  for i = 1:nfla*nflc
    tr = rfl(i)-jfl(i);
    wr0 = [1 tr tr^2 tr^3]*mx;
    tz = zfl(i)-ifl(i);
    wz0 = [1 tz tz^2 tz^3]*mx;    
    wb = reshape(wz0'*wr0,1,16);
    k = ifl(i)+nz*(jfl(i)-1);
    fl.test(i,:) = wb*dpsibarzrdx(k+i16,:)/dpsibarlimit;
    fl.teste(i) = -wb*dpsibarzrde(k+i16)/dpsibarlimit;
  end
  flmax = max(abs(fl.teste));
  % dpsizrde can be much larger than psizr_err,
  % in such cases only a fraction should be removed
  de = min(1,1/flmax); % fraction of flux error to remove

  % Check evolution speed
  % flspeed = max(abs(d.fl.test*A*x));
  % When gamma becomes infinite the calculated value is zero
  % but flspeed may still be high
  % flspeed can be slowed down with:
  % lstari2 = pinv([lstar;wn*d.fl.test]);
  % lstari = lstari2(1:nx,1:nx);
  % where suitable wn is found that makes flspeed*tstep < 1

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % bdtest detects jumps in boundary-defining point %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  dpsibd = nan(nbdtest,1);
  dpsibdde = nan(nbdtest,1);
  for i = 2:bd.n
    dpsibdde(i) = bd.w(i,:)*dpsizrde(bd.ii(i,:))-dpsibryde;
    dpsibd(i) = bd.psi(i)-psibry;
  end
  % Don't allow removal of flux error to cause jump of boundary-defining point
  bdmax = max(dpsibdde*de./dpsibd); % bdtest for flux error removal
  de = de/max(1,2*bdmax); % Ensure -dpsibdde*de./dpsibd <= 0.5
  dpsibdc = dpsibd - dpsibdde*de; % Flux error removed from dpsibd
  for i = 2:bd.n
    dpsibddx = bd.w(i,:)*dpsizrdx(bd.ii(i,:),:)-dpsibrydx;
    bd.test(i,:) = -dpsibddx/dpsibdc(i);
    bd.teste(i,:) = -dpsibdde(i)/dpsibd(i);
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Outputs, plasma specific %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  if isfield(iy,'rbdef')
    y(iy.rbdef) = rbdef;
    C(iy.rbdef,:) = drbdefdx;
    if bdmax1 > 1
      iE(iy.rbdef) = 0;
    end
  end
  if isfield(iy,'zbdef')
    y(iy.zbdef) = zbdef;
    C(iy.zbdef,:) = dzbdefdx;
    if bdmax1 > 1
      iE(iy.zbdef) = 0;
    end
  end
  if isfield(iy,'psibry')
    y(iy.psibry) = psibry;
    C(iy.psibry,:) = dpsibrydx;
  end
  if isfield(iy,'psimag')
    y(iy.psimag) = psimag;
    C(iy.psimag,:) = dpsimagdx;
  end
  if isfield(iy,'R8')
    y(iy.R8) = b.R8;
    C(iy.R8,:) = b.resp.R8p*dpsizrdx;
  end
  if isfield(iy,'Z8')
    y(iy.Z8) = b.Z8;
    C(iy.Z8,:) = b.resp.Z8p*dpsizrdx;
  end
  if isfield(iy,'rsurf')
    y(iy.rsurf) = b.rsurf;
    C(iy.rsurf,:) = b.resp.rsurfp*dpsizrdx;
  end
  if isfield(iy,'zsurf')
    y(iy.zsurf) = b.zsurf;
    C(iy.zsurf,:) = b.resp.zsurfp*dpsizrdx;
  end
  if isfield(iy,'aminor')
    y(iy.aminor) = b.aminor;
    C(iy.aminor,:) = b.resp.aminorp*dpsizrdx;
  end
  if isfield(iy,'bminor')
    y(iy.bminor) = b.bminor;
    C(iy.bminor,:) = b.resp.bminorp*dpsizrdx;
  end
  if isfield(iy,'elong')
    y(iy.elong) = b.elong;
    C(iy.elong,:) = b.resp.elongp*dpsizrdx;
  end
  if isfield(iy,'triu')
    y(iy.triu) = b.triu;
    C(iy.triu,:) = b.resp.triup*dpsizrdx;
  end
  if isfield(iy,'tril')
    y(iy.tril) = b.tril;
    C(iy.tril,:) = b.resp.trilp*dpsizrdx;
  end
  if isfield(iy,'squo')
    y(iy.squo) = b.squo;
    C(iy.squo,:) = b.resp.squop*dpsizrdx;
  end
  if isfield(iy,'squi')
    y(iy.squi) = b.squi;
    C(iy.squi,:) = b.resp.squip*dpsizrdx;
  end
  if isfield(iy,'sqli')
    y(iy.sqli) = b.sqli;
    C(iy.sqli,:) = b.resp.sqlip*dpsizrdx;
  end
  if isfield(iy,'sqlo')
    y(iy.sqlo) = b.sqlo;
    C(iy.sqlo,:) = b.resp.sqlop*dpsizrdx;
  end
  if isfield(iy,'drsep')
    y(iy.drsep) = b.drsep;
    C(iy.drsep,:) = b.resp.drsepp*dpsizrdx;
  end
  if isfield(iy,'Cl')
    y(iy.Cl) = b.Cl;
    C(iy.Cl,:) = b.resp.Clp*dpsizrdx;
  end
  if isfield(iy,'Vtot')
    y(iy.Vtot) = b.Vtot;
    dVtotdx = 2*pi*sum(b.resp.dRA.b(1:b.n-1))*dpsibrydx;
    for i = 1:b.n-1
      dVtotdx = dVtotdx + 2*pi*b.resp.dRA.w(i,:)*dpsizrdx(b.resp.I36(i,:),:);
    end
    C(iy.Vtot,:) = dVtotdx;
  end
  if isfield(iy,'Atot')
    y(iy.Atot) = b.Atot;
    dAtotdx = sum(b.resp.dA0.b(1:b.n-1))*dpsibrydx;
    for i = 1:b.n-1
      dAtotdx = dAtotdx + b.resp.dA0.w(i,:)*dpsizrdx(b.resp.I36(i,:),:);
    end
    C(iy.Atot,:) = dAtotdx;
  end
  if isfield(iy,'Ltot')
    y(iy.Ltot) = b.Ltot;
    dLtotdx = sum(b.resp.dAR.b(1:b.n-1))*dpsibrydx;
    for i = 1:b.n-1
      dLtotdx = dLtotdx + b.resp.dAR.w(i,:)*dpsizrdx(b.resp.I36(i,:),:);
    end
    C(iy.Ltot,:) = dLtotdx;
  end
  if isfield(iy,'rb')
    y(iy.rb) = rb;
    C(iy.rb,:) = drbdx;
  end
  if isfield(iy,'zb')
    y(iy.zb) = zb;
    C(iy.zb,:) = dzbdx;
  end
  if isfield(iy,'rcur')
    y(iy.rcur) = rcur;
    C(iy.rcur,:) = drcurdx;
  end
  if isfield(iy,'zcur')
    y(iy.zcur) = zcur;
    C(iy.zcur,:) = dzcurdx;
  end
  if isfield(iy,'rcurdot')
    C(iy.rcurdot,:) = drcurdx*A;
    D(iy.rcurdot,:) = drcurdx*B;
    y(iy.rcurdot) = C(iy.rcurdot,:)*x;
  end
  if isfield(iy,'zcurdot')
    C(iy.zcurdot,:) = dzcurdx*A;
    D(iy.zcurdot,:) = dzcurdx*B;
    y(iy.zcurdot) = C(iy.zcurdot,:)*x;
  end
  if isfield(iy,'fluxerror')
    y(iy.fluxerror) = max(abs(psizr_err(:)/dba));
  end
  if isfield(iy,'isdiverted')
    y(iy.isdiverted) = double(lim == 0);
  end
  if isfield(iy,'psihalo') & isfield(c,'jhalo')
    y(iy.psihalo) = psih;
    C(iy.psihalo,:) = dpsihdx;
  end
  if isfield(iy,'psipla')
    y(iy.psipla) = psipla;
    C(iy.psipla,:) = dpsipladx;
  end
  if isfield(iy,'dpsiplade')
    y(iy.dpsiplade) = dpsiplade;
  end
  if isfield(iy,'cpasma')
    y(iy.cpasma) = cpasma;
    C(iy.cpasma,:) = dcpasmadx;
  end
  if isfield(iy,'Ip')
    y(iy.Ip) = cpasma/1e6;
    C(iy.Ip,:) = dcpasmadx/1e6;
  end
  if isfield(iy,'li')
    y(iy.li) = li;
    C(iy.li,:) = dlidx;
  end
  if isfield(iy,'y2')
    y(iy.y2) = y2;
    C(iy.y2,:) = dy2dx;
  end
  if isfield(iy,'y2n')
    y(iy.y2n) = y2n;
    C(iy.y2n,:) = dy2ndx;
  end
  if isfield(iy,'betap')
    y(iy.betap) = betap;
    C(iy.betap,:) = dbetapdx;
  end
  if isfield(iy,'Wth')
    y(iy.Wth) = Wth;
    C(iy.Wth,:) = dWthdx;
  end
  if isfield(iy,'L1t')
    y(iy.L1t) = p.L1t;
  end
  if isfield(iy,'psit')
    y(iy.psit) = p.psit;
    C(iy.psit,:) = p.resp.psitfp*dfpoldx;
  end
  if isfield(iy,'qpsi')
    y(iy.qpsi) = p.q;
  end
  if isfield(iy,'vind')
    y(iy.vind) = 0;
    C(iy.vind,:) = p.resp.psisps*dpsizrdx*A + p.resp.psisfp*dfpoldx*A;
    D(iy.vind,:) = p.resp.psisps*dpsizrdx*B + p.resp.psisfp*dfpoldx*B;
  end
  if isfield(iy,'jpar')
    y(iy.jpar) = p.jpar;
    C(iy.jpar,:) = p.resp.jFps*dpsizrdx;
  end

  % Done with objects for plasma case

else % No plasma

  pcurrt = zeros(nz,nr);
  cpasma = 0;
  li = nan;
  betap = 0;
  rcur = nan;
  zcur = nan;
  y2 = nan;
  rb = nan(nb,1);
  zb = nan(nb,1);
  dpsizrdx = dpsizrdxe(:,1:nx);
  dpsimagdxe = zeros(1,nx+1);
  dpsibrydxe = zeros(1,nx+1);
  dpresdxe = zeros(np,nx+1);
  dpprimedxe = zeros(np,nx+1);
  dfpoldxe = zeros(np,nx+1);
  dffprimdxe = zeros(np,nx+1);
  dfprimedxe = zeros(np,nx+1);
  djparFdxe = zeros(np,nx+1);
  dpcurrtdxe = zeros(ngg,nx+1);
  dpreszrdxe = zeros(ngg,nx+1);
  dWzrdxe = zeros(ngg,nx+1);
  dB2Vcdxe = zeros(ngg,nx+1);
  dcpasmadxe = zeros(1,nx+1);
  drsurfdxe = zeros(1,nx+1);
  dzsurfdxe = zeros(1,nx+1);
  daminordxe = zeros(1,nx+1);
  dbminordxe = zeros(1,nx+1);
  delongdxe = zeros(1,nx+1);
  dtriudxe = zeros(1,nx+1);
  dtrildxe = zeros(1,nx+1);
  dsquodxe = zeros(1,nx+1);
  dsquidxe = zeros(1,nx+1);
  dsqlidxe = zeros(1,nx+1);
  dsqlodxe = zeros(1,nx+1);
  dCldxe = zeros(1,nx+1);
  dVtotdxe = zeros(1,nx+1);
  dAtotdxe = zeros(1,nx+1);
  dLtotdxe = zeros(1,nx+1);
  dbtsurfdxe = zeros(1,nx+1);
  dWthdxe = zeros(1,nx+1);
  dbp2flxdxe = zeros(1,nx+1);
  dpsipladxe = zeros(1,nx+1);
  dbetapdxe = zeros(1,nx+1);
  dbetatdxe = zeros(1,nx+1);
  dbetandxe = zeros(1,nx+1);
  dB2Vdxe = zeros(1,nx+1);
  dlidxe = zeros(1,nx+1);
  drcurdxe = zeros(1,nx+1);
  dzcurdxe = zeros(1,nx+1);
  dy2dxe = zeros(1,nx+1);
  dy2ndxe = zeros(1,nx+1);
  drmaxisdxe = zeros(1,nx+1);
  dzmaxisdxe = zeros(1,nx+1);
  drbdxe = zeros(nb,nx+1);
  dzbdxe = zeros(nb,nx+1);

  dpsicdx = dpsicdxe(:,1:nx);
  dpsivdx = dpsivdxe(:,1:nx);
  dpsihdx = dpsihdxe(:,1:nx);
  lstar = [dpsicdx; dpsivdx; dpsihdx; zeros(nsp+nsf,nic+niv+nih) eye(nsp+nsf)];
  lstari = inv(lstar);
  lstari(ix.sf,:) = lstari(ix.sf,:)*0.01;
  lstari(ix.sp,:) = lstari(ix.sp,:)*0.01;
  A = -lstari*Rhat;
  B = lstari*Vhat;
  [vecs,vals] = eigsort(A);            
  gamma = vals(1);
  drcurdv = drcurdx*vecs(:,1);
  dzcurdv = dzcurdx*vecs(:,1);
  de = 1;

  if isfield(iy,'isdiverted')
    y(iy.isdiverted) = -1;
  end
  
  % Vacuum solution for p
  p = gsprofiles(e,b);
  
  % Done with objects for no-plasma case

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs, plasma or no plasma %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

if isfield(iy,'x')
  y(iy.x) = x;
  C(iy.x,:) = eye(nx);
end
if isfield(iy,'psizr')
  y(iy.psizr) = psizr(:);
  C(iy.psizr,:) = dpsizrdx;
end
if isfield(iy,'pcurrt')
  y(iy.pcurrt) = pcurrt(:);
  C(iy.pcurrt,:) = dpcurrtdx;
end
if isfield(iy,'it') % mu0*it=2*pi*R*Bvac => it = 5e6*R*Bvac
  y(iy.it) = 5e6*(c.fb+c.vf*sf);
  C(iy.it,ix.sf) = 5e6*c.vf;
end
if isfield(iy,'ic')
  y(iy.ic) = ic;
  C(iy.ic,ix.ic) = eye(nic);
end
if isfield(iy,'iv')
  y(iy.iv) = iv;
  C(iy.iv,ix.iv) = eye(niv);
end
if isfield(iy,'ih')
  y(iy.ih) = ih;
  C(iy.ih,ix.ih) = eye(nih);
end
if isfield(iy,'sp')
  nsp = numel(iy.sp);
  y(iy.sp) = sp;
  C(iy.sp,ix.sp) = eye(nsp);
end
if isfield(iy,'sf')
  y(iy.sf) = sf;
  C(iy.sf,ix.sf) = eye(nsf);
end
if isfield(iy,'fl') % flux loop signals
  C(iy.fl,:) = c.mpl'*dpcurrtdx;
  C(iy.fl,ix.ic) = C(iy.fl,ix.ic) + c.mlc;
  C(iy.fl,ix.iv) = C(iy.fl,ix.iv) + c.mlv;
  y(iy.fl) = c.mlc*ic+c.mlv*iv+c.mpl'*pcurrt(:);
end
if isfield(iy,'lv') % Loop voltage signals
  C(iy.lv,:) = -c.mph'*dpcurrtdx*A - [c.mhc c.mhv]*A(1:nic+niv,:);
  D(iy.lv,:) = -c.mph'*dpcurrtdx*B - [c.mhc c.mhv]*B(1:nic+niv,:);
  y(iy.lv) = C(iy.lv,:)*x;
end
if isfield(iy,'bp')
  C(iy.bp,:) = c.gpb'*dpcurrtdx;
  C(iy.bp,ix.ic) = C(iy.bp,ix.ic) + c.gbc;
  C(iy.bp,ix.iv) = C(iy.bp,ix.iv) + c.gbv;
  y(iy.bp) = c.gbc*ic+c.gbv*iv+c.gpb'*pcurrt(:);
end
if isfield(iy,'rog')
  nc = size(c.Pcc,1);
  nv = size(c.Pvc,1);
  C(iy.rog,:) = c.rldata(:,end)*dcpasmadx;
  C(iy.rog,ix.ic) = C(iy.rog,ix.ic) + c.rldata(:,1:nc)*c.Pcc;
  C(iy.rog,ix.iv) = C(iy.rog,ix.iv) + c.rldata(:,nc+(1:nv))*c.Pvc;
  y(iy.rog) = c.rldata(:,1:nc)*c.Pcc*ic + ...
              c.rldata(:,nc+(1:nv))*c.Pvc*iv + ...
	      c.rldata(:,end)*cpasma;
end
if isfield(iy,'gamma')
  y(iy.gamma) = real(gamma);
end
if isfield(iy,'drcurdv')
  y(iy.drcurdv) = real(drcurdv);
end
if isfield(iy,'dzcurdv')
  y(iy.dzcurdv) = real(dzcurdv);
end
if isfield(iy,'nupdate')
  y(iy.nupdate) = nupdate;
end
if isfield(iy,'ncalc')
  y(iy.ncalc) = ncalc;
end

% Calculate y0 so that y = y0 + C*x
y0 = y-C*x;

% Calculate E to keep old signals initially and morph to new
% y + E = yold
if plasma
  E(iE) = yold(iE) - y(iE);
end

% Archive dynamics
d.A = A;
d.B = B;
d.C = C;
d.D = D;
d.E = E;
d.xdot0 = xdot0;
d.y0 = y0;
d.gamma = real(gamma);
d.fl = fl;
d.bd = bd;
d.de = de;
d.fmax = 1;

% Archive response


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response to state vector x

r.dpsizrdx = dpsizrdxe(:,1:nx);
r.dpsimagdx = dpsimagdxe(1:nx);
r.dpsibrydx = dpsibrydxe(1:nx);
r.dpsicdx = dpsicdxe(:,1:nx);
r.dpsivdx = dpsivdxe(:,1:nx);

r.dpresdx = dpresdxe(:,1:nx);
r.dpprimedx = dpprimedxe(:,1:nx);
r.dfpoldx = dfpoldxe(:,1:nx);
r.dffprimdx = dffprimdxe(:,1:nx);
r.dfprimedx = dfprimedxe(:,1:nx);
r.djparFdx = djparFdxe(:,1:nx);
r.dpcurrtdx = dpcurrtdxe(:,1:nx);

r.dcpasmadx = dcpasmadxe(1:nx);
r.dWthdx = dWthdxe(1:nx);
r.dlidx = dlidxe(1:nx);
r.dbetapdx = dbetapdxe(1:nx);
r.dbetatdx = dbetatdxe(1:nx);
r.dbetandx = dbetandxe(1:nx);
r.dpsipladx = dpsipladxe(1:nx);
r.dbtsurfdx = dbtsurfdxe(1:nx);
r.drcurdx = drcurdxe(1:nx);
r.dzcurdx = dzcurdxe(1:nx);
r.dy2dx = dy2dxe(1:nx);
r.dy2ndx = dy2ndxe(1:nx);
r.drmaxisdx = drmaxisdxe(1:nx);
r.dzmaxisdx = dzmaxisdxe(1:nx);

r.drsurfdx = drsurfdxe(1:nx);
r.dzsurfdx = dzsurfdxe(1:nx);
r.daminordx = daminordxe(1:nx);
r.dbminordx = dbminordxe(1:nx);
r.delongdx = delongdxe(1:nx);

r.dtriudx = dtriudxe(1:nx);
r.dtrildx = dtrildxe(1:nx);

r.dsquodx = dsquodxe(1:nx);
r.dsquidx = dsquidxe(1:nx);
r.dsqlidx = dsqlidxe(1:nx);
r.dsqlodx = dsqlodxe(1:nx);

r.dCldx = dCldxe(1:nx);
r.dVtotdx = dVtotdxe(1:nx);
r.dAtotdx = dAtotdxe(1:nx);
r.dLtotdx = dLtotdxe(1:nx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Response to flux error "e"

r.dpsizrde = dpsizrdxe(:,nx+1);
r.dpsimagde = dpsimagdxe(nx+1);
r.dpsibryde = dpsibrydxe(nx+1);
r.dpsicde = dpsicdxe(:,nx+1);
r.dpsivde = dpsivdxe(:,nx+1);

r.dpresde = dpresdxe(:,nx+1);
r.dpprimede = dpprimedxe(:,nx+1);
r.dfpolde = dfpoldxe(:,nx+1);
r.dffprimde = dffprimdxe(:,nx+1);
r.dfprimede = dfprimedxe(:,nx+1);
r.djparFde = djparFdxe(:,nx+1);
r.dpcurrtde = dpcurrtdxe(:,nx+1);

r.dcpasmade = dcpasmadxe(nx+1);
r.dWthde = dWthdxe(nx+1);
r.dlide = dlidxe(nx+1);
r.dbetapde = dbetapdxe(nx+1);
r.dbetatde = dbetatdxe(nx+1);
r.dbetande = dbetandxe(nx+1);
r.dpsiplade = dpsipladxe(nx+1);
r.dbtsurfde = dbtsurfdxe(nx+1);
r.drcurde = drcurdxe(nx+1);
r.dzcurde = dzcurdxe(nx+1);
r.dy2de = dy2dxe(nx+1);
r.dy2nde = dy2ndxe(nx+1);
r.drmaxisde = drmaxisdxe(nx+1);
r.dzmaxisde = dzmaxisdxe(nx+1);

r.drsurfde = drsurfdxe(nx+1);
r.dzsurfde = dzsurfdxe(nx+1);
r.daminorde = daminordxe(nx+1);
r.dbminorde = dbminordxe(nx+1);
r.delongde = delongdxe(nx+1);

r.dtriude = dtriudxe(nx+1);
r.dtrilde = dtrildxe(nx+1);

r.dsquode = dsquodxe(nx+1);
r.dsquide = dsquidxe(nx+1);
r.dsqlide = dsqlidxe(nx+1);
r.dsqlode = dsqlodxe(nx+1);

r.dClde = dCldxe(nx+1);
r.dVtotde = dVtotdxe(nx+1);
r.dAtotde = dAtotdxe(nx+1);
r.dLtotde = dLtotdxe(nx+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Remember new plasma response
% if best fs is for very same x then replace it, otherwise replace worst fs
if max(abs(d.fl.x - dd(ii).fl.x)) > 0
  [~, ii] = min(fs); % replace worst since d and dd(ii) are for different x's
end
dd(ii) = d;
ee(ii).x = x;
ee(ii).psizr = psizr;
ee(ii).rmaxis = rmaxis;
ee(ii).zmaxis = zmaxis;
ee(ii).psimag = psimag;
ee(ii).psibry = psibry;
ee(ii).cpasma = cpasma;
ee(ii).li = li;
ee(ii).betap = betap;
ee(ii).alcfs = b.aminor;
ee(ii).zcur = zcur;
ee(ii).rb = rb;
ee(ii).zb = zb;
rr(ii).dpsizrdx = dpsizrdxe(:,1:nx);
rr(ii).dpsizrde = dpsizrdxe(:,nx+1);
rr(ii).drmaxisdx = drmaxisdxe(1:nx);
rr(ii).dzmaxisdx = dzmaxisdxe(1:nx);
rr(ii).dpsimagdx = dpsimagdxe(1:nx);
rr(ii).dpsibrydx = dpsibrydxe(1:nx);
rr(ii).dcpasmadx = dcpasmadxe(1:nx);
rr(ii).dlidx = dlidxe(1:nx);
rr(ii).dbetapdx = dbetapdxe(1:nx);
rr(ii).dzcurdx = dzcurdxe(1:nx);
rr(ii).drbdx = drbdxe(:,1:nx);
rr(ii).dzbdx = dzbdxe(:,1:nx);


if ncalc > 190 & ncalc < 220
%  save(['ncalc' num2str(ncalc)])
end
