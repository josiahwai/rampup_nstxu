function e = gsinit(c,init)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE: e = gsinit(c,init)
%
%  PURPOSE: Create initial data for gseq, gsdesign, gsevolve
%
%  INPUTS: config, the output from gsconfig.m
%          init, structure with initial equilibrium quantities:
%                psizr is required,
%                rg, zg required if init and config grids differ
%                ic (TokSys coil currents) or cc (EFIT coil currents)
%                iv (TokSys vessel currents) or vc (EFIT vessel currents)
%                psimag or rmaxis & zmaxis, psibry or rbbbs & zbbbs, 
%                sp or pres or pprime,
%                sf or fpol or ffprim & rzero & bzero
%
%  OUTPUTS: e, initial equilibrium
%
%  RESTRICTIONS: The initial equilibrium likely won't be perfectly 
%                converged unless made by gsdesign
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%
%  WRITTEN BY:  Anders Welander ON 2016-11-14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu0 = 4e-7*pi;
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % for interpolation

if isfield(init,'psizr')
  psizr = init.psizr;
else
  error('psizr must be a field in init')
end

if isfield(init,'rg') & isfield(init,'zg')
  if numel(c.rg)~=numel(init.rg) | numel(c.zg)~=numel(init.zg) | ...
    max(abs(c.rg(:)-init.rg(:)))/(c.rg(end)-c.rg(1)) > 1e-4 | ...
    max(abs(c.zg(:)-init.zg(:)))/(c.zg(end)-c.zg(1)) > 1e-4
    psizr = gs_interp2(init.rg,init.zg,init.psizr,c.rg,c.zg);
  end
else
  [nz, nr] = size(psizr);
  if nz ~= c.nz | nr ~= c.nr % this should never happen
    % Assume the grid in init covers the same area but of different density
    rg = linspace(c.rg(1),c.rg(end),nr);
    zg = linspace(c.zg(1),c.zg(end),nz)';
    psizr = gs_interp2(rg,zg,init.psizr,c.rg,c.zg);
  end
end


if isfield(init,'ic')
  ic = c.Pcci*init.ic(:);
else
  ic = zeros(c.nic,1);
end

if isfield(init,'iv')
  iv = c.Pvci*init.iv(:);
else
  iv = zeros(c.niv,1);
end

if isfield(init,'ih')
  ih = init.ih;
else
  ih = zeros(c.nih,1);
end


%%%%%%%%%%%%%%%%%% INITIALIZE PROFILES %%%%%%%%%%%%%%%%%%

% Create initial values for sf
if strcmp(upper(c.tokamak),'ITER')
  rzero =  6.2;
  bzero = -5.3;
elseif strcmp(upper(c.tokamak),'EAST')
  rzero = 1.965;
  bzero = 2.000;
elseif strcmp(upper(c.tokamak),'D3D')
  rzero =  1.6955;
  bzero = -1.6955;
else
  rzero = mean(c.rg);
  bzero = rzero;
end
if isfield(init,'sf') & numel(init.sf) == c.nsf
  sf = init.sf(:);
elseif isfield(init,'fpol') & max(abs(init.fpol)) > 0 & ...
  isfield(init,'psimag') & isfield(init,'psibry')
  g = mu0*(init.psibry-init.psimag)/2/pi;
  f = init.fpol(:).^2/2-init.fpol(end)^2/2;
  n = length(f);
  x = linspace(0,1,n)'; % x = psibar
  k = ones(n,1);        % to contain indices to spline regions
  for i = 1:c.nkf
    k(x > c.psikf(i)) = i;
  end
  x = x*ones(1,c.nsf); % x = psibar
  sf = pinv(c.mf0(k,:)+c.mf1(k,:).*x+c.mf2(k,:).*x.^2+c.mf3(k,:).*x.^3)*f/g;
  if c.constraints == 0
    sf(end) = init.fpol(end);
  end
elseif isfield(init,'ffprim')
  f = init.ffprim(:)/mu0;
  n = length(f);
  x = linspace(0,1,n)'; % x = psibar
  k = ones(n,1);        % to contain indices to spline regions
  for i = 1:c.nkf
    k(x > c.psikf(i)) = i;
  end
  x = x*ones(1,c.nsf); % x = psibar
  sf = pinv(c.mf1(k,:) + 2*c.mf2(k,:).*x + 3*c.mf3(k,:).*x.^2)*f;
  if c.constraints == 0
    if isfield(init,'rzero') & isfield(init,'bzero')
      sf(end) = init.rzero*init.bzero;
    else
      disp(['Warning gsinit: toroidal field unknown, ', ...
	'init should contain fpol or rzero, bzero. '])
	disp(['*************** Defaulting to rzero*bzero = ', num2str(rzero), ...
	'*' num2str(bzero), ' = ', num2str(rzero*bzero), ' ***************'])
      sf(end) = rzero*bzero;
    end
  end
end

% Create initial values for sp
if isfield(init,'sp') & numel(init.sp) == c.nsp
  sp = init.sp(:);
elseif isfield(init,'pres') & max(abs(init.pres)) > 0 & ...
  isfield(init,'psimag') & isfield(init,'psibry')
  g = (init.psibry-init.psimag)/2/pi;
  f = init.pres(:)-init.pres(end);
  n = length(f);
  x = linspace(0,1,n)'; % x = psibar
  k = ones(n,1);        % to contain indices to spline regions
  for i = 1:c.nkp
    k(x > c.psikp(i)) = i;
  end
  x = x*ones(1,c.nsp);
  sp = pinv(c.mp0(k,:)+c.mp1(k,:).*x+c.mp2(k,:).*x.^2+c.mp3(k,:).*x.^3)*f/g;
  if c.constraints == 0
    sp(end) = init.pres(end);
  end
elseif isfield(init,'pprime')
  f = init.pprime(:);
  n = length(f);
  x = linspace(0,1,n)'; % x = psibar
  k = ones(n,1);        % to contain indices to spline regions
  for i = 1:c.nkp
    k(x > c.psikp(i)) = i;
  end
  x = x*ones(1,c.nsp);
  sp = pinv(c.mp1(k,:) + 2*c.mp2(k,:).*x + 3*c.mp3(k,:).*x.^2)*f;
end

x = [ic; iv; ih; sf; sp];

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
nsf = c.nsf;
nsp = c.nsp;
nx = c.nx;
ny = c.ny;
nu = c.nu;
ix = c.ix;
pb = c.pb;
vp = c.vp;
mp0 = c.mp0;
mp1 = c.mp1;
mp2 = c.mp2;
mp3 = c.mp3;
fb = c.fb;
vf = c.vf;
mf0 = c.mf0;
mf1 = c.mf1;
mf2 = c.mf2;
mf3 = c.mf3;
mcc = c.mcc;
mcv = c.mcv;
mvv = c.mvv;
mpc = c.mpc;
mpv = c.mpv;
mgg = c.mgg;
psikp = c.psikp;
psikf = c.psikf;
nbdtest = c.nbdtest;
nfla = c.nfla;
nflc = c.nflc;
nfltest = c.nfltest;
dpsibarlimit = c.dpsibarlimit;
evolve_option = c.evolve_option;
Rhat = c.Rhat;
Vhat = c.Vhat;

  
% Analyze topology of psizr
b = gsboundary(c,psizr);

% Shorter names for topology data
plasma = b.plasma;
igg = b.igg;
egg = b.egg;
n = b.n;
rmaxis = b.rmaxis;
zmaxis = b.zmaxis;
psimag = b.psimag;
iia = b.resp.iia;
wa = b.resp.wa;
psibry = b.psibry;
iib = b.resp.iib;
wb = b.resp.wb;
yc = b.psizrc;
I36 = b.resp.I36;
ycw = b.resp.psizrc.w;
ycb = b.resp.psizrc.b;
ycr = b.psizrcr;
ycrw = b.resp.psizrcr.w;
ycrb = b.resp.psizrcr.b;
ycz = b.psizrcz;
yczw = b.resp.psizrcz.w;
yczb = b.resp.psizrcz.b;
RAc = b.RAcell;
RAcw = b.resp.RAcell.w;
RAcb = b.resp.RAcell.b;
ARc = b.ARcell;
ARcw = b.resp.ARcell.w;
ARcb = b.resp.ARcell.b;
rgc = b.rgc;
rgcw = b.resp.rgc.w;
rgcb = b.resp.rgc.b;
zgc = b.zgc;
zgcw = b.resp.zgc.w;
zgcb = b.resp.zgc.b;
dC = b.dC;
dCw = b.resp.dC.w;
dCb = b.resp.dC.b;
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
Cl = b.Cl;
Vtot = b.Vtot;
Atot = b.Atot;
Ltot = b.Ltot;
i49 = b.resp.di7x7';

% Flux from external currents
psizr_app = reshape(mpc*ic + mpv*iv, nz, nr);

if plasma

  % flag grid cells containing plasma
  iplasma = igg | egg;

  dba = psibry-psimag;

  % Normalized flux at grid points
  psibarzr = (psizr-psimag)/dba;

  % Normalized flux at current centers within cells, ync
  ync = (yc-psimag)/dba;

  % Index spline knots for ync
  ikp = ones(nz,nr);
  for i = 2:nkp
    ikp(ync > psikp(i)) = i; % knots on grid for pres
  end
  ikf = ones(nz,nr);
  for i = 2:nkf
    ikf(ync > psikf(i)) = i; % knots on grid for fpol
  end

  % Coefficients for all spline knots
  pb = pb + vp*sp;
  p0 = mp0*sp;
  p1 = mp1*sp;
  p2 = mp2*sp;
  p3 = mp3*sp;
  fb = fb + vf*sf;
  f0 = mf0*sf;
  f1 = mf1*sf;
  f2 = mf2*sf;
  f3 = mf3*sf;

  % pres, fpol^2/2 depend on g
  g = dba/(2*pi);

  % Pressure at cell current centers, rgc, zgc
  preszr = g*(p0(ikp)+p1(ikp).*ync+p2(ikp).*ync.^2+p3(ikp).*ync.^3) + pb;
  preszr(~iplasma) = 0;

  % Thermal energy within cells
  Wzr = 3*pi*RAc.*preszr;

  % Volume integrals of Bpol^2 for cells
  B2Vc = ARc.*(ycr.^2 + ycz.^2)/(2*pi);

  % fpol^2/2 at rgc, zgc
  hfp2zr = g*(f0(ikf)+f1(ikf).*ync+f2(ikf).*ync.^2+f3(ikf).*ync.^3) + fb^2/2;

  % pcurrt = Integral[(R*pprime+ffprim/mu0/R)*dAcell]
  % Ignore R,Z dependency in pprime, ffprim, use values at rgc, zgc
  ppzr = p1(ikp) + 2*p2(ikp).*ync + 3*p3(ikp).*ync.^2; 
  fmzr = f1(ikf) + 2*f2(ikf).*ync + 3*f3(ikf).*ync.^2;
  pcurrt = RAc.*ppzr + ARc.*fmzr;

  % psizr_pla = Integral(m(R,Z)*i(R,Z))*dA
  % Ignore R, Z dependencies within cells
  psizr_pla = reshape(mgg*pcurrt(:),nz,nr);

  % Convergence error, psizr_err
  psizr_err = psizr-psizr_app-psizr_pla;

  % Conductor circuit fluxes
  psic = mcc *ic + mcv*iv + mpc'*pcurrt(:);
  psiv = mcv'*ic + mvv*iv + mpv'*pcurrt(:);

  % Scalar parameters
  cpasma = sum(pcurrt(:));
  btsurf = fb/rsurf;
  Wth = sum(Wzr(:));
  B2V = sum(B2Vc(:));
  Cl = sum(dC(1:n-1));
  bp2flx = (mu0*cpasma/Cl)^2;
  psipla = yc(:)'*pcurrt(:)/cpasma;
  betap = 4/3*mu0*Wth/Vtot/bp2flx;
  ipnorm = abs(cpasma/1e6/aminor/btsurf);
  betat = 400/3*Wth/Vtot/btsurf^2*mu0;
  betan = betat/ipnorm;
  li = B2V/Vtot/bp2flx; % Volume-averaged Bp^2 / bp2flx
  rcur = sum(rgc(:).*pcurrt(:))/cpasma;
  zcur = sum(zgc(:).*pcurrt(:))/cpasma;

else % No plasma
  psizr_pla = zeros(nz,nr);
  psizr_err = zeros(nz,nr);
  pcurrt = zeros(nz,nr);
  preszr = zeros(nz,nr);
  Wzr = zeros(nz,nr);
  B2Vc = zeros(nz,nr);
  psic = mcc *ic + mcv*iv;
  psiv = mcv'*ic + mvv*iv;
  cpasma = 0;
  rsurf = nan;
  aminor = 0;
  btsurf = nan;
  Wth = 0;
  Cl = 0;
  bp2flx = nan;
  psipla = nan;
  Vtot = 0;
  betap = 0;
  ipnorm = 0;
  betat = 0;
  betan = 0;
  B2V = 0;
  li = nan;
  rcur = nan;
  zcur = nan;

end
% equilibrium structure, e
nbbbs = 2*nr-1;
e = struct('x', x, 'rg', rg, 'zg', zg, ...
  'psizr', psizr, 'psizr_app', psizr_app, 'psizr_pla', zeros(nz,nr), ...
  'psizr_err',zeros(nz,nr), 'psimag', nan, 'psibry', nan, ...
  'rbbbs', nan(nbbbs,1), 'zbbbs', nan(nbbbs,1), 'nbbbs', nbbbs, ...
  'rmaxis', rmaxis, 'zmaxis', zmaxis, 'rbdef', nan, 'zbdef', nan, ...
  'pres', zeros(1,nr)+c.pb+c.vp*sp, 'pprime', zeros(1,nr), ...
  'fpol', zeros(1,nr)+c.fb+c.vf*sf, 'ffprim', zeros(1,nr), ...
  'pcurrt', zeros(nz,nr), 'preszr', zeros(nz,nr), 'Wzr', zeros(nz,nr), ...
  'B2Vc', zeros(nz,nr), 'psic', psic, 'psiv', psiv, 'cpasma', 0, ...
  'rsurf', nan, 'aminor', 0, 'btsurf', nan, 'Wth', 0, 'Cl', 0, ...
  'bp2flx', 0, 'psipla', nan, 'Vtot', 0, 'betap', 0, 'ipnorm', 0, ...
  'betat', 0, 'betan', 0, 'B2V', 0, 'li', nan, 'rcur', nan, 'zcur', nan);

% Archive updated equilibrium
e.x = x;
e.psizr = psizr;
e.psizr_app = psizr_app;
e.psizr_pla = psizr_pla;
e.psizr_err = psizr_err;
e.psimag = psimag;
e.psibry = psibry;
e.pcurrt = pcurrt;
e.preszr = preszr;
e.Wzr = Wzr;
e.B2Vc = B2Vc;
e.psic = psic;
e.psiv = psiv;
e.cpasma = cpasma;
e.rsurf = rsurf;
e.aminor = aminor;
e.btsurf = btsurf;
e.Wth = Wth;
e.Cl = Cl;
e.bp2flx = bp2flx;
e.psipla = psipla;
e.Vtot = Vtot;
e.betap = betap;
e.ipnorm = ipnorm;
e.betat = betat;
e.betan = betan;
e.B2V = B2V;
e.li = li;
e.rcur = rcur;
e.zcur = zcur;

