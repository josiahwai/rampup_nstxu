function c = gsconfig(tok,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   c = gsconfig(tok,opt)
%
%  PURPOSE: Configure parameters used by gseq, gsdesign, gsupdate, GSevolve
%
%  INPUTS: tok, structure with TokSys tokamak info, tok_data_struct
%          opt, options
%          For complete information type   edit gsconfig_template.m
%
%  OUTPUTS: c, configuration data for gs codes
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  WRITTEN BY:  Anders Welander  ON 2016-11-06
%
%  MODIFICATION HISTORY: 2017-12-11 Changed the one input config to 
%    two inputs tok, opt but still allowing all fields to be in one input
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
  error('Need an input with TokSys description of tokamak, tok_data_struct')
end

config = tok;
if nargin > 1
  f = fields(opt);
  for i = 1:numel(f)
    nam = char(f(i));
    config = setfield(config,nam,getfield(opt,nam));
  end
end
mu0 = 4e-7*pi;

% Verbose if set will report about every recognized field in config
% and warn about unrecognized fields and omitted fields
if isfield(config,'verbose')
  verbose = config.verbose;
else
  verbose = 0;
end

if isfield(config,'tokamak')
  tokamak = config.tokamak;
else
  tokamak = 'Unknown';
  if verbose > 0
    disp('The name of this tokamak is not known, the field tokamak is missing')
  end
end

%%%%%%%%%%%%%%%%% CONFIGURE GRID %%%%%%%%%%%%%%%%%
if isfield(config,'rg')
  rg = config.rg(:)';
else
  error('Grid info required. Field rg (and maybe others) missing in config.')
end
nr = length(rg);
dr = (rg(nr)-rg(1))/(nr-1);
if isfield(config,'zg')
  zg = config.zg(:);
else
  error('Grid info required. Field zg (and maybe others) missing in config.')
end
nz = length(zg);
dz = (zg(nz)-zg(1))/(nz-1);
ngg = nr*nz;
rgg = ones(nz,1)*rg;
zgg = zg*ones(1,nr);
Ag = dr*dz;
RA = rgg*Ag;
AR = (log(rgg+dr/2)-log(rgg-dr/2))*dz;
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % Used for interpolation
ir4 = (-1:2)'*nz; % For finding flux on horizontal grid lines
iz4 = (-1:2)';    % For finding flux on vertical grid lines
ir8 = [iz4-nz;iz4+nz]; % For finding derivatives w.r.t. r on vertical grid lines
iz8 = [ir4- 1;ir4+ 1]; % For finding derivatives w.r.t. z on horizontal grid lines
i16 = reshape(iz4*ones(1,4)+ones(4,1)*ir4',1,16);
%%%%%%%%%%%%%%%%% GRID CONFIGURED %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% CONFIGURE LIMITER %%%%%%%%%%%%%%%%%
if isfield(config,'limdata') % Tolerate format variations that cause no ambiguity
  if size(config.limdata,1) < size(config.limdata,2)
    limdata = config.limdata';
  else
    limdata = config.limdata;
  end
  if any(limdata(:,2) < 0)
    Rlim = limdata(:,1);
    Zlim = limdata(:,2);
  else
    Zlim = limdata(:,1);
    Rlim = limdata(:,2);
  end
elseif isfield(config,'Rlim') & isfield(config,'Zlim')
  Rlim = config.Rlim;
  Zlim = config.Zlim;
else
  error('Either limdata or Rlim, Zlim needed in config with limiter geometry')
end
if Rlim(end) ~= Rlim(1) | Zlim(end) ~= Zlim(1)
  Rlim(end+1) = Rlim(1);
  Zlim(end+1) = Zlim(1);
end
% Remove unnecessary points
ilim = Rlim==Rlim;
j = 1;
for i = 2:length(Rlim)-1
  % Check if i is on a straight line between j and i+1
  dR = Rlim(i+1)-Rlim(j);
  dZ = Zlim(i+1)-Zlim(j);
  if abs(dR) > abs(dZ)
    r = Rlim(i);
    z = Zlim(j)+dZ/dR*(Rlim(i)-Rlim(j));
    onstraightline = abs(z-Zlim(i))/dz < 0.001;
  else
    z = Zlim(i);
    r = Rlim(j)+dR/dZ*(Zlim(i)-Zlim(j));
    onstraightline = abs(r-Rlim(i))/dr < 0.001;
  end
  if onstraightline | ...
    abs(Rlim(i)-Rlim(j))/dr < 0.001 & abs(Zlim(i)-Zlim(j))/dz < 0.001
    ilim(i) = 0;
  else
    j = i;
  end
end
nlim = length(Rlim);


% For quickly finding inside of limiter

limpoly.r = Rlim;
limpoly.z = Zlim;
limpoly.n = nlim;

% Unique z values
limpoly.unique.z = unique(limpoly.z);
limpoly.unique.n = length(limpoly.unique.z);

% Find ranges in r inside limpoly at unique z
% nranges will hold maximum number of ranges
nranges = 1;
% rranges are valid just above unique z
limpoly.rranges = nan(limpoly.unique.n,2*nranges);
limpoly.drdz = nan(limpoly.unique.n,2*nranges);
for i = 1:limpoly.unique.n
  z = limpoly.unique.z(i);
  for j = 1:limpoly.n-1
    dzp = limpoly.z(j+1)-limpoly.z(j);
    if limpoly.z(j) <= z & z < limpoly.z(j+1) || ...
       limpoly.z(j) > z & z >= limpoly.z(j+1)
      % interpolate to find r at z
      drp = limpoly.r(j+1)-limpoly.r(j);
      drdz = drp/dzp;
      f = (z-limpoly.z(j));
      r = limpoly.r(j) + f*drdz;
      % Put r & drdz in first free entry
      k = 1;
      while ~isnan(limpoly.rranges(i,k))
        k = k+1;
	if k > 2*nranges
	  nranges = nranges + 1;
	  limpoly.rranges(:,k:2*nranges) = nan;
	  limpoly.drdz(:,k:2*nranges) = nan;
	end
      end
      limpoly.rranges(i,k) = r;
      limpoly.drdz(i,k) = drdz;
    end % End of < = > tests
  end % End of j loop
  [limpoly.rranges(i,:), kranges] = sort(limpoly.rranges(i,:));
  limpoly.drdz(i,:) = limpoly.drdz(i,kranges);
  for j = 1:2:size(limpoly.rranges,2)-1 % Sorting won't work for zero ranges
    if limpoly.rranges(i,j) == limpoly.rranges(i,j+1)
      limpoly.drdz(i,j:j+1) = sort(limpoly.drdz(i,j:j+1)); % Fix for zero range
    end
  end
end % End of i loop
limpoly.nranges = nranges;
% End of processing limpoly


% Create limiter coordinates rl,zl, used for finding touch points
R = 2*(Rlim-rg(1))/dr+1;
Z = 2*(Zlim-zg(1))/dz+1;
if abs(R(1)-round(R(1))) < 1e-4
  R([1 end]) = round(R(1));
end
if abs(Z(1)-round(Z(1))) < 1e-4
  Z([1 end]) = round(Z(1));
end
rl = R(1);
zl = Z(1);
for i = 1:nlim-1
  if R(i+1) > R(i)
    r1 = floor(R(i)+1):ceil(R(i+1)-1);
  else
    r1 = ceil(R(i)-1):-1:floor(R(i+1)+1);
  end
  x1 = (r1-R(i))/(R(i+1)-R(i));
  z1 = Z(i)+x1*(Z(i+1)-Z(i));
  if Z(i+1) > Z(i)
    z2 = floor(Z(i)+1):ceil(Z(i+1)-1);
  else
    z2 = ceil(Z(i)-1):-1:floor(Z(i+1)+1);
  end
  x2 = (z2-Z(i))/(Z(i+1)-Z(i));
  r2 = R(i)+x2*(R(i+1)-R(i));
  x = [x1 x2];
  r = [r1 r2];
  z = [z1 z2];
  [x,k] = sort(x);
  rl(end+[1:length(x)]) = r(k);
  zl(end+[1:length(x)]) = z(k);
  rl(end+1) = R(i+1);
  zl(end+1) = Z(i+1);
end
nl = length(rl);
% Scale to floating index
rl = 1+(rl-1)/2;
zl = 1+(zl-1)/2;

if any(rl < 1 | rl > nr | zl < 1 | zl > nz)
  disp('Warning gsconfigure: Limiter extends outside grid')
end
drl = diff(rl)';
dzl = diff(zl)';
dl = sqrt(drl.^2+dzl.^2);
irl = min(nr-2,max(2,floor(min([rl(1:nl-1); rl(2:nl)]))))';
izl = min(nz-2,max(2,floor(min([zl(1:nl-1); zl(2:nl)]))))';
trl = rl(1:nl-1)'-irl;
tzl = zl(1:nl-1)'-izl;
trl2 = rl(2:nl)'-irl;
tzl2 = zl(2:nl)'-izl;
kl = izl+(irl-1)*nz;
il = [kl-(nz+1)   kl-nz     kl-(nz-1)   kl-(nz-2) ...
      kl-1        kl        kl+1        kl+2      ...
      kl+(nz-1)   kl+nz     kl+(nz+1)   kl+(nz+2) ...
      kl+(2*nz-1) kl+(2*nz) kl+(2*nz+1) kl+(2*nz+2)];

r = [trl.^0 trl.^1 trl.^2 trl.^3]*mx;
z = [tzl.^0 tzl.^1 tzl.^2 tzl.^3]*mx;
wl = ...
  [z(:,1).*r(:,1) z(:,2).*r(:,1) z(:,3).*r(:,1) z(:,4).*r(:,1) ...
   z(:,1).*r(:,2) z(:,2).*r(:,2) z(:,3).*r(:,2) z(:,4).*r(:,2) ...
   z(:,1).*r(:,3) z(:,2).*r(:,3) z(:,3).*r(:,3) z(:,4).*r(:,3) ...
   z(:,1).*r(:,4) z(:,2).*r(:,4) z(:,3).*r(:,4) z(:,4).*r(:,4)];

r = [0*trl trl.^0 2*trl 3*trl.^2]*mx;
z = [tzl.^0 tzl.^1 tzl.^2 tzl.^3]*mx;
wlr = ...
  [z(:,1).*r(:,1) z(:,2).*r(:,1) z(:,3).*r(:,1) z(:,4).*r(:,1) ...
   z(:,1).*r(:,2) z(:,2).*r(:,2) z(:,3).*r(:,2) z(:,4).*r(:,2) ...
   z(:,1).*r(:,3) z(:,2).*r(:,3) z(:,3).*r(:,3) z(:,4).*r(:,3) ...
   z(:,1).*r(:,4) z(:,2).*r(:,4) z(:,3).*r(:,4) z(:,4).*r(:,4)];

r = [trl.^0 trl.^1 trl.^2 trl.^3]*mx;
z = [0*tzl tzl.^0 2*tzl 3*tzl.^2]*mx;
wlz = ...
  [z(:,1).*r(:,1) z(:,2).*r(:,1) z(:,3).*r(:,1) z(:,4).*r(:,1) ...
   z(:,1).*r(:,2) z(:,2).*r(:,2) z(:,3).*r(:,2) z(:,4).*r(:,2) ...
   z(:,1).*r(:,3) z(:,2).*r(:,3) z(:,3).*r(:,3) z(:,4).*r(:,3) ...
   z(:,1).*r(:,4) z(:,2).*r(:,4) z(:,3).*r(:,4) z(:,4).*r(:,4)];

% Calculate derivative toward higher limiter index
wld1 = drl./dl*ones(1,16).*wlr + dzl./dl*ones(1,16).*wlz;

r = [0*trl2 trl2.^0 2*trl2 3*trl2.^2]*mx;
z = [tzl2.^0 tzl2.^1 tzl2.^2 tzl2.^3]*mx;
wlr = ...
  [z(:,1).*r(:,1) z(:,2).*r(:,1) z(:,3).*r(:,1) z(:,4).*r(:,1) ...
   z(:,1).*r(:,2) z(:,2).*r(:,2) z(:,3).*r(:,2) z(:,4).*r(:,2) ...
   z(:,1).*r(:,3) z(:,2).*r(:,3) z(:,3).*r(:,3) z(:,4).*r(:,3) ...
   z(:,1).*r(:,4) z(:,2).*r(:,4) z(:,3).*r(:,4) z(:,4).*r(:,4)];

r = [trl2.^0 trl2.^1 trl2.^2 trl2.^3]*mx;
z = [0*tzl2 tzl2.^0 2*tzl2 3*tzl2.^2]*mx;
wlz = ...
  [z(:,1).*r(:,1) z(:,2).*r(:,1) z(:,3).*r(:,1) z(:,4).*r(:,1) ...
   z(:,1).*r(:,2) z(:,2).*r(:,2) z(:,3).*r(:,2) z(:,4).*r(:,2) ...
   z(:,1).*r(:,3) z(:,2).*r(:,3) z(:,3).*r(:,3) z(:,4).*r(:,3) ...
   z(:,1).*r(:,4) z(:,2).*r(:,4) z(:,3).*r(:,4) z(:,4).*r(:,4)];

% Calculate derivative at higher limiter index coming from lower
wld2 = drl./dl*ones(1,16).*wlr + dzl./dl*ones(1,16).*wlz;

% Check if the left or right side is inside when walking along limiter
inside = inpolygon((rl(1:nl-1)+rl(2:nl))/2+diff(zl)/1e9, ...
                  (zl(1:nl-1)+zl(2:nl))/2-diff(rl)/1e9,rl,zl);
if sum(inside) == 0
  turnin = [0 -1;+1 0];
else
  turnin = [0 +1;-1 0];
end

anglel = angle(diff(rl([1:nl 2])).*diff(rl([nl-1 1:nl]))+...
  diff(zl([1:nl 2])).*diff(zl([nl-1 1:nl])) + ...
  1i*(turnin(1,2)*diff(zl([1:nl 2])).*diff(rl([nl-1 1:nl]))+...
      turnin(2,1)*diff(rl([1:nl 2])).*diff(zl([nl-1 1:nl]))));
      
concavel = anglel > 0.001; % Corner protrudes into machine
convexl = anglel < -0.001;

if 0
  % Confirmation plot
  clf
  hold on
  for i = 1:limpoly.unique.n
    z2 = limpoly.unique.z(i)+[0 0];
    for k = 1:limpoly.nranges
      kk = 2*k-[1 0];
      r2 = limpoly.rranges(i,kk);
      plot(r2,z2,'LineWidth',3)
    end
  end
  plot(limpoly.r,limpoly.z,'r')
end
%%%%%%%%%%%%%%%%% LIMITER CONFIGURED %%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%% CONFIGURE SPLINES %%%%%%%%%%%%%%%%%

if strcmp(upper(tokamak),'ITER')
  rzero =  6.2;
  bzero = -5.3;
elseif strcmp(upper(tokamak),'D3D')
  rzero =  1.6955;
  bzero = -1.6955;
elseif strcmp(upper(tokamak),'EAST')
  rzero = 1.965;
  bzero = 2.000;
elseif strcmp(upper(tokamak),'KSTAR')
  rzero =  1.8;
  bzero =  2.39;
else
  rzero = mean(rg);
  bzero = rzero;
end

if isfield(config,'constraints')
  constraints = config.constraints;
elseif isfield(config,'constraints')
  constraints = config.constraints;
else
  constraints = 0;
end

% Normalized flux for nr equally spaced values from axis to boundary
psibar = linspace(0,1,nr)';

% pres and fpol^2/2 described by splined third degree polynomials
if isfield(config,'psikn') & length(config.psikn) > 1 & ~any(isnan(config.psikn))
  psikn = config.psikn; % psibar for knot junctions
elseif isfield(config,'psikn') & length(config.psikn) > 1 & ~any(isnan(config.psikn))
  psikn = config.psikn; % psibar for knot junctions
elseif isfield(config,'nkn')
  psikn = linspace(0,1,config.nkn+1); % default is equally spaced values
else
  psikn = [0 1];
end
if isfield(config,'psikp') & length(config.psikp) > 1 & ~any(isnan(config.psikp))
  psikp = config.psikp;
else
  psikp = psikn;
end
if isfield(config,'psikf') & length(config.psikf) > 1 & ~any(isnan(config.psikf))
  psikf = config.psikf;
else
  psikf = psikn;
end
if constraints > 0
  psikp = linspace(0,1,round((nr+1)/2));
  psikf = linspace(0,1,round((nr+1)/2));
end
nkp = length(psikp)-1; % number of third degree polynomials in spline
nsp = nkp+3; % Number of states for unconstrained pressure profile
nkf = length(psikf)-1; % number of third degree polynomials in spline
nsf = nkf+3; % Number of states for unconstrained fpol profile

% mp0 is matrix to calculate pprime coefficient 0
% vp is vector to calculate pressure at boundary from states sp
% pb is pressure at boundary if no state contains it, otherwise 0
% p0 = mp0*sp; p1 = mp1*sp; p2 = mp2*sp; p3 = mp3*sp; 
% presb = pb + vp*sp;
% pres = (psibry-psimag)/2/pi*(p0+p1*psibar+p2*psibar^2+p3*psibar^3) + presb
% pprime = p1+2*p2*psibar+3*p3*psibar^2
% For unconstrained profiles:
% sp(1:nkn) are values of p3 for each interval
% sp(nkn+1) is value of p2 for last interval
% sp(nkn+2) is value of p1 for last interval
% sp(nkn+3) is pressure at boundary
mp0 = zeros(nkp,nkp+3); mp0(nkp,nkp+[0:2]) = -1;
mp1 = zeros(nkp,nkp+3); mp1(nkp,nkp+2) = 1;
mp2 = zeros(nkp,nkp+3); mp2(nkp,nkp+1) = 1;
mp3 = zeros(nkp,nkp+3); mp3(nkp,nkp+0) = 1;
for i = nkp-1:-1:1
  x = psikp(i+1);
  mp3(i,i) = 1;
  % Match second derivative at knot i, i.e. 2*c2+6*c3*x
  d3 = mp3(i+1,:)-mp3(i,:);
  mp2(i,:) = mp2(i+1,:) + 3*x*d3;
  % Match first derivative at knot i, i.e. c1+2*c2*x+3*c3*x^2
  d2 = mp2(i+1,:)-mp2(i,:);
  mp1(i,:) = mp1(i+1,:) + 2*x*d2 + 3*x^2*d3;
  % Match value at knot i, i.e. c0+c1*x+c2*x^2+c3*x^3
  d1 = mp1(i+1,:)-mp1(i,:);
  mp0(i,:) = mp0(i+1,:) + x*d1 + x^2*d2 + x^3*d3;
end
pb = 0;
vp = [zeros(1,nkp+2) 1];
ikp = ones(nr,1); % Spline region at each of the psibar values
for i = 2:nkp
  ikp(psibar > psikp(i)) = i;
end

% mf0 is matrix to calculate ffprim/mu0 coefficient 0
% vf is vector to calculate fpol at boundary from states sf
% fb is fpol at boundary if no state contains it, otherwise 0
% f0 = mf0*sf; f1 = mf1*sf; f2 = mf2*sf; f3 = mf3*sf; 
% fpolb = fb + vf*sf;
% fpol^2/2 = (psibry-psimag)/2/pi*(f0+f1*psibar+f2*psibar^2+f3*psibar^3) + fpolb^2/2
% ffprim = f1+2*f2*psibar+3*f3*psibar^2
% For unconstrained profiles:
% sf(1:nkn) are values of f3 for each interval
% sf(nkn+1) is value of f2 for last interval
% sf(nkn+2) is value of f1 for last interval
% sf(nkn+3) is fpol at boundary
mf0 = zeros(nkf,nkf+3); mf0(nkf,nkf+[0:2]) = -1;
mf1 = zeros(nkf,nkf+3); mf1(nkf,nkf+2) = 1;
mf2 = zeros(nkf,nkf+3); mf2(nkf,nkf+1) = 1;
mf3 = zeros(nkf,nkf+3); mf3(nkf,nkf+0) = 1;
for i = nkf-1:-1:1
  x = psikf(i+1);
  mf3(i,i) = 1;
  % Match second derivative at knot i, i.e. 2*c2+6*c3*x
  d3 = mf3(i+1,:)-mf3(i,:);
  mf2(i,:) = mf2(i+1,:) + 3*x*d3;
  % Match first derivative at knot i, i.e. c1+2*c2*x+3*c3*x^2
  d2 = mf2(i+1,:)-mf2(i,:);
  mf1(i,:) = mf1(i+1,:) + 2*x*d2 + 3*x^2*d3;
  % Match value at knot i, i.e. c0+c1*x+c2*x^2+c3*x^3
  d1 = mf1(i+1,:)-mf1(i,:);
  mf0(i,:) = mf0(i+1,:) + x*d1 + x^2*d2 + x^3*d3;
end
fb = 0;
vf = [zeros(1,nkf+2) 1];
ikf = ones(nr,1); % Spline region at each of the psibar values
for i = 2:nkf
  ikf(psibar > psikf(i)) = i;
end

% If profiles are constrained to 3 degrees of freedom (tied to cpasma, li, betap)
% then states for pres(1), fpol(1) are removed, sp reduced to 1 and sf to 2 states

defaultp = false;
if isfield(config,'pres0')
  pb0 = config.pres0(end);
  f = config.pres0(:)-pb0;
  n = length(f);
  x = linspace(0,1,n)';
  k = ones(n,1);
  for i = 2:nkp
    k(x > psikp(i)) = i;
  end
  x = x*ones(1,nsp);
  sp0 = pinv(mp0(k,:)+mp1(k,:).*x+mp2(k,:).*x.^2+mp3(k,:).*x.^3)*f;
else
  pb0 = 0;
  if isfield(config,'pprime0')
    f = config.pprime0(:);
  else
    f = ones(nr,1);
    defaultp = true;
  end
  n = length(f);
  x = linspace(0,1,n)';
  k = ones(n,1);
  for i = 2:nkp
    k(x > psikp(i)) = i;
  end
  x = x*ones(1,nsp);
  sp0 = pinv(mp1(k,:)+2*mp2(k,:).*x+3*mp3(k,:).*x.^2)*f;
end
sf0 = zeros(nsf,1);
defaultf = false;
if isfield(config,'fpol0')
  fb0 = config.fpol0(end);
  f = config.fpol0(:).^2/2-fb0^2/2;
  n = length(f);
  x = linspace(0,1,n)';
  k = ones(n,1);
  for i = 2:nkf
    k(x > psikf(i)) = i;
  end
  x = x*ones(1,nsf);
  sf0 = pinv(mf0(k,:)+mf1(k,:).*x+mf2(k,:).*x.^2+mf3(k,:).*x.^3)*f;
else
  fb0 = rzero*bzero;
  if isfield(config,'ffprim0')
    f = config.ffprim0(:);
  else
    f = ones(nr,1);
    defaultf = true;
  end
  n = length(f);
  x = linspace(0,1,n)';
  k = ones(n,1);
  for i = 2:nkf
    k(x > psikf(i)) = i;
  end
  x = x*ones(1,nsf);
  sf0 = pinv(mf1(k,:)+2*mf2(k,:).*x+3*mf3(k,:).*x.^2)*f;
end
x = psibar*ones(1,nsf);
sg0 = pinv(mf1(ikf,:)+2*mf2(ikf,:).*x+3*mf3(ikf,:).*x.^2)*(1-psibar)/mu0;

if constraints > 0
  nsp = 1;
  % Formula for coefficients will still be: p0 = mp0*sp, etc
  mp0 = mp0*sp0;
  mp1 = mp1*sp0;
  mp2 = mp2*sp0;
  mp3 = mp3*sp0;
  % Formula for boundary will still be: pres(psibar=1) = pb + vp*s
  pb = pb0;
  vp = 0;
  nsf = 2;
  % Formula for coefficients will still be: f0 = mf0*sf, etc
  mf0 = mf0*[sf0 sg0];
  mf1 = mf1*[sf0 sg0];
  mf2 = mf2*[sf0 sg0];
  mf3 = mf3*[sf0 sg0];
  % Formula for boundary will still be: fpol(psibar=1) = fb + vf*s
  fb = fb0;
  vf = [0 0];
end

%%%%%%%%%%%%%%%%% SPLINES CONFIGURED %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% CONFIGURE MUTUALS %%%%%%%%%%%%%%%%%
if isfield(config,'mcc')
  nc = size(config.mcc,1);
else
  error('Field mcc (and maybe others) are missing from config')
end
if isfield(config,'mvv')
  nv = size(config.mvv,1);
else
  error('Field mvv (and maybe others) are missing from config')
end

% TokSys coil currents = Pcc * circuit currents
if isfield(config,'Pcc')
  if size(config.Pcc,1) ~= nc
    error(['size(config.Pcc,1) must equal ' num2str(nc) 10 ...
      'TokSys coil currents = config.Pcc * currents in circuits'])
  end
  if isfield(config,'cccirc')
    disp(['Warning gsconfig: Circuits were specified with Pcc ' ...
      'and cccirc, Pcc will be used, cccirc will be ignored'])
  end
  if isfield(config,'buscode')
    disp(['Warning gsconfig: Circuits were specified with Pcc ' ...
      'and buscode, Pcc will be used, buscode will be ignored'])
  end
  Pcc = config.Pcc;
  nic = size(Pcc,2);
elseif isfield(config,'cccirc')
  nic = max(abs(config.cccirc));
  Pcc = zeros(nc,nic);
  if length(config.cccirc(:)) ~= nc
    error(['Length of cccirc differs from size(config.mcc,1). ' ...
      'Tokamak is ' upper(tokamak) ', with ' num2str(nc) ...
      ' coils but length(config.cccirc) = ' num2str(length(config.cccirc))])
  end
  for i = 1:nc
    j = abs(config.cccirc(i));
    Pcc(i,j) = sign(config.cccirc(i));
  end
  if isfield(config,'buscode')
    disp(['Warning gsconfig: Circuits were specified with both ' ...
      'cccirc and buscode, cccirc will be used, buscode will be ignored'])
  end
elseif isfield(config,'buscode')
  Pcc = eye(nc);
  ind = find(config.buscode);
  Pcc(min(ind),ind) = -1;
  Pcc = Pcc(:,[1:ind(1)-1 ind(1)+1:nc]);
  nic = nc-1;
else
  Pcc = eye(nc);
  nic = nc;
end
Pcci = pinv(Pcc);

% TokSys vessel currents = Pvc * circuit vessel currents
if isfield(config,'Pvc')
  if size(config.Pvc,1) ~= nv
    error(['size(config.Pvc,1) must equal ' num2str(nv) 10 ...
      'TokSys vessel currents = config.Pvc * currents in vessel circuits'])
  end
  if isfield(config,'vccirc')
    disp(['Warning gsconfig: Vessel circuits were specified with Pvc ' ...
      'and vccirc, Pvc will be used, vccirc will be ignored'])
  end
  Pvc = config.Pvc;
  niv = size(Pvc,2);
elseif isfield(config,'vccirc')
  niv = max(abs(config.vccirc));
  Pvc = zeros(nv,niv);
  if length(config.vccirc(:)) ~= nv
    error(['Length of vccirc differs from size(config.mvv,1). ' ...
      'Tokamak is ' upper(tokamak) ', with ' num2str(nv) ...
      ' vessel elements but length(config.vccirc) = ' num2str(length(config.vccirc))])
  end
  for i = 1:nv
    j = abs(config.vccirc(i));
    Pvc(i,j) = sign(config.vccirc(i));
  end
else
  Pvc = eye(nv);
  niv = nv;
end
Pvci = pinv(Pvc);

% Read in conductor system
mcc = Pcc'*config.mcc*Pcc;
if isfield(config,'mcv')
  mcv = Pcc'*config.mcv*Pvc;
else
  error('Field mcv (and maybe others) are missing from config')
end
mvv = Pvc'*config.mvv*Pvc;
if isfield(config,'mpc')
  mpc = config.mpc*Pcc;
else
  'room for code to calculate mpc'
end
if isfield(config,'mpv')
  mpv = config.mpv*Pvc;
else
  'room for code to calculate mpv'
end
if isfield(config,'mpp')
  mpp = config.mpp;
else
  'room for code to calculate mpp'
end

% Resistances for all conductors
if isfield(config,'ress')
  ress = config.ress(:);
  resc = ress(1:nc);
  resv = ress(nc+(1:nv));
elseif isfield(config,'resc') & isfield(config,'resv')
  resc = config.resc(:);
  resv = config.resv(:);
else
  error('Conductor resistances resc and resv (or ress) needed in config')
end
P = [Pcc zeros(nc,niv); zeros(nv,nic) Pvc];
rss = P'*diag([resc; resv])*P;

% Inductances for all conductors
mss = [mcc mcv; mcv' mvv];

% Extra resistance for all conductors
if isfield(config,'Rext')
  if min(size(config.Rext)) == 1
    Rext = diag(config.Rext);
  else
    Rext = config.Rext;
  end
  rssx = P'*Rext*P;
else
  rssx = zeros(nic+niv);
end

% Extra inductance for all conductors
if isfield(config,'Lext')
  if min(size(config.Lext)) == 1
    Lext = diag(config.Lext);
  else
    Lext = config.Lext;
  end
  mssx = P'*Lext*P;
else
  mssx = zeros(nic+niv);
end

% power supply voltage = (rss+rssx)*is + (mss+mssx)*isdot + plasma influence

mgg = zeros(ngg);
for j = 1:nz
  izshift = (1+abs(j-(1:nz))-(1:nz))'*ones(1,nr);
  mgg(:,j+(0:nr-1)*nz) = mpp((1:ngg)'+izshift(:),:);
end

% Flux loops
if isfield(config,'mlc')
  mlc = config.mlc*Pcc;
else
  mlc = zeros(0,nic);
end
if isfield(config,'mlv')
  mlv = config.mlv*Pvc;
else
  mlv = zeros(0,niv);
end
if isfield(config,'mpl')
  mpl = config.mpl;
else
  mpl = zeros(ngg,0);
end
nfl = size(mlc,1);

% Voltage loops
if isfield(config,'mhc')
  mhc = config.mhc*Pcc;
else
  mhc = zeros(0,nic);
end
if isfield(config,'mhv')
  mhv = config.mhv*Pvc;
else
  mhv = zeros(0,niv);
end
if isfield(config,'mph')
  mph = config.mph;
else
  mph = zeros(ngg,0);
end
nlv = size(mhc,1);

% Magnetic probes
if isfield(config,'gbc')
  gbc = config.gbc*Pcc;
else
  gbc = zeros(0,nic);
end
if isfield(config,'gbv')
  gbv = config.gbv*Pvc;
else
  gbv = zeros(0,niv);
end
if isfield(config,'gpb')
  gpb = config.gpb;
else
  gpb = zeros(ngg,0);
end
nbp = size(gbc,1);

% Rogowski loops
if isfield(config,'rldata')
  if size(config.rldata,1) == nc+nv+1 % This is what it should be
    rldata = config.rldata'; % The transpose is used in gs codes
  elseif size(config.rldata,2) == nc+nv+1 % Transposed? We'll take it!
    rldata = config.rldata;
  else
    rldata = zeros(0,nc+nv+1);
  end
else
  rldata = zeros(0,nc+nv+1);
end
rldata = [rldata(:,1:nc)*Pcc rldata(:,nc+(1:nv))*Pvc rldata(:,end)];
nrog = size(rldata,1);


%%%%%%%%%%%%%%%%% MUTUALS CONFIGURED %%%%%%%%%%%%%%%%%

if isfield(config,'limits')
  limits = config.limits;
else
  limits = [];
end
if ~isfield(limits,'ic')
  if strcmp(upper(tokamak),'D3D')
     % Ampere limits for E&F coils
    limits.ic = [2e3 2e3 ...
      59 59 59 59 54 124 124 59 89 ...
      59 59 59 59 54 124 124 59 89]'*[-100 100];
  end
  if strcmp(upper(tokamak),'EAST')
    % Limits for PF1 to PF14 and IC1, IC2, taken from PCS settings:
    limits.ic = [12.8*ones(6,1);11.6*ones(4,1);10.2*ones(4,1);5;5]*...
      [-1000 1000];
  end
  if strcmp(upper(tokamak),'KSTAR')
    % Limits taken from PCS settings for shot 5516 for coils in the order:
    % PF1-7U, PF1-7L, IVCU, IRCU, IVCL, IRCL
    limits.ic = [7.5 8.5 9 9 7 3.5 3.5 7.5 8.5 9 9 7 3.5 3.5 3 8 3 8]'*...
      [-1000 1000];
  end
  if strcmp(upper(tokamak),'ITER')
    % Limits taken from ITER_data_2010-v3.3.xls
    % PF1 PF2 PF3 PF4 PF5 PF6 CS3L CS2L CS1L CS1U CS2U CS3U VS3U VS3L
    % The high B limits, 0.4 K sub cooling PF6
    limits.bc = [65 50 50 50 60 70 130 130 130 130 130 130 nan nan]'*[-0.10 0.10];
    limits.ic = [41 50 50 50 33 41  40  40  40  40  40  40  10  10]'*[-1000 1000];
    % The low B limits, 0.4 K sub cooling PF6
    limits.bc = [64 48 48 48 57 68 126 126 126 126 126 126 nan nan]'*[-0.10 0.10];
    limits.ic = [48 55 55 55 52 52  45  45  45  45  45  45 10  10]'*[-1000 1000];
    % The high B limits
    limits.bc = [65 50 50 50 60 65 130 130 130 130 130 130 nan nan]'*[-0.10 0.10];
    limits.ic = [41 50 50 50 33 41  40  40  40  40  40  40 10  10]'*[-1000 1000];
    % The low B limits
    limits.bc = [64 48 48 48 57 64 126 126 126 126 126 126 nan nan]'*[-0.10 0.10];
    limits.ic = [48 55 55 55 52 48  45  45  45  45  45  45 10  10]'*[-1000 1000];
  end
end

% Halo model, 0=no model, 1=scalar model, 2=to be developed
if isfield(config,'halo')
  halo = config.halo;
else
  halo = 0; % 0 = no halo
end
% jhalo(nz,nr) is halo current density when halo = 1
if isfield(config,'jhalo')
  jhalo = config.jhalo;
else
  jhalo = 1; % uniform current density in halo
end
if halo == 0
  nih = 0;
else
  nih = 1;
end
[Acell, RAcell, ARcell] = polycellint(rg,zg,Rlim,Zlim);

% For contouring
if isfield(config,'nb') % Number of different points on a contour
  nb = config.nb;
else
  nb = 2*nr-2;
end
if isfield(config,'psibarc')
  psibarc = config.psibarc;
else
  psibarc = psibar;
end
np = length(psibarc);

% Indices in the state vector, x
ix.ic = [1:nic];
ix.iv = [1:niv]+nic;
ix.ih = [1:nih]+nic+niv;
ix.sf = [1:nsf]+nic+niv+nih;
ix.sp = [1:nsp]+nic+niv+nih+nsf;
nx = nic+niv+nih+nsf+nsp;

if isfield(config,'evolve_option')
  evolve_option  = config.evolve_option;
else
  evolve_option  = 0;
end

% Actuators are connected to the ic, ih, sf, sp states
if evolve_option == 2
  iu.ic = [1:nic];
  iu.ih = [1:nih]+nic;
  iu.sf = [1:np]+nic+nih;
  iu.sp = [1:nsp]+nic+nih+np;
  nu = nic+nih+np+nsp;
else
  iu.ic = [1:nic];
  iu.ih = [1:nih]+nic;
  iu.sf = [1:nsf]+nic+nih;
  iu.sp = [1:nsp]+nic+nih+nsf;
  nu = nic+nih+nsf+nsp;
end

% Maximum number of traced boundary points
nbmax = 4*(nr+nz);

% Smallest minor radius for which all current is inside last-closed-flux-surface
if isfield(config,'amin')
  amin = config.amin;
else
  amin = (max(Rlim)-min(Rlim))/4;
end

% max number of old d,e,r,t to remember in gsupdate
if isfield(config,'nn')
 nn = config.nn;
else
 nn = 9;
end

% max number of points to monitor for change of boundary-defining point
nbdtest = 25;

% number of points to monitor for change of normalized flux in plasma
nfla = 16;
nflc = 4;
nfltest = nfla*nflc;

if isfield(config,'dpsibarlimit')
 dpsibarlimit = config.dpsibarlimit;
else
 dpsibarlimit = 5e-3;
end

% Configure diagnostic points rdp, zdp
if isfield(config,'rdp') & isfield(config,'zdp')
  lrdp = length(config.rdp(:));
  lzdp = length(config.zdp(:));
  ndp = max(lrdp, lzdp)*(min(lrdp, lzdp) > 0);
  if lrdp == ndp
    rdp = config.rdp(:);
  else
    rdp = config.rdp(1)+zeros(ndp,1);
  end
  if lzdp == ndp
    zdp = config.zdp(:);
  else
    zdp = config.zdp(1)+zeros(ndp,1);
  end
  if lrdp > 1 & lzdp > 1 & lrdp ~= lzdp
    disp('Warning gsconfig: The sizes of diagnostic points rdp and zdp differ!')
    disp('They need to be two equally sized arrays or one scalar and one array.')
  end
else
  rdp = [];
  zdp = [];
end
ndp = length(rdp);
mdc = zeros(ndp,nc);
mdv = zeros(ndp,nv);
mpd = zeros(ngg,ndp);
grdc = zeros(ndp,nc);
grdv = zeros(ndp,nv);
grpd = zeros(ngg,ndp);
gzdc = zeros(ndp,nc);
gzdv = zeros(ndp,nv);
gzpd = zeros(ngg,ndp);
idp =  rdp > rg(1) & rdp < rg(nr) & zdp > zg(1) & zdp < zg(nz);
for j = 1:nic
  [mdc(idp,j) grdc(idp,j) gzdc(idp,j)] = gs_interp2(rg,zg,mpc(:,j),rdp(idp),zdp(idp));
end
for j = 1:niv
  [mdv(idp,j) grdv(idp,j) gzdv(idp,j)] = gs_interp2(rg,zg,mpv(:,j),rdp(idp),zdp(idp));
end
for j = 1:ngg
  [mpd(j,idp) grpd(j,idp) gzpd(j,idp)] = gs_interp2(rg,zg,mgg(:,j),rdp(idp),zdp(idp));
end

% Configure diagnostic vector of x-points
if isfield(config,'nxpoints')
  nxpoints = config.nxpoints;
end
if ~exist('nxpoints','var') | isempty(nxpoints)
  nxpoints = 2;
end
if nxpoints < 1
  nxpoints = 1;
end
if isfield(config,'xpointorder')
  xpointorder = config.xpointorder;
end
if ~exist('xpointorder','var') | isempty(xpointorder)
  xpointorder = 1;
end
if xpointorder < 1 || xpointorder > 4
  disp('Warning gsconfig: xpointorder must be between 1 and 4, setting it to 1.')
  xpointorder = 1;
end

% Configure outputs
need_mrmz = 0;
if isfield(config,'outputs')
  outputs = config.outputs;
end
if ~exist('outputs','var') | isempty(outputs)
  outputs = '';
end
outputs(outputs == 0) = 32;
iy = [];
ny  = 0; % Will hold number of outputs = length(y)
ds = []; % Will hold descriptions for the fields of iy
iE = false; % flag indices to nonzero elements of E
for i = 1:size(outputs,1)
  str = deblank(outputs(i,:));
  if isfield(iy,str)
    disp(['Warning in gsconfig: output ', str, ...
     ' appears more than once in outputs, only first entry used.'])
  elseif strcmp(str,'x')
    ds.x = 'state vector';
    iy.x = ny+(1:nx);
 iE(iy.x) = false;
    ny = ny+nx;
  elseif strcmp(str,'psizr')
    ds.psizr = 'Poloidal flux on the grid';
    iy.psizr = ny+(1:ngg);
 iE(iy.psizr) = false;
    ny = ny+ngg;
  elseif strcmp(str,'pcurrt')
    ds.pcurrt = 'Poloidal current within grid cells';
    iy.pcurrt = ny+(1:ngg);
 iE(iy.pcurrt) = false;
    ny = ny+ngg;
  elseif strcmp(str,'it')
    ds.it = 'Current in toroidal field coil';
    iy.it = ny+1;
 iE(iy.it) = false;
    ny = ny+1;
  elseif strcmp(str,'ic')
    ds.ic = 'Coil circuit currents producing flux = c.mpc*ic';
    iy.ic = ny+(1:nic);
 iE(iy.ic) = false;
    ny = ny+nic;
  elseif strcmp(str,'iv')
    ds.iv = 'Vessel circuit currents producing flux = c.mpv*iv';
    iy.iv = ny+(1:niv);
 iE(iy.iv) = false;
    ny = ny+niv;
  elseif strcmp(str,'ih')
    ds.ih = 'Halo currents';
    iy.ih = ny+(1:nih);
 iE(iy.ih) = false;
    ny = ny+nih;
  elseif strcmp(str,'sf')
    ds.sf = 'Spline vector for fpol^2/2';
    iy.sf = ny+(1:nsf);
 iE(iy.sf) = false;
    ny = ny+nsf;
  elseif strcmp(str,'sp')
    ds.sp = 'Spline vector for pres';
    iy.sp = ny+(1:nsp);
 iE(iy.sp) = false;
    ny = ny+nsp;
  elseif strcmp(str,'fl')
    ds.fl = 'Flux loop signals as defined by config fields mlc, mlv, mpl';
    if isfield(config,'mlc') & isfield(config,'mlv') & isfield(config,'mpl')
      iy.fl = ny+(1:nfl);
   iE(iy.fl) = true;
      ny = ny+nfl;
    else
      disp(['Warning gsconfig: ', ...
       'Failed to include flux loops (fl) among outputs! ', ...
       'At least one of fields mlc, mlv, mpl missing.'])
    end
  elseif strcmp(str,'lv')
    ds.lv = 'Loop voltage signals as defined by config fields mhc, mhv, mph';
    if isfield(config,'mhc') & isfield(config,'mhv') & isfield(config,'mph')
      iy.lv = ny+(1:nlv);
   iE(iy.lv) = false;
      ny = ny+nlv;
    else
      disp(['Warning gsconfig: ', ...
       'Failed to include flux loops for loop voltage (lv) among outputs! ', ...
       'At least one of fields mhc, mhv, mph missing.'])
    end
  elseif strcmp(str,'bp')
    ds.bp = 'Magnetic probe signals as defined by config fields gbc, gbv, gpb';
    if isfield(config,'gbc') & isfield(config,'gbv') & isfield(config,'gpb')
      iy.bp = ny+(1:nbp);
   iE(iy.bp) = true;
      ny = ny+nbp;
    else
      disp(['Warning gsconfig: ', ...
       'Failed to include magnetic probes (bp) in outputs! ', ...
       'At least one of fields gbc, gbv, gpb missing from config.'])
    end
  elseif strcmp(str,'rog')
    ds.rog = 'Rogowski signals as defined by field rldata';
    if isfield(config,'rldata')
      if size(config.rldata,1) == nc+nv+1 % This is what it should be
        rldata = config.rldata'; % The transpose is used in gs codes
      elseif size(config.rldata,2) == nc+nv+1 % Transposed? We'll take it!
        rldata = config.rldata;
      else
        rldata = zeros(0,nc+nv+1);
	disp(['Warning gsconfig: ', ....
	 'Failed to include Rogowski loops (rog) in outputs! ', ...
	 'The field rldata must have dimensions [nc+nv+1,nrl].'])
      end
      nrog = size(rldata,1);
      iy.rog = ny+(1:nrog);
   iE(iy.rog) = true;
      ny = ny+nrog;
    else
      disp(['Warning gsconfig: ', ....
       'Failed to include Rogowski loops (rog) among outputs! ', ...
       'The field rldata is missing.'])
    end
  elseif strcmp(str,'psidp')
    ds.psidp = 'Flux at diagnostic points rdp, zdp';
    if ndp > 0
      iy.psidp = ny+(1:ndp);
   iE(iy.psidp) = true;
      ny = ny+ndp;
    else
      disp(['Warning gsconfig: ', ...
       'Failed to include flux at diagnostic points (psidp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing.'])
    end
  elseif strcmp(str,'brdp')
    ds.brdp = 'Radial magnetic field at diagnostic points rdp, zdp';
    if ndp > 0
      iy.brdp = ny+(1:ndp);
   iE(iy.brdp) = true;
      ny = ny+ndp;
    else
      disp(['Warning gsconfig: ', ...
       'Failed to include radial magnetic field at diagnostic points (brdp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing.'])
    end
  elseif strcmp(str,'bzdp')
    ds.bzdp = 'Vertical magnetic field at diagnostic points rdp, zdp';
    if ndp > 0
      iy.bzdp = ny+(1:ndp);
   iE(iy.bzdp) = true;
      ny = ny+ndp;
    else
      disp(['Warning gsconfig: ', ...
       'Failed to include vertical magnetic field at diagnostic points (bzdp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing.'])
    end
  elseif strcmp(str,'btdp')
    ds.btdp = 'Toroidal magnetic field at diagnostic points rdp, zdp';
    if ndp > 0
      iy.btdp = ny+(1:ndp);
   iE(iy.btdp) = true;
      ny = ny+ndp;
    else
      disp(['Warning gsconfig: ', ...
       'Failed to include vertical magnetic field at diagnostic points (btdp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing.'])
    end
  elseif strcmp(str,'gapdp')
    ds.gapdp = 'Distance from rdp, zdp to boundary in direction toward axis';
    if ndp > 0
      iy.gapdp = ny+(1:ndp);
   iE(iy.gapdp) = true;
      ny = ny+ndp;
    else
      disp(['Warning gsconfig: ', ...
       'Failed to include gaps from rdp, zdp to boundary (gapdp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing.'])
    end
  elseif strcmp(str,'gaps')
    ds.gaps = 'Gaps from r,z along gr,gz to separatrix, gapspec=[r,z,gr,gz]';
    if ngaps > 0
      iy.gaps = ny+(1:ngaps);
   iE(iy.gaps) = true;
      ny = ny+ngaps;
    else
      disp(['Warning gsconfig: ', ...
       'Failed to include gaps among outputs! ', ...
       'No gaps were specified in field gapspec'])
    end
  elseif strcmp(str,'ys')
    ds.ys = 'Flux at all conductors';
    iy.ys = ny+(1:nc+nv);
 iE(iy.ys) = true;
    ny = ny+nc+nv;
  elseif strcmp(str,'rx')
    if     xpointorder == 1
      ds.rx = ['Radius for ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to normalized flux'];
    elseif xpointorder == 2
      ds.rx = ['Radius for ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to height'];
    elseif xpointorder == 3
      ds.rx = ['Radius for ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to normalized flux, then tracked'];
    elseif xpointorder == 4
      ds.rx = ['Radius for ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to height, then tracked'];
    end
    iy.rx = ny+(1:nxpoints);
 iE(iy.rx) = true;
    ny = ny+nxpoints;
  elseif strcmp(str,'zx')
    if     xpointorder == 1
      ds.zx = ['Height for ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to normalized flux'];
    elseif xpointorder == 2
      ds.zx = ['Height for ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to height'];
    elseif xpointorder == 3
      ds.zx = ['Height for ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to normalized flux, then tracked'];
    elseif xpointorder == 4
      ds.zx = ['Height for ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to height, then tracked'];
    end
    iy.zx = ny+(1:nxpoints);
 iE(iy.zx) = true;
    ny = ny+nxpoints;
  elseif strcmp(str,'psix')
    if     xpointorder == 1
      ds.psix = ['The flux at ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to normalized flux'];
    elseif xpointorder == 2
      ds.psix = ['The flux at ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to height'];
    elseif xpointorder == 3
      ds.psix = ['The flux at ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to normalized flux, then tracked'];
    elseif xpointorder == 4
      ds.psix = ['The flux at ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to height, then tracked'];
    end
    iy.psix = ny+(1:nxpoints);
 iE(iy.psix) = true;
    ny = ny+nxpoints;
  elseif strcmp(str,'rcur')
    ds.rcur = 'Radius of current centroid';
    iy.rcur = ny+1;
 iE(iy.rcur) = true;
    ny = ny+1;
  elseif strcmp(str,'zcur')
    ds.zcur = 'Height of current centroid';
    iy.zcur = ny+1;
 iE(iy.zcur) = true;
    ny = ny+1;
  elseif strcmp(str,'rcurdot')
    ds.rcurdot = 'Radial speed of current centroid';
    iy.rcurdot = ny+1;
 iE(iy.rcurdot) = true;
    ny = ny+1;
  elseif strcmp(str,'zcurdot')
    ds.zcurdot = 'Vertical speed of current centroid';
    iy.zcurdot = ny+1;
 iE(iy.zcurdot) = true;
    ny = ny+1;
  elseif strcmp(str,'cpasma')
    ds.cpasma = 'Total toroidal plasma current';
    iy.cpasma = ny+1;
 iE(iy.cpasma) = true;
    ny = ny+1;
  elseif strcmp(str,'Ip')
    ds.Ip = 'Total toroidal plasma current [MA]';
    iy.Ip = ny+1;
 iE(iy.Ip) = true;
    ny = ny+1;
  elseif strcmp(str,'aminor')
    ds.aminor = 'Horizontal minor radius of plasma';
    iy.aminor = ny+1;
 iE(iy.aminor) = true;
    ny = ny+1;
  elseif strcmp(str,'rb')
    ds.rb = 'Radius of boundary points';
    iy.rb = ny+(1:nb);
 iE(iy.rb) = false;
    ny = ny+nb;
  elseif strcmp(str,'zb')
    ds.zbbbs = 'Height of boundary points';
    iy.zb = ny+(1:nb);
 iE(iy.zb) = false;
    ny = ny+nb;
  elseif strcmp(str,'rmaxis')
    ds.rmaxis = 'Radius of magnetic axis';
    iy.rmaxis = ny+1;
 iE(iy.rmaxis) = true;
    ny = ny+1;
  elseif strcmp(str,'zmaxis')
    ds.zmaxis = 'Height of magnetic axis';
    iy.zmaxis = ny+1;
 iE(iy.zmaxis) = true;
    ny = ny+1;
  elseif strcmp(str,'psimag')
    ds.psimag = 'Flux at magnetic axis';
    iy.psimag = ny+1;
 iE(iy.psimag) = true;
    ny = ny+1;
  elseif strcmp(str,'rbdef')
    ds.rbdef = 'Radius of point that defines the boundary';
    iy.rbdef = ny+1;
 iE(iy.rbdef) = false;
    ny = ny+1;
  elseif strcmp(str,'zbdef')
    ds.zbdef = 'Height of point that defines the boundary';
    iy.zbdef = ny+1;
 iE(iy.zbdef) = false;
    ny = ny+1;
  elseif strcmp(str,'psibry')
    ds.psibry = 'Flux at boundary';
    iy.psibry = ny+1;
 iE(iy.psibry) = true;
    ny = ny+1;
  elseif strcmp(str,'isdiverted')
    ds.isdiverted = '-1 = vacuum, 0 = limited, 1 = diverted';
    iy.isdiverted = ny+1;
 iE(iy.isdiverted) = false;
    ny = ny+1;
  elseif strcmp(str,'li')
    ds.li = 'Normalized inductance';
    iy.li = ny+1;
 iE(iy.li) = true;
    ny = ny+1;
  elseif strcmp(str,'y2')
    ds.y2 = 'Second moment of the current distribution';
    iy.y2 = ny+1;
 iE(iy.y2) = true;
    ny = ny+1;
  elseif strcmp(str,'y2n')
    ds.y2n = 'Normalized second moment of the current distribution';
    iy.y2n = ny+1;
 iE(iy.y2n) = true;
    ny = ny+1;
  elseif strcmp(str,'betap')
    ds.betap = 'Poloidal beta';
    iy.betap = ny+1;
 iE(iy.betap) = true;
    ny = ny+1;
  elseif strcmp(str,'betan')
    ds.betan = 'Normalized toroidal beta';
    iy.betan = ny+1;
 iE(iy.betan) = true;
    ny = ny+1;
  elseif strcmp(str,'Wth')
    ds.Wth = 'Thermal energy';
    iy.Wth = ny+1;
 iE(iy.Wth) = true;
    ny = ny+1;
  elseif strcmp(str,'psihalo')
    ds.psihalo = 'Halo flux (current-weighted average)';
    iy.psihalo = ny+1;
 iE(iy.psihalo) = true;
    ny = ny+1;
  elseif strcmp(str,'psipla')
    ds.psipla = 'Plasma flux';
    iy.psipla = ny+1;
 iE(iy.psipla) = true;
    ny = ny+1;
  elseif strcmp(str,'dpsiplade')
    ds.dpsiplade = '-psipla response to flux error correction';
    iy.dpsiplade = ny+1;
 iE(iy.dpsiplade) = false;
    ny = ny+1;
  elseif strcmp(str,'bp2flx')
    ds.bp2flx = '(mu0*cpasma/Cl)^2';
    iy.bp2flx = ny+1;
 iE(iy.bp2flx) = true;
    ny = ny+1;
  elseif strcmp(str,'psiplaapp')
    ds.psiplaapp = 'Plasma flux from conductors';
    iy.psiplaapp = ny+1;
 iE(iy.psiplaapp) = true;
    ny = ny+1;
  elseif strcmp(str,'vplas')
    ds.vplas = 'Plasma voltage from conductors';
    iy.vplas = ny+1;
 iE(iy.vplas) = true;
    ny = ny+1;
  elseif strcmp(str,'Lpla')
    ds.Lpla = 'Plasma self inductance';
    iy.Lpla = ny+1;
 iE(iy.Lpla) = false;
    ny = ny+1;
  elseif strcmp(str,'lconn')
    ds.lconn = 'Representative connection length';
    iy.lconn = ny+1;
 iE(iy.lconn) = true;
    ny = ny+1;
  elseif strcmp(str,'R8')
    ds.R8 = 'R for points used to calculate shape parameters';
    iy.R8 = ny+[1:8];
 iE(iy.R8) = true;
    ny = ny+8;
  elseif strcmp(str,'Z8')
    ds.Z8 = 'Z for points used to calculate shape parameters';
    iy.Z8 = ny+[1:8];
 iE(iy.Z8) = true;
    ny = ny+8;
  elseif strcmp(str,'rsurf')
    ds.rsurf = 'R for geometric center';
    iy.rsurf = ny+1;
 iE(iy.rsurf) = true;
    ny = ny+1;
  elseif strcmp(str,'zsurf')
    ds.zsurf = 'Z for geometric center';
    iy.zsurf = ny+1;
 iE(iy.zsurf) = true;
    ny = ny+1;
  elseif strcmp(str,'aminor')
    ds.aminor = 'Half the radial width of the plasma';
    iy.aminor = ny+1;
 iE(iy.aminor) = true;
    ny = ny+1;
  elseif strcmp(str,'bminor')
    ds.bminor = 'Half the vertical height of the plasma';
    iy.bminor = ny+1;
 iE(iy.bminor) = true;
    ny = ny+1;
  elseif strcmp(str,'elong')
    ds.elong = 'Elongation = bminor/aminor';
    iy.elong = ny+1;
 iE(iy.elong) = true;
    ny = ny+1;
  elseif strcmp(str,'triu')
    ds.triu = 'Upper triangularity';
    iy.triu = ny+1;
 iE(iy.triu) = true;
    ny = ny+1;
  elseif strcmp(str,'tril')
    ds.tril = 'Lower triangularity';
    iy.tril = ny+1;
 iE(iy.tril) = true;
    ny = ny+1;
  elseif strcmp(str,'squo')
    ds.squo = 'Upper outer squareness';
    iy.squo = ny+1;
 iE(iy.squo) = true;
    ny = ny+1;
  elseif strcmp(str,'squi')
    ds.squi = 'Upper inner squareness';
    iy.squi = ny+1;
 iE(iy.squi) = true;
    ny = ny+1;
  elseif strcmp(str,'sqli')
    ds.sqli = 'Lower inner squareness';
    iy.sqli = ny+1;
 iE(iy.sqli) = true;
    ny = ny+1;
  elseif strcmp(str,'sqlo')
    ds.sqlo = 'Lower outer squareness';
    iy.sqlo = ny+1;
 iE(iy.sqlo) = true;
    ny = ny+1;
  elseif strcmp(str,'drsep')
    ds.drsep = '(flux difference to competing bdef)/(outboard midplane gradient)';
    iy.drsep = ny+1;
 iE(iy.drsep) = true;
    ny = ny+1;
  elseif strcmp(str,'Cl')
    ds.Cl = 'Contour length = length of the boundary';
    iy.Cl = ny+1;
 iE(iy.Cl) = true;
    ny = ny+1;
  elseif strcmp(str,'Vtot')
    ds.Vtot = 'Total plasma volume';
    iy.Vtot = ny+1;
 iE(iy.Vtot) = true;
    ny = ny+1;
  elseif strcmp(str,'Atot')
    ds.Atot = 'Area of plasma cross section';
    iy.Atot = ny+1;
 iE(iy.Atot) = true;
    ny = ny+1;
  elseif strcmp(str,'Ltot')
    ds.Ltot = 'Area-integral of 1/R';
    iy.Ltot = ny+1;
 iE(iy.Ltot) = true;
    ny = ny+1;
  elseif strcmp(str,'rhot')
    ds.rhot = 'Square root of normalized toroidal flux';
    iy.rhot = ny+(1:nr);
 iE(iy.rhot) = true;
    ny = ny+nr;
    calculate_profiles = 1;
    rhot = zeros(nr,1);
    drhotdx = zeros(nr,nx);
  elseif strcmp(str,'jtav')
    ds.jtav = 'Flux-surface averaged toroidal current density profile';
    iy.jtav = ny+(1:nr);
 iE(iy.jtav) = true;
    ny = ny+nr;
    calculate_profiles = 1;
    jtav = zeros(nr,1);
    djtavdx = zeros(nr,nx);
  elseif strcmp(str,'vres')
    ds.vres = 'Resistive voltage profile';
    iy.vres = ny+(1:np);
 iE(iy.vres) = false;
    ny = ny+np;
    calculate_helical_voltage = 1;
    vres = zeros(np,1);
    dvresdx = zeros(np,nx);
  elseif strcmp(str,'vind')
    ds.vind = 'Inductive voltage profile';
    iy.vind = ny+(1:np);
 iE(iy.vind) = false;
    ny = ny+np;
  elseif strcmp(str,'jpar')
    ds.jpar = 'Profile of parallel current';
    iy.jpar = ny+(1:np);
 iE(iy.jpar) = false;
    ny = ny+np;
  elseif strcmp(str,'psit')
    ds.psit = 'Profile of toroidal flux';
    iy.psit = ny+(1:np);
 iE(iy.psit) = false;
    ny = ny+np;
  elseif strcmp(str,'qpsi')
    ds.qpsi = 'Safety factor q';
    iy.qpsi = ny+(1:np);
 iE(iy.qpsi) = false;
    ny = ny+np;
  elseif strcmp(str,'L1t')
    ds.L1t = 'Average field line lengths through 1 toroidal turn';
    iy.L1t = ny+(1:np);
 iE(iy.L1t) = false;
    ny = ny+np;
  elseif strcmp(str,'fluxerror')
    ds.fluxerror = 'Error in normalized flux, < 1e-2 for valid Grad-Shafranov solutions';
    iy.fluxerror = ny+1;
 iE(iy.fluxerror) = false;
    ny = ny+1;
  elseif strcmp(str,'time')
    ds.time = 'Simulation time (used for plots and time-dependent configuration data)';
    iy.time = ny+1;
 iE(iy.time) = false;
    ny = ny+1;
  elseif strcmp(str,'gamma')
    ds.gamma = 'Growth rate of most unstable mode';
    iy.gamma = ny+1;
 iE(iy.gamma) = false;
    ny = ny+1;
  elseif strcmp(str,'drcurdv')
    ds.drcurdv = 'Response of rcur to most unstable eigen vector';
    iy.drcurdv = ny+1;
 iE(iy.drcurdv) = false;
    ny = ny+1;
  elseif strcmp(str,'dzcurdv')
    ds.dzcurdv = 'Response of zcur to most unstable eigen vector';
    iy.dzcurdv = ny+1;
 iE(iy.dzcurdv) = false;
    ny = ny+1;
  elseif strcmp(str,'nupdate')
    ds.nupdate = 'Number of calls to gsupdate';
    iy.nupdate = ny+1;
 iE(iy.nupdate) = false;
    ny = ny+1;
  elseif strcmp(str,'ncalc')
    ds.ncalc = 'Number of response calculation';
    iy.ncalc = ny+1;
 iE(iy.ncalc) = false;
    ny = ny+1;
  elseif strcmp(str,'frc')
    ds.frc = 'Radial forces on coils';
    iy.frc = ny+(1:ncoil);
 iE(iy.frc) = true;
    ny = ny+ncoil;
    need_mrmz = 1;
    frc = zeros(ncoil,1);
    dfrcdx = zeros(ncoil,nx);
  elseif strcmp(str,'fzc')
    ds.fzc = 'Vertical forces on coils';
    iy.fzc = ny+(1:ncoil);
 iE(iy.fzc) = true;
    ny = ny+ncoil;
    need_mrmz = 1;
    fzc = zeros(ncoil,1);
    dfzcdx = zeros(ncoil,nx);
  else 
    disp(['Warning gsconfig: output ' str ' not recognized'])
  end
end
if ny == 0
  ds.psizr = 'poloidal flux on the grid';
  iy.psizr = ny+(1:ngg);
  iE(iy.psizr) = true;
  ny = ny+ngg;
end
iy.descriptions = ds;

if need_mrmz
  if isfield(config,'mrcc') & isfield(config,'mzcc') & ...
     isfield(config,'mrpc') & isfield(config,'mzpc')
    mrcc = config.mrcc;
    mzcc = config.mzcc;
    mrpc = config.mrpc;
    mzpc = config.mzpc;
  elseif isfield(config,'fcdata') & isfield(config,'fcnturn')
    fcdata = config.fcdata;
    fcnturn = config.fcnturn;
    mrcc = zeros(nc);
    mzcc = zeros(nc);
    mrpc = zeros(ngg,nc);
    mzpc = zeros(ngg,nc);
    nfc = size(fcdata,2);
    nec = nc-nfc;
    for i = 1:nfc
      fprintf(...
'Creating objects for coil force calculations: coil %d of %d...',i+nec,nc)
      ri = fcdata(2,i)+[-1 -1 1 1 -1]/2*fcdata(4,i);
      zi = fcdata(1,i)+[-1 1 1 -1 -1]/2*fcdata(3,i);
      for j = 1:nfc
	rj = fcdata(2,j)+[-1 -1 1 1 -1]/2*fcdata(4,j);
	zj = fcdata(1,j)+[-1 1 1 -1 -1]/2*fcdata(3,j);
	[ms, mrs, mzs] = mpolygon2polygon(ri,zi,rj,zj);
	%MCC(i,j) = fcnturn(i)*fcnturn(j)*ms;
	mrcc(i+nec,j+nec) = fcnturn(i)*fcnturn(j)*mrs;
	mzcc(i+nec,j+nec) = fcnturn(i)*fcnturn(j)*mzs;
      end
      % For validation
      %[ms2, mrs2, mzs2] = mpolygon2point(ri,zi,rgg(:),zgg(:),2);
      %MPC(:,i) = fcnturn(i)*ms2;
      %GBR2C(:,i) = -fcnturn(i)*mzs2./rgg(:)/2/pi;
      %GBZ2C(:,i) = +fcnturn(i)*mrs2./rgg(:)/2/pi;
      [ms, mrs, mzs] = mpolygon2point(ri,zi,rgg(:),zgg(:),1);
      mrpc(:,i+nec) = fcnturn(i)*mzs;
      mzpc(:,i+nec) = fcnturn(i)*mrs;
      fprintf(char(8+zeros(1,80)))
    end
    fprintf('                                                               ')
    fprintf(char(8+zeros(1,80)))
  else
    disp('Warning gsconfig: mrcc, mzcc, mrpc, mzpc can not be calculated')
    disp('The fields fcdata and fcnturn are missing')
  end
end

lstarx = [mssx zeros(nic+niv,nih+nsf+nsp); zeros(nih+nsf+nsp,nx)];
Rhat = [rss+rssx zeros(nic+niv,nih+nsf+nsp); zeros(nih+nsf+nsp,nx)];
Vhat = [eye(nic) zeros(nic,nu-nic); ...  % Coil power supplies
        zeros(niv,nu); ...      % Vessel, zero voltage
        zeros(nih+nsf+nsp,nic) eye(nih+nsf+nsp,nu-nic)]; % halo, sf, sp

clear c d

d.tokamak = 'Name of tokamak';
c.tokamak = tokamak;

d.rg = 'R of grid points [m]';
c.rg = rg;
d.nr = 'Number of grid points radially';
c.nr = nr;
d.dr = 'Distance radially between grid points [m]';
c.dr = dr;

d.zg = 'Z of grid points [m]';
c.zg = zg;
d.nz = 'Number of grid points in Z-direction';
c.nz = nz;
d.dz = 'Distance in Z-direction between grid points [m]';
c.dz = dz;

d.rgg = 'R for all nz*nr grid points';
c.rgg = rgg;
d.zgg = 'Z for all nz*nr grid points';
c.zgg = zgg;
d.ngg = 'Number of grid points = nz*nr';
c.ngg = ngg;

d.Ag = 'Area of a grid cell [m2]';
c.Ag = Ag;
d.RA = 'Integral(R*dA) over cell';
c.RA = RA;
d.AR = 'Integral(dA/R) over cell';
c.AR = AR;

d.i16 = 'index differences to 16 grid points';
c.i16 = i16;

d.Rlim = 'R of limiter for points at corners';
c.Rlim = Rlim;
d.Zlim = 'Z of limiter for points at corners';
c.Zlim = Zlim;
d.limpoly = 'Variables used to determine if point is inside limiter';
c.limpoly = limpoly;
d.rl = 'R of limiter at corners and lines of a doubly dense grid';
c.rl = rl;
d.zl = 'Z of limiter at corners and lines of a doubly dense grid';
c.zl = zl;
d.nl = 'length(rl)';
c.nl = nl;
d.il = 'Indices used for interpolation to rl, zl';
c.il = il;
d.wl = 'Flux at rl,zl = sum(wl.*psizr(il),2)';
c.wl = wl;
d.wld1 = 'wld1(i,:) calculates derivative at i along vector from i to i+1';
c.wld1 = wld1;
d.wld2 = 'wld1(i,:) calculates derivative at i+1 along vector from i to i+1';
c.wld2 = wld2;
d.drl = 'diff(rl)';
c.drl = drl;
d.dzl = 'diff(zl)';
c.dzl = dzl;
d.dl = 'sqrt(drl.^2+dzl.^2)';
c.dl = dl;
d.irl = 'r-indices';
c.irl = irl;
d.trl = 'rl(1:nl-1)''-irl';
c.trl = trl;
d.tzl = 'zl(1:nl-1)''-izl';
c.tzl = tzl;
d.izl = 'z-indices';
c.izl = izl;
d.kl = 'Collapsed indices';
c.kl = kl;
d.turnin = 'Matrix to vector along limiter into machine';
c.turnin = turnin;
d.concavel = 'Limiter corner points into machine';
c.concavel = concavel;

d.Pcc = 'coil currents = Pcc * circuit currents';
c.Pcc = Pcc;
d.Pcci = 'circuit currents = Pcci * coil currents';
c.Pcci = Pcci;
d.Pvc = 'vessel currents = Pvc * vessel circuit currents';
c.Pvc = Pvc;
d.Pvci = 'vessel circuit currents = Pvci * vessel currents';
c.Pvci = Pvci;

d.constraints = 'Reduces nsp to 1 and nsf to 2 if > 0';
c.constraints = constraints;
d.psikf = 'ffprim/mu0 spline knot positions';
c.psikf = psikf;
d.nkf = 'length(psikf)-1';
c.nkf = nkf;
d.mf0 = 'ffprim/mu0 coefficients f0 = mf0*sf';
c.mf0 = mf0;
d.mf1 = 'ffprim/mu0 coefficients f1 = mf1*sf';
c.mf1 = mf1;
d.mf2 = 'ffprim/mu0 coefficients f2 = mf2*sf';
c.mf2 = mf2;
d.mf3 = 'ffprim/mu0 coefficients f3 = mf3*sf';
c.mf3 = mf3;
d.fb = 'fpol at boundary = fb + vf*sf';
c.fb = fb;
d.vf = 'fpol at boundary = fb + vf*sf';
c.vf = vf;
d.psikp = 'pprime spline knot positions';
c.psikp = psikp;
d.nkp = 'length(psikp)-1';
c.nkp = nkp;
d.mp0 = 'pprime coefficients p0 = mp0*sp';
c.mp0 = mp0;
d.mp1 = 'pprime coefficients p1 = mp1*sp';
c.mp1 = mp1;
d.mp2 = 'pprime coefficients p2 = mp2*sp';
c.mp2 = mp2;
d.mp3 = 'pprime coefficients p3 = mp3*sp';
c.mp3 = mp3;
d.pb = 'pres at boundary = pb + vp*sp';
c.pb = pb;
d.vp = 'pres at boundary = pb + vp*sp';
c.vp = vp;

d.defaultf = 'True if constrained ffprim/mu0 is the default';
c.defaultf = defaultf;
d.defaultp = 'True if constrained pprime is the default = ones(nr,1)';
c.defaultp = defaultp;

d.mcc = 'Mutuals between coil circuits and coil circuits';
c.mcc = mcc;
d.mcv = 'Mutuals between coil circuits and vessel circuits';
c.mcv = mcv;
d.mvv = 'Mutuals between vessel circuits and vessel circuits';
c.mvv = mvv;
d.mpc = 'Mutuals between grid points and coil circuits';
c.mpc = mpc;
d.mpv = 'Mutuals between grid points and vessel circuits';
c.mpv = mpv;
d.mgg = 'Mutuals between grid points and grid points';
c.mgg = mgg;

d.nfl = 'Number of flux loops';
c.nfl = nfl;
d.mlc = 'Mutuals between flux loops and coils';
c.mlc = mlc;
d.mlv = 'Mutuals between flux loops and vessel';
c.mlv = mlv;
d.mpl = 'Mutuals between flux loops and grid';
c.mpl = mpl;

d.nlv = 'Number of voltage loops';
c.nlv = nlv;
d.mhc = 'Mutuals between voltage loops and coils';
c.mhc = mhc;
d.mhv = 'Mutuals between voltage loops and vessel';
c.mhv = mhv;
d.mph = 'Mutuals between voltage loops and grid';
c.mph = mph;

d.nbp = 'Number of magnetic probes';
c.nbp = nbp;
d.gbc = 'Greens between magnetic probes and coils';
c.gbc = gbc;
d.gbv = 'Greens between magnetic probes and vessel';
c.gbv = gbv;
d.gpb = 'Greens between magnetic probes and grid';
c.gpb = gpb;

d.nrog = 'Number of Rogowski loops';
c.nrog = nrog;
d.rldata = 'Fraction of coil, vessel, plasma inside Rogowski loops';
c.rldata = rldata;

d.halo = '0=no halo, 1=scalar halo model';
c.halo = halo;
if halo > 0
% This field is excluded if halo == 0 because the compiler
% for matlab function blocks in Simulink is a piece of crap
% that must be tricked into doing the right thing
% In this case gsupdate must check if the field jhalo exists
% and include halo currents only if that is the case
d.jhalo = 'Current density used with halo = 1';
c.jhalo = jhalo;
end
d.Acell = 'Area inside limiter';
c.Acell = Acell;
d.RAcell = 'Integral[R*dA] inside limiter';
c.RAcell = RAcell;
d.ARcell = 'Integral[dA/R] inside limiter';
c.ARcell = ARcell;

d.nic = 'Number of circuits with coils (and a power supply)';
c.nic = nic;
d.niv = 'Number of circuits with vessel elements (no power supply)';
c.niv = niv;
d.nih = 'Number of states describing halo currents';
c.nih = nih;
d.nsf = 'Number of states for ffprim/mu0 and boundary fpol';
c.nsf = nsf;
d.nsp = 'Number of states for pprime and boundary pressure';
c.nsp = nsp;
d.nx = 'Total number of states = nic+niv+nih+nsp+nsf';
c.nx = nx;
d.nu = 'Number of actuators';
c.nu = nu;
d.ny = 'Number of outputs';
c.ny = ny;

d.ix = 'Contains indices to ic, iv, ih, sp, sf in state vector x';
c.ix = ix;
d.iu = 'Contains indices to ic, ih, sp, sf in actuator vector u';
c.iu = iu;
d.iy = 'Contains indices for diagnostics in output vector y';
c.iy = iy;
d.iE = 'Flag nonzero elements of E';
c.iE = iE;

d.nb = 'Number of different points on a contour';
c.nb = nb;
d.np = 'Number of contours';
c.np = np;
d.psibarp = 'Normalized flux for contours';
c.psibarp = linspace(0,1,np)';

d.nbmax = 'Fixed number of points returned by gscontour22 (padding is nans)';
c.nbmax = nbmax;
d.amin = 'Current exists outside last-closed-flux-surface if its aminor<amin';
c.amin = amin;
d.dpsibarlimit = 'Maximum change before updating response';
c.dpsibarlimit = dpsibarlimit;
d.nn = 'Max number of previous plasma responses remembered by gsupdate';
c.nn = nn;
d.nfla = 'Number of angles from axis to fl points';
c.nfla = nfla;
d.nflc = 'Number of contours with fl points';
c.nflc = nflc;
d.nfltest = 'Number of fluxes checked for validity of linear approximation';
c.nfltest = nfltest;
d.nbdtest = 'Number of points checked for jump in boundary-defining point';
c.nbdtest = nbdtest;

d.evolve_option = 'Option for how to evolve equilibrium';
c.evolve_option = evolve_option;
d.lstarx = 'Inductance in cables between power supplies and coils';
c.lstarx = lstarx;
d.Rhat = 'Generalized resistance in evolution equation';
c.Rhat = Rhat;
d.Vhat = 'Generalized voltage in evolution equation';
c.Vhat = Vhat;

if isfield(config,'ecdata')
  d.ecdata = 'E coil geometry';
  c.ecdata = config.ecdata;
end
if isfield(config,'fcdata')
  d.fcdata = 'F coil geometry';
  c.fcdata = config.fcdata;
end
if isfield(config,'vvdata')
  d.vvdata = 'Vessel element geometry';
  c.vvdata = config.vvdata;
end

d.plots = 'Settings for plots';
c.plots.tplot = [];
c.plots.dtplot = 0;
c.plots.updates = false;

d.now = 'Time stamp used as identity for this configuration';
c.now = now;

c.info = d;
