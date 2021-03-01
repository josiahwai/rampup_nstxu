function b = gsboundary(c,psizr,ra,za,dnx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   b = gsboundary(c,psizr,ra,za,dnx);
%
%  PURPOSE: Find boundary and related quantities needed to 
%           calculate flux from plasma
%
%  INPUTS:  c, configuration data made by gsconfig.m
%           psizr, flux on the grid c.rg, c.zg
%           ra, za, [m] pick axis closest to ra, za, if several
%           dnx, use nearx method when contour is within dnx of
%             an x-point [max=1, default=0.5 floating index units]
%           Nearx keeps boundary speed finite near x-points using
%           an interpolation that gives nonzero field at x-points
%
%  OUTPUTS: b, information about plasma boundary

%  RESTRICTIONS: Interpolation between traced boundary points is
%    inaccurate if boundary direction turns a lot between points
%    and fails where turn > 90 degrees, issues are unlikely when
%    plasma covers more than 4 grid cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Anders Welander ON 2016-11-13
%
%  MODIFICATION HISTORY:
%  Discovered Oct 13 2017 that the numerical landmines that cause
%  lstari in gsupdate to explode with high eigenvalues are caused
%  by boundary points in close proximity to the x-point. This also 
%  affects cond(eye(ngg)-jg). The solution is the NEARX method
%  that addresses both this issue and the issue with making double
%  equilibria in gsdesign.
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinates ra, za (if given) help to locate the right axis
if nargin < 3
  ra = nan;
end
if nargin < 4
  za = nan;
end

% dnx is a selected distance from boundary to x-points where nearx treatment kicks in
% Modify the boundary where x-points are within the distance dnx
if nargin < 5
  dnx = 0.5; % floating points units
end

% For cubic Hermite interpolation
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;

% Grid
nz = c.nz;
nr = c.nr;
i16 = c.i16; % index differences to 4x4 points used for interpolation
ngg = c.ngg; % nr*nz
RA = c.RA;
Ag = c.Ag;
AR = c.AR;
rg = c.rg;
dr = c.dr;
zg = c.zg;
dz = c.dz;

% Number of traced boundary points
nbmax = c.nbmax;

% Number of possible boundary-defining points to save
nbdtest = c.nbdtest;

% Used to quickly determine if point is inside limiter
limpoly = c.limpoly;

% Approximate location of axis
ra = (ra-rg(1))/dr+1; % Switch to grid cell units
za = (za-zg(1))/dz+1;
if isnan(ra)
  ra = nr/2;
end
if isnan(za)
  za = nz/2;
end

% For response calculations
S = ones(1,16); % Sixteen
T = ones(1,36); % Thirty six
F = ones(1,49); % Forty nine

% The 4x4 central indices in 6x6 matrix of indices
ic = [8 9 10 11 14 15 16 17 20 21 22 23 26 27 28 29];

% index shifts to 7x7 cells around an origin
di7x7 = reshape((-3:3)'*ones(1,7)+ones(7,1)*nz*(-3:3),49,1);

% Info about the boundary-defining (bdef) point and points that can become bdef
bd.n = 0;
bd.r = nan(1,nbdtest);
bd.z = nan(1,nbdtest);
bd.w = nan(nbdtest,16);
bd.ii = ones(nbdtest,16);
bd.psi = nan(1,nbdtest);
bd.lim = zeros(1,nbdtest);
bd.wrr = zeros(nbdtest,16);
bd.wrz = zeros(nbdtest,16);
bd.wzz = zeros(nbdtest,16);
bd.drdpsi = zeros(nbdtest,16);
bd.dzdpsi = zeros(nbdtest,16);
bd.info.n = 'number of potential boundary-defining points';
bd.info.r = 'r of boundary-defining point and n-1 positions it can jump';
bd.info.z = 'z of boundary-defining point and n-1 positions it can jump';
bd.info.w = 'weights such that sum(w.*psizr(ii),2) = flux at points r, z';
bd.info.ii = 'indices such that sum(w.*psizr(ii),2) = flux at points r, z';
bd.info.psi = 'flux at points r, z';
bd.info.lim = '0 for nulls, touch points between lim:lim+1 in c.rl,c.zl';
bd.info.wrr = 'sum(wrr.*psizr(ii),2) = d(flux)/dr^2 at points r, z';
bd.info.wrz = 'sum(wrz.*psizr(ii),2) = d(flux)/dr/dz at points r, z';
bd.info.wzz = 'sum(wzz.*psizr(ii),2) = d(flux)/dz^2 at points r, z';
bd.info.drdpsi = 'r response to psizr, dr = sum(drdpsi.*dpsizr(ii),2)';
bd.info.dzdpsi = 'z response to psizr, dz = sum(dzdpsi.*dpsizr(ii),2)';


% For interpolation
xs = linspace(0,1,nr);

% Contour with nr-2 interpolated points between contoured points
R = nan((nr-1)*nbmax+1,1);
Z = nan((nr-1)*nbmax+1,1);

% Variables of size nbmax are uppercase (except l)

% 1:st index of cell with boundary segment
IC = ones(nbmax,1);

% 2:nd index of cell with boundary segment
JC = ones(nbmax,1);

% *b is response of * to flux change at the boundary-defining point
% *w is response of * to flux change in 6x6 surrounding grid points
% Example: dR0 = sum(R0w.*dpsizr(I36),2) + R0b*dpsibry

% Indices for use with variables of size [nbmax,36]
I36 = ones(nbmax,36);

% zeros and ones for use with nbmax-sized variables
O = zeros(nbmax,1);
l = ones(nbmax,1);

Yrw = zeros(nbmax,36);
Yzw = zeros(nbmax,36);
DRw = zeros(nbmax,36);
DZw = zeros(nbmax,36);
T0w = zeros(nbmax,36);
T1w = zeros(nbmax,36);
dCw = zeros(nbmax,36);

Yrb = zeros(nbmax,1);
Yzb = zeros(nbmax,1);
DRb = zeros(nbmax,1);
DZb = zeros(nbmax,1);
T0b = zeros(nbmax,1);
T1b = zeros(nbmax,1);
dCb = zeros(nbmax,1);

% Contoured points R0, Z0
R0 = nan(nbmax,1); % R of contour points
Z0 = nan(nbmax,1); % Z of contour points
R0b = zeros(nbmax,1);
Z0b = zeros(nbmax,1);
R0w = zeros(nbmax,36);
Z0w = zeros(nbmax,36);

Yr = nan(nbmax,1); % d(flux)/dr for points R0, Z0
Yz = nan(nbmax,1); % d(flux)/dz for points R0, Z0
Yrr = nan(nbmax,1); % d2(flux)/dr2 for points R0, Z0
Yrz = nan(nbmax,1); % d2(flux)/drdz for points R0, Z0
Yzz = nan(nbmax,1); % d2(flux)/dz2 for points R0, Z0

DR = nan(nbmax,1); % diff(R([1:end 1]))
DZ = nan(nbmax,1); % diff(Z([1:end 1]))

% for interpolation between contoured points
R1 = nan(nbmax,1);
Z1 = nan(nbmax,1);
R2 = nan(nbmax,1);
Z2 = nan(nbmax,1);
R3 = nan(nbmax,1);
Z3 = nan(nbmax,1);
R1b = zeros(nbmax,1);
R2b = zeros(nbmax,1);
R3b = zeros(nbmax,1);
Z1b = zeros(nbmax,1);
Z2b = zeros(nbmax,1);
Z3b = zeros(nbmax,1);
R1w = zeros(nbmax,36);
R2w = zeros(nbmax,36);
R3w = zeros(nbmax,36);
Z1w = zeros(nbmax,36);
Z2w = zeros(nbmax,36);
Z3w = zeros(nbmax,36);

% Boundary direction, 0=after the point, 1=before the point, going counter-clockwise
Br0 = nan(nbmax,1);
Bz0 = nan(nbmax,1);
Br1 = nan(nbmax,1);
Bz1 = nan(nbmax,1);
Br0b = nan(nbmax,1);
Bz0b = nan(nbmax,1);
Br1b = nan(nbmax,1);
Bz1b = nan(nbmax,1);
Br0w = nan(nbmax,36);
Bz0w = nan(nbmax,36);
Br1w = nan(nbmax,36);
Bz1w = nan(nbmax,36);

T0 = nan(nbmax,1); % Field line deviation from DR,DZ at [1:end]
T1 = nan(nbmax,1); % Field line deviation from DR,DZ at [2:end 1]
H = false(nbmax,1); % flags points on horizontal grid lines
V = false(nbmax,1); % flags points on vertical grid lines 
N = false(nbmax,1);% flags points rn, zn
NX = false(nbmax,1);% flags nearx points

dC = zeros(nbmax,1);

drsep = nan;
Cl = 0;
Vtot = 0;
Atot = 0;
Ltot = 0;
drsepp = zeros(1,ngg);
Clp = zeros(1,ngg);
Vtotp = zeros(1,ngg);
Atotp = zeros(1,ngg);
Ltotp = zeros(1,ngg);

% Allocate memory for [nz,nr] sized variables
Acell = zeros(nz,nr);  % plasma-covered area of cells
RAcell = zeros(nz,nr); % Integral(R*dA) over plasma within cells
ZAcell = zeros(nz,nr); % Integral(Z*dA) over plasma within cells
ARcell = zeros(nz,nr); % Integral(dA/R) over plasma within cells

Acellb = zeros(ngg,1);
Acellw = zeros(ngg,49);
RAcellb = zeros(ngg,1);
RAcellw = zeros(ngg,49);
ZAcellb = zeros(ngg,1);
ZAcellw = zeros(ngg,49);
ARcellb = zeros(ngg,1);
ARcellw = zeros(ngg,49);

rgc = ones(nz,1)*rg;   % R for center of current within cell
zgc = zg*ones(1,nr);   % Z for center of current within cell
rgcb = zeros(ngg,1);
zgcb = zeros(ngg,1);
rgcw = zeros(ngg,49);
zgcw = zeros(ngg,49);

igg = false(nz,nr); % will be 1 for grid points covered by plasma
egg = false(nz,nr); % will be 1 for grid cells cut by boundary
% psizrc is flux at rgc,zgc, also calculating derivatives w.r.t. r,z
psizrc = psizr;
psizrcr = [-psizr(:,3)+4*psizr(:,2)-3*psizr(:,1) ...
           psizr(:,3:end)-psizr(:,1:end-2) ...
	   psizr(:,end-2)-4*psizr(:,end-1)+3*psizr(:,end)]/2/dr;
psizrcz = [-psizr(3,:)+4*psizr(2,:)-3*psizr(1,:); ...
           psizr(3:end,:)-psizr(1:end-2,:); ...
	   psizr(end-2,:)-4*psizr(end-1,:)+3*psizr(end,:)]/2/dz;
psizrcrr = zeros(nz,nr); % Will hold d2(psizrc)/drdr for egg
psizrcrz = zeros(nz,nr); % Will hold d2(psizrc)/drdz for egg
psizrczz = zeros(nz,nr); % Will hold d2(psizrc)/dzdz for egg

psizrcw = zeros(ngg,49);
psizrcw(:,25) = 1;

psizrcb = zeros(ngg,1);

psizrcrw = zeros(ngg,49);
psizrcrw(:,18) = -1/dr/2;
psizrcrw(:,32) = 1/dr/2;

psizrcrb = zeros(ngg,1);

psizrczw = zeros(ngg,49);
psizrczw(:,24) = -1/dz/2;
psizrczw(:,26) = 1/dz/2;

psizrczb = zeros(ngg,1);
  
dA0 = zeros(nbmax,1);
dRA = zeros(nbmax,1);
dZA = zeros(nbmax,1);
dAR = zeros(nbmax,1);
  
dA0w = zeros(nbmax,36);
dRAw = zeros(nbmax,36);
dZAw = zeros(nbmax,36);
dARw = zeros(nbmax,36);

dA0b = zeros(nbmax,1);
dRAb = zeros(nbmax,1);
dZAb = zeros(nbmax,1);
dARb = zeros(nbmax,1);

dA0e = zeros(nbmax,1);
dRAe = zeros(nbmax,1);
dZAe = zeros(nbmax,1);
dARe = zeros(nbmax,1);

dA0ew = zeros(nbmax,36);
dRAew = zeros(nbmax,36);
dZAew = zeros(nbmax,36);
dARew = zeros(nbmax,36);

dA0eb = zeros(nbmax,1);
dRAeb = zeros(nbmax,1);
dZAeb = zeros(nbmax,1);
dAReb = zeros(nbmax,1);

% Shape parameters, initialized with vacuum solution
R8 = nan(8,1); % points begin with outboard midplane
Z8 = nan(8,1); % walk around in 45 degree steps
rsurf = nan; % R for the geometric center
zsurf = nan; % Z for the geometric center
aminor = 0;  % half of the radial width of the plasma
bminor = 0;  % half of the height of the plasma
elong = nan; % elongation = bminor/aminor
tril = nan;  % Lower triangularity
triu = nan;  % Upper triangularity
squo = nan;  % Upper outer squareness
squi = nan;  % Upper inner squareness
sqli = nan;  % Lower inner squareness
sqlo = nan;  % Lower outer squareness

% Response of shape parameters to psizr, vacuum solution
R8p = zeros(8,ngg);
Z8p = zeros(8,ngg);
rsurfp = zeros(1,ngg);
zsurfp = zeros(1,ngg);
aminorp = zeros(1,ngg);
bminorp = zeros(1,ngg);
elongp = zeros(1,ngg);
trilp = zeros(1,ngg);
triup = zeros(1,ngg);
squop = zeros(1,ngg);
squip = zeros(1,ngg);
sqlip = zeros(1,ngg);
sqlop = zeros(1,ngg);

% Find nulls in the flux
nulls = gsnulls(psizr,4*nbdtest);

% Find points where plasma possibly touches the limiter
touches = gstouches(c,psizr);

% Find magnetic axis and flag nulls inside the limiter
wa = nan(1,16);
warr = nan(1,16);
warz = nan(1,16);
wazz = nan(1,16);
iia = ones(16,1);
rmaxis = nan;
zmaxis = nan;
emaxis = 1;
kmaxis = 1;
psimag = nan;
drmaxisdpsi = zeros(1,16);
dzmaxisdpsi = zeros(1,16);
curdir = 0;
invessel = true(1,nbdtest);
d2min = inf;
for i = 1:nulls.count
  r = (nulls.r(i)-1)*dr+rg(1);
  z = (nulls.z(i)-1)*dz+zg(1);
  f = false;
  for j = 1:limpoly.unique.n-1
    if limpoly.unique.z(j) < z & z <= limpoly.unique.z(j+1)
      for k = 1:limpoly.nranges
        d = z - limpoly.unique.z(j);
        rmin = limpoly.rranges(j,2*k-1)+limpoly.drdz(j,2*k-1)*d;
        rmax = limpoly.rranges(j,2*k  )+limpoly.drdz(j,2*k  )*d;
	f = f | rmin < r & r < rmax;
      end
    end
  end
  invessel(i) = f;
  d2 = (nulls.r(i)-ra)^2+(nulls.z(i)-za)^2;
  if f & nulls.type(i) == 'O' & d2 < d2min
    d2min = d2;
    wa = nulls.w(i,:);
    warr = nulls.wrr(i,:)/dr/dr;
    warz = nulls.wrz(i,:)/dr/dz;
    wazz = nulls.wzz(i,:)/dz/dz;
    iia = nulls.ii(i,:)';
    drmaxisdpsi = dr*nulls.drdy(i,:);
    dzmaxisdpsi = dz*nulls.dzdy(i,:);
    rmaxis = nulls.r(i);
    zmaxis = nulls.z(i);
    ur = nulls.ur(i);
    uz = nulls.uz(i);
    vr = +uz;
    vz = -ur;
    yarr = warr*psizr(iia);
    yarz = warz*psizr(iia);
    yazz = wazz*psizr(iia);
    % Find principal axis for physics units (gsnulls has it for grid units)
    u = 0;
    for k = 1:3
      ru = cos(u);
      zu = sin(u);
      D3u = (yazz-yarr)*ru*zu+yarz*(ru^2-zu^2); % d(d2)/dv should be 0
      D4u = (yazz-yarr)*(ru^2-zu^2)-4*yarz*ru*zu;
      u = u-D3u/D4u;
    end
    yauu = yarr*ru^2+2*yarz*ru*zu+yazz*zu^2;
    yavv = yarr*zu^2-2*yarz*ru*zu+yazz*ru^2;
    emaxis = sqrt(abs(yarr/yazz));
    kmaxis = sqrt(abs(yauu/yavv));
    psimag = nulls.y(i);
    curdir = -sign(nulls.yuu(i));
  end
end

% Memory allocation for possible boundary-defining points
npbd = max(200,8*nbdtest);
pbd.n = 0;
pbd.r = nan(1,npbd);
pbd.z = nan(1,npbd);
pbd.w = nan(npbd,16);
pbd.ii = ones(npbd,16);
pbd.psi = nan(1,npbd);
pbd.lim = zeros(1,npbd);
pbd.wrr = zeros(npbd,16);
pbd.wrz = zeros(npbd,16);
pbd.wzz = zeros(npbd,16);
pbd.wrrr = zeros(npbd,16);
pbd.wrrz = zeros(npbd,16);
pbd.wrzz = zeros(npbd,16);
pbd.wzzz = zeros(npbd,16);
pbd.drdpsi = zeros(npbd,16);
pbd.dzdpsi = zeros(npbd,16);
pbd.priority = inf(1,npbd);

% Store nulls that can potentially define the boundary in pbd
for i = 1:nulls.count
  if invessel(i)                    & ...
     nulls.type(i) == 'X'           & ...
     curdir*(psimag-nulls.y(i)) > 0 & ...
     (nulls.r(i)-rmaxis)^2+(nulls.z(i)-zmaxis)^2 > 0.25 & ...
     pbd.n < npbd
    pbd.n = pbd.n+1;
    pbd.r(pbd.n) = nulls.r(i);
    pbd.z(pbd.n) = nulls.z(i);
    pbd.w(pbd.n,:) = nulls.w(i,:);
    pbd.ii(pbd.n,:) = nulls.ii(i,:);
    pbd.psi(pbd.n) = nulls.y(i);
    pbd.wrr(pbd.n,:) = nulls.wrr(i,:);
    pbd.wrz(pbd.n,:) = nulls.wrz(i,:);
    pbd.wzz(pbd.n,:) = nulls.wzz(i,:);
    pbd.wrrr(pbd.n,:) = nulls.wrrr(i,:);
    pbd.wrrz(pbd.n,:) = nulls.wrrz(i,:);
    pbd.wrzz(pbd.n,:) = nulls.wrzz(i,:);
    pbd.wzzz(pbd.n,:) = nulls.wzzz(i,:);
    pbd.drdpsi(pbd.n,:) = nulls.drdy(i,:);
    pbd.dzdpsi(pbd.n,:) = nulls.dzdy(i,:);
    if (nulls.r(i)-rmaxis)^2+(nulls.z(i)-zmaxis)^2 > 25
      pbd.priority(pbd.n) = -curdir*nulls.y(i);
    else
      pbd.priority(pbd.n) = 1e9;
    end
  end
end

% Store possible touch points in pbd
for i = 1:touches.count
  if pbd.n < npbd & curdir*(psimag-touches.psi(i)) > 0
    pbd.n = pbd.n+1;
    pbd.r(pbd.n) = touches.r(i);
    pbd.z(pbd.n) = touches.z(i);
    pbd.w(pbd.n,:) = touches.w(i,:);
    pbd.ii(pbd.n,:) = touches.ii(i,:);
    pbd.psi(pbd.n) = touches.psi(i);
    pbd.lim(pbd.n) = touches.limiter(i);
    pbd.drdpsi(pbd.n,:) = touches.drdpsi(i,:);
    pbd.dzdpsi(pbd.n,:) = touches.dzdpsi(i,:);
    pbd.priority(pbd.n) = -curdir*touches.psi(i);
  end
end

% Sort the points to begin with point that has flux closest to psimag
[~,kpbd] = sort(pbd.priority);
pbd.r = pbd.r(kpbd);
pbd.z = pbd.z(kpbd);
pbd.w = pbd.w(kpbd,:);
pbd.ii = pbd.ii(kpbd,:);
pbd.psi = pbd.psi(kpbd);
pbd.lim = pbd.lim(kpbd);
pbd.wrr = pbd.wrr(kpbd,:);
pbd.wrz = pbd.wrz(kpbd,:);
pbd.wzz = pbd.wzz(kpbd,:);
pbd.wrrr = pbd.wrrr(kpbd,:);
pbd.wrrz = pbd.wrrz(kpbd,:);
pbd.wrzz = pbd.wrzz(kpbd,:);
pbd.wzzz = pbd.wzzz(kpbd,:);
pbd.drdpsi = pbd.drdpsi(kpbd,:);
pbd.dzdpsi = pbd.dzdpsi(kpbd,:);

% Try to make a closed contour using the points in pbd as starting points
% It has to wrap around the magnetic axis ra, za and be inside the vessel
% When successful the boundary has been found
n = 0; % number of boundary points, value is updated when boundary found below
nn = 0; % Number of points in R, Z, updated if boundary exists
ipbd = 1; % index in pbd that holds boundary-defining point, updated below
lim = 0; % index in c.rl,c.zl that boundary touches, 0 means no touching, updated below
boundary_within_vessel = false; % if gscontour22 fails this flag won't be set
rbdef = nan;
zbdef = nan;
if ~isnan(rmaxis)
  % A magnetic axis exists and therefore a boundary must exist
  % (but it may be small and gscontour22 might still fail in unforeseen special cases)
  for i = 1:pbd.n
    ipbd = i;
    rbdef = pbd.r(ipbd);
    zbdef = pbd.z(ipbd);
    lim = pbd.lim(ipbd);
    [R0,Z0,n] = gscontour22(psizr,rbdef,zbdef,nbmax,rmaxis,zmaxis);
    if ~isnan(R0(1))
      % Check that the boundary is within the vessel
      boundary_within_vessel = true;
      for m = 2:n-1
	r = (R0(m)-1)*dr+rg(1);
	z = (Z0(m)-1)*dz+zg(1);
	f = false;
	for j = 1:limpoly.unique.n-1
	  if limpoly.unique.z(j) < z & z <= limpoly.unique.z(j+1)
	    for k = 1:limpoly.nranges
              d = z - limpoly.unique.z(j);
              rmin = limpoly.rranges(j,2*k-1)+limpoly.drdz(j,2*k-1)*d;
              rmax = limpoly.rranges(j,2*k  )+limpoly.drdz(j,2*k  )*d;
	      % Add +dr/1e9 below to ensure test works when corner touches plasma
	      f = f | rmin < r+dr/1e9 & r < rmax+dr/1e9;
	    end
	  end
	end
	% f is true if point m is inside vessel
	boundary_within_vessel = boundary_within_vessel & f;
      end
      % Make R, Z always go counter-clockwise
      if sum((R0(1:n-1)+R0(2:n)).*(Z0(2:n)-Z0(1:n-1))) < 0
	R0(1:n) = R0(n:-1:1);
	Z0(1:n) = Z0(n:-1:1);
      end
      if boundary_within_vessel
        break
      end
    end
  end
end

% Reject points 1:ipbd-1
% The flux can not be changing monotonically from those points to the axis 
% or they are outside the limiter
% This means they can not become boundary-defining by slight change in psizr
% Hence they are rejected
if ipbd > 1
  kpbd = [ipbd:npbd npbd+zeros(1,ipbd-1)];
  pbd.n = pbd.n-ipbd+1;
  pbd.r = pbd.r(kpbd);
  pbd.z = pbd.z(kpbd);
  pbd.w = pbd.w(kpbd,:);
  pbd.ii = pbd.ii(kpbd,:);
  pbd.psi = pbd.psi(kpbd);
  pbd.lim = pbd.lim(kpbd);
  pbd.wrr = pbd.wrr(kpbd,:);
  pbd.wrz = pbd.wrz(kpbd,:);
  pbd.wzz = pbd.wzz(kpbd,:);
  pbd.wrrr = pbd.wrrr(kpbd,:);
  pbd.wrrz = pbd.wrrz(kpbd,:);
  pbd.wrzz = pbd.wrzz(kpbd,:);
  pbd.wzzz = pbd.wzzz(kpbd,:);
  pbd.drdpsi = pbd.drdpsi(kpbd,:);
  pbd.dzdpsi = pbd.dzdpsi(kpbd,:);
  ipbd = 1;
end

% bd contains the nbdtest most important points from pbd
% and can be used to check when the bdef point will jump to new place
bd.n = min(nbdtest,pbd.n);
bd.r = pbd.r(1:nbdtest);
bd.z = pbd.z(1:nbdtest);
bd.w = pbd.w(1:nbdtest,:);
bd.ii = pbd.ii(1:nbdtest,:);
bd.psi = pbd.psi(1:nbdtest);
bd.lim = pbd.lim(1:nbdtest);
bd.wrr = pbd.wrr(1:nbdtest,:);
bd.wrz = pbd.wrz(1:nbdtest,:);
bd.wzz = pbd.wzz(1:nbdtest,:);
bd.wrrr = pbd.wrrr(1:nbdtest,:);
bd.wrrz = pbd.wrrz(1:nbdtest,:);
bd.wrzz = pbd.wrzz(1:nbdtest,:);
bd.wzzz = pbd.wzzz(1:nbdtest,:);
bd.drdpsi = pbd.drdpsi(1:nbdtest,:);
bd.dzdpsi = pbd.dzdpsi(1:nbdtest,:);

% Flux at boundary-defining point and hence all of boundary
wb = bd.w(1,:);
iib = bd.ii(1,:)';
wbrr = bd.wrr(1,:);
wbrz = bd.wrz(1,:);
wbzz = bd.wzz(1,:);
wbrrr = bd.wrrr(1,:);
wbrrz = bd.wrrz(1,:);
wbrzz = bd.wrzz(1,:);
wbzzz = bd.wzzz(1,:);
psibry = bd.psi(1);
drbdefdpsi = bd.drdpsi(1,:);
dzbdefdpsi = bd.dzdpsi(1,:);

% A plasma exists if there is an axis and boundary within the vessel and the area is big enough
plasma = ~isnan(rmaxis) & ...
         ~isnan(R0(1)) & ...
	 boundary_within_vessel & ...
	 polyarea(R0(1:n),Z0(1:n)) > 4; % Methods here require plasma to cover at least 4 cells

% Indices to point 3,3 in 6x6 matrix of points around each of R0,Z0
IR = max(3,min(nr-3,floor(R0))); % J = floor(R0)
IZ = max(3,min(nz-3,floor(Z0))); % I = floor(Z0)

s = sign(psimag-psibry);

if plasma

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                              NEARX METHOD                                 %

  % Modify the boundary near x-points to prevent the speed from going infinite
  % when the boundary is close to the x-point. This method also allows a more
  % accurate treatment of double-null plasmas since extra boundary points are
  % inserted near all x-points, not just the one that defines the boundary

  % Flag x-points within dnx of boundary
  nearxregion = nulls.type == 'X'; % initial value

  % If x-point regions overlap then only the closest is nearx if within dnx of boundary
  d2x = (nulls.r-rmaxis).^2+(nulls.z-zmaxis).^2;
  isbdef = nulls.r == rbdef & nulls.z == zbdef;  
  for i = 1:nulls.count
    for j = i+1:nulls.count
      if nearxregion(i) & nearxregion(j) & ...
	(nulls.r(i)-nulls.r(j))^2+(nulls.z(i)-nulls.z(j))^2 < 4
	if d2x(i) > d2x(j) & ~isbdef(i) | isbdef(j)
	  nearxregion(i) = 0;
	else
	  nearxregion(j) = 0;
	end
      end
    end
  end

  % Check all nulls to find those that are close to the boundary
  for i = 1:nulls.count

    if nearxregion(i) % This is an x-point

      % position and flux
      rx = nulls.r(i);
      zx = nulls.z(i);
      psix = nulls.y(i);

      % R1, R2, R3 have temporary use in this for loop
      R2 = (R0-rx).^2+(Z0-zx).^2;
      R1 = sqrt(R2);
      [R3,K] = sort(R1);

      if R3(1) < 2*dnx % This x-point is near enough that it might be a "nearx" point

	% The perimeter point rp, zp is the point on a circle around the x-point
	% where the normalized flux is minimized and most inside or nearest the plasma

	% initial value of unit vector pointing toward rp, zp
	up = [rmaxis-rx, zmaxis-zx]/sqrt((rx-rmaxis)^2+(zx-zmaxis)^2);

	% Rotate vector 'up' until rp, zp is at the minimum in normalized flux
	da = pi/16; % initial angle increment
	nflip = 0; % number of times the rotation has flipped direction
	for k = 1:32
	  up = [up(1)*cos(da)-sin(da)*up(2), sin(da)*up(1)+cos(da)*up(2)];
	  ub = [-up(2), up(1)];
	  rp = rx + dnx*up(1);
	  zp = zx + dnx*up(2);
	  irp = floor(rp);
	  izp = floor(zp);
	  tr = rp-irp;
	  tz = zp-izp;
	  y44 = psizr(izp-1:izp+2,irp-1:irp+2);
	  wr0 = [1 tr tr^2 tr^3]*mx;
	  wz0 = [1 tz tz^2 tz^3]*mx;
	  wr1 = [0 1 2*tr 3*tr^2]*mx;
	  wz1 = [0 1 2*tz 3*tz^2]*mx;
	  wr2 = [0 0 2 6*tr]*mx;
	  wz2 = [0 0 2 6*tz]*mx;
	  psip = wz0*y44*wr0';
	  psipr = wz0*y44*wr1';
	  psipz = wz1*y44*wr0';
	  psipa = ub*[psipr;psipz]*dnx;
	  psiprr = wz0*y44*wr2';
	  psiprz = wz1*y44*wr1';
	  psipzz = wz2*y44*wr0';
	  if s*psipa*da < 0 % Heading the wrong way
	    da = -da/2;
	    nflip = nflip+1;
	  end
	  if nflip > 2
	    % Time for Newton-Rhapson to make psipa = 0
	    upa = [-up(1)*sin(da)-cos(da)*up(2), cos(da)*up(1)-sin(da)*up(2)];
	    uba = [-upa(2), upa(1)];
	    psipar = ub*[psiprr;psiprz]*dnx;
	    psipaz = ub*[psiprz;psipzz]*dnx;
	    psipa2 = ub*[psipar;psipaz]*dnx + uba*[psipr;psipz]*dnx;
	    if abs(psipa2) > 1e-6 % stick to -da/2 flips in the rare event of small psipa2
	      da = -psipa/psipa2;
	    end
	    if abs(da) < 1e-14
	      break
	    end
	  end
	end

	if (psip-psimag)/(psibry-psimag) < 1

	  % perimeter point rp, zp is inside the plasma since normalized flux < 1
	  % so use the NEARX method
	  % Contour enters region rho < dnx at r1, z1
	  % Nearest approach to x is at rn, zn (point along vector [rp-rx,zp-zx])
	  % Contour exits region rho < dnx at r2, z2

	  % Indices		
	  irx = nulls.j(i);
	  izx = nulls.i(i);
	  irnearx = max(3,min(nr-3,irx));
	  iznearx = max(3,min(nz-3,izx));
	  iinearx = iznearx+[-2:3]'*ones(1,6) + ones(6,1)*nz*(irnearx+[-3:2]);

	  % rx, zx, psix
	  rxw = zeros(1,36);
	  zxw = zeros(1,36);
	  rxw(izx-iznearx+(irx-irnearx)*6+ic) = nulls.drdy(i,:);
	  zxw(izx-iznearx+(irx-irnearx)*6+ic) = nulls.dzdy(i,:);
	  psixw = zeros(1,36);
	  psixw(ic) = nulls.w(i,:);

	  % rp, zp, psip, psipr, psipz, psiprho
	  wp = zeros(6,6);
	  wpr = zeros(6,6);
	  wpz = zeros(6,6);
	  wp(izp-iznearx+[2:5],irp-irnearx+[2:5]) = wz0'*wr0;
	  wpr(izp-iznearx+[2:5],irp-irnearx+[2:5]) = wz0'*wr1;
	  wpz(izp-iznearx+[2:5],irp-irnearx+[2:5]) = wz1'*wr0;
	  % d((rp-rx)^2+(zp-zx)^2)/dy = 0	  
	  % d((zx-zp)*psipr + (rp-rx)*psipz)/dy = 0
	  % (rp-rx)*(rpw-rxw) + (zp-zx)*(zpw-zxw) = 0	  
	  % (zx-zp)*(wpr+psiprr*rpw+psiprz*zpw) + (zxw-zpw)*psipr + ...
	  % (rp-rx)*(wpz+psiprz*rpw+psipzz*zpw) + (rpw-rxw)*psipz = 0
	  % rpw =  rxw - (zp-zx)/(rp-rx)*(zpw-zxw)	  
	  % (zx-zp)*wpr + zxw*psipr + (rp-rx)*wpz - rxw*psipz + ...
	  % ((zx-zp)*psiprz+(rp-rx)*psipzz-psipr)*zpw + ...
	  % ((zx-zp)*psiprr+(rp-rx)*psiprz+psipz)*rpw = 0
	  % rpw =  rxw + (zp-zx)/(rp-rx)*zxw - (zp-zx)/(rp-rx)*zpw	  
	  dum1w = (zx-zp)*wpr(:)' + zxw*psipr + (rp-rx)*wpz(:)' - rxw*psipz + ...
	          ((zx-zp)*psiprr+(rp-rx)*psiprz+psipz)*(rxw + (zp-zx)/(rp-rx)*zxw);	  	  
	  dum2 = ((zx-zp)*psiprz+(rp-rx)*psipzz-psipr) + ...
		 ((zx-zp)^2/(rp-rx)*psiprr+(zx-zp)*psiprz+psipz*(zx-zp)/(rp-rx));		 
          zpw = -dum1w/dum2;
	  rpw =  rxw + (zp-zx)/(rp-rx)*zxw - (zp-zx)/(rp-rx)*zpw;
	  psipw = wp(:)' + psipr*rpw + psipz*zpw;
	  psiprw = wpr(:)' + psiprr*rpw + psiprz*zpw;
	  psipzw = wpz(:)' + psiprz*rpw + psipzz*zpw;
	  psiprhow = ((rp-rx)*psiprw+(rpw-rxw)*psipr + (zp-zx)*psipzw+(zpw-zxw)*psipz)/dnx;	  

	  % Find where the boundary ENTERS the nearx region
	  % Rotate vector u1 until r1, z1 is at flux = the boundary flux psibry
	  u1 = up; % initial value of unit vector pointing toward r2, z2
	  da = pi/16; % initial angle increment (positive value makes it go toward entry point)
	  nflip = 0; % number of times the rotation has flipped direction
	  for k = 1:32
	    u1 = [u1(1)*cos(da)-sin(da)*u1(2), sin(da)*u1(1)+cos(da)*u1(2)]; % rotate u1
	    ub = [-u1(2), u1(1)];
	    r1 = rx + dnx*u1(1);
	    z1 = zx + dnx*u1(2);
	    ir1 = floor(r1);
	    iz1 = floor(z1);
	    y44 = psizr(iz1-1:iz1+2,ir1-1:ir1+2);
	    tr = r1-ir1;
	    tz = z1-iz1;
	    wr0 = [1 tr tr^2 tr^3]*mx;
	    wz0 = [1 tz tz^2 tz^3]*mx;
	    wr1 = [0 1 2*tr 3*tr^2]*mx;
	    wz1 = [0 1 2*tz 3*tz^2]*mx;
	    wr2 = [0 0 2 6*tr]*mx;
	    wz2 = [0 0 2 6*tz]*mx;
	    psi1 = wz0*y44*wr0';
	    psi1r = wz0*y44*wr1';
	    psi1z = wz1*y44*wr0';
	    psi1a = ub*[psi1r;psi1z]*dnx;
	    psi1rr = wz0*y44*wr2';
	    psi1rz = wz1*y44*wr1';
	    psi1zz = wz2*y44*wr0';
	    if s*(psi1-psibry)*da < 0 % Heading the wrong way
	      da = -da/2;
	      nflip = nflip+1;
	    end
	    if nflip > 1
	      % Time for Newton-Rhapson to make psi1 = psibry
	      if abs(psi1a) > 1e-6 % stick to -da/2 flips in the rare event of small psi1a
		da = -(psi1-psibry)/psi1a;
	      end
	      if abs(da) < 1e-14
		break
	      end
	    end
	  end
	  w1 = zeros(6,6);
	  w1(iz1-iznearx+[2:5],ir1-irnearx+[2:5]) = wz0'*wr0;
	  w1r = zeros(6,6);
	  w1r(iz1-iznearx+[2:5],ir1-irnearx+[2:5]) = wz0'*wr1;
	  w1z = zeros(6,6);
	  w1z(iz1-iznearx+[2:5],ir1-irnearx+[2:5]) = wz1'*wr0;
	  % Response of r1, z1:
	  % d((r1-rx)^2+(z1-zx)^2)/dy = 0
	  % d(psi1-psibry)/dy = 0	  
	  % (r1-rx)*(r1w-rxw) + (z1-zx)*(z1w-zxw) = 0	  
	  % w1 + psi1r*r1w + psi1z*z1w = 0
	  % psi1r*r1b + psi1z*z1b = 1	  
	  % r1w = rxw -(z1-zx)/(r1-rx)*(z1w-zxw)
	  dumw = w1(:)'+psi1r*(rxw+(z1-zx)/(r1-rx)*zxw);
	  dum = psi1z-psi1r*(z1-zx)/(r1-rx);
	  z1w = -dumw/dum;
	  z1b = 1/dum;
	  r1w = rxw -(z1-zx)/(r1-rx)*(z1w-zxw);
	  r1b = -(z1-zx)/(r1-rx)*z1b;
	  psi1w = w1(:)' + psi1r*r1w + psi1z*z1w;
	  psi1b = psi1r*r1b + psi1z*z1b;
	  psi1rw = w1r(:)' + psi1rr*r1w + psi1rz*z1w;
	  psi1rb = psi1rr*r1b + psi1rz*z1b;
	  psi1zw = w1z(:)' + psi1rz*r1w + psi1zz*z1w;
	  psi1zb = psi1rz*r1b + psi1zz*z1b;


	  % Find where the boundary EXITS the nearx region
	  % Rotate vector u2 until r2, z2 is at flux = the boundary flux psibry	  
	  u2 = up; % initial value of unit vector pointing toward r2, z2
	  da = -pi/16; % initial angle increment (negative value makes it go toward exit point)
	  nflip = 0; % number of times the rotation has flipped direction
	  for k = 1:32
	    u2 = [u2(1)*cos(da)-sin(da)*u2(2), sin(da)*u2(1)+cos(da)*u2(2)]; % rotate u2
	    ub = [-u2(2), u2(1)];
	    r2 = rx + dnx*u2(1);
	    z2 = zx + dnx*u2(2);
	    ir2 = floor(r2);
	    iz2 = floor(z2);
	    tr = r2-ir2;
	    tz = z2-iz2;
	    wr0 = [1 tr tr^2 tr^3]*mx;
	    wz0 = [1 tz tz^2 tz^3]*mx;
	    wr1 = [0 1 2*tr 3*tr^2]*mx;
	    wz1 = [0 1 2*tz 3*tz^2]*mx;
	    wr2 = [0 0 2 6*tr]*mx;
	    wz2 = [0 0 2 6*tz]*mx;
	    y44 = psizr(iz2-1:iz2+2,ir2-1:ir2+2);
	    psi2 = wz0*y44*wr0';
	    psi2r = wz0*y44*wr1';
	    psi2z = wz1*y44*wr0';
	    psi2a = ub*[psi2r;psi2z]*dnx;
	    psi2rr = wz0*y44*wr2';
	    psi2rz = wz1*y44*wr1';
	    psi2zz = wz2*y44*wr0';
	    if s*(psi2-psibry)*da > 0 % Heading the wrong way
	      da = -da/2;
	      nflip = nflip+1;
	    end
	    if nflip > 1
	      % Time for Newton-Rhapson to make psi2 = psibry
	      if abs(psi2a) > 1e-6 % stick to -da/2 flips in the rare event of small psi2a
		da = -(psi2-psibry)/psi2a;
	      end
	      if abs(da) < 1e-14
		break
	      end
	    end
	  end
	  w2 = zeros(6,6);
	  w2(iz2-iznearx+[2:5],ir2-irnearx+[2:5]) = wz0'*wr0;
	  w2r = zeros(6,6);
	  w2r(iz2-iznearx+[2:5],ir2-irnearx+[2:5]) = wz0'*wr1;
	  w2z = zeros(6,6);
	  w2z(iz2-iznearx+[2:5],ir2-irnearx+[2:5]) = wz1'*wr0;
	  dumw = w2(:)'+psi2r*(rxw+(z2-zx)/(r2-rx)*zxw);
	  dum = psi2z-psi2r*(z2-zx)/(r2-rx);
	  z2w = -dumw/dum;
	  z2b = 1/dum;
	  r2w = rxw -(z2-zx)/(r2-rx)*(z2w-zxw);
	  r2b = -(z2-zx)/(r2-rx)*z2b;
	  psi2w = w2(:)' + psi2r*r2w + psi2z*z2w;
	  psi2b = psi2r*r2b + psi2z*z2b;
	  psi2rw = w2r(:)' + psi2rr*r2w + psi2rz*z2w;
	  psi2rb = psi2rr*r2b + psi2rz*z2b;
	  psi2zw = w2z(:)' + psi2rz*r2w + psi2zz*z2w;
	  psi2zb = psi2rz*r2b + psi2zz*z2b;

	  % The modified interpolation scheme only affects point rn, zn
	  % When psizr is changing so that it is approaching the x-point
	  % it will match the position and speed of regular interpolation
	  % at point rp, zp and will then first move slightly faster and then
	  % slower toward the end where speed goes infinite with regular interpolation
	  % Treat contour point nearest x (rn,zn) as if dpsi = cx(1)*rho + cx(2)*rho^3

	  % The gradient at point p
	  psiprho = up*[psipr;psipz];

	  % With modified interpolation:
	  % psip = cx(1)*dnx + cx(2)*dnx^3
	  % psiprho = cx(1) + 3*cx(2)*dnx^2

	  % Make psip and psiprho with modified interpolation match regular interpolation	  
	  cx = [dnx dnx^3; 1 3*dnx^2]\[psip-psix; psiprho];

	  % Response of cx
	  cx1 = cx(1); % cx1 = 3/2*(psip-psix)/dnx-psiprho/2
	  cx1w = 3/2/dnx*(psipw-psixw)-psiprhow/2;
	  cx2 = cx(2); % cx2 = (psiprho*dnx-(psip-psix))/2/dnx^3
	  cx2w = (psiprhow*dnx-(psipw-psixw))/2/dnx^3;

	  % psibry-psix = cx(1)*rhn+cx(2)*rhn^3, find rhn
	  % Check if there is an extreme in the interval: cx(1)+3*cx(2)*rhn^2 = 0
	  rhn2 = -cx(1)/3/cx(2);
	  if rhn2 > 0 & rhn2 < dnx^2
	    rhn = sqrt(rhn2);
	    psin = psix + cx(1)*rhn + cx(2)*rhn^3;
	    if (psix-psibry)*(psin-psibry) <= 0 % psibry is in interval 0:rhn
	      rhn = rhn/2;
	    else % psibry is after rhn but before dnx
	      rhn = (rhn+dnx)/2;
	    end
	  else
	    rhn = dnx/2;
	  end
	  % Initial value of rhn has been created for search
	  for k = 1:7
	    psin = psix + cx(1)*rhn + cx(2)*rhn^3; % Flux at rn, zn
	    psingrad = cx(1) + 3*cx(2)*rhn^2; % Flux gradient
	    drhn = (psibry-psin)/psingrad; % Newton-Rhapson change toward solution
	    rhn = rhn+drhn;
	  end
	  % Response of rhn
	  rhnw = -(psixw + cx1w*rhn + cx2w*rhn^3)/(cx1+3*rhn^2*cx2);
	  rhnb = 1/(cx1 + 3*rhn^2*cx2);

	  % rn, zn
	  rn = rx + rhn*up(1); 
	  zn = zx + rhn*up(2);

	  % Response of rn, zn
	  rnw = (1-rhn/dnx)*rxw + rhn/dnx*rpw + (rp-rx)/dnx*rhnw;
	  rnb = (rp-rx)/dnx*rhnb;
	  znw = (1-rhn/dnx)*zxw + rhn/dnx*zpw + (zp-zx)/dnx*rhnw;
	  znb = (zp-zx)/dnx*rhnb;

	  % Unit vector pointing along boundary at rp, zp
	  ubp = s*[psipz, -psipr]/sqrt(psipr^2+psipz^2);

	  % Response of ubp
	  ubpr = ubp(1);
	  ubprw = s*psipzw/sqrt(psipr^2+psipz^2) - ...
	    ubpr/(psipr^2+psipz^2)*(psipr*psiprw+psipz*psipzw);
	  ubpz = ubp(2);
	  ubpzw = -s*psiprw/sqrt(psipr^2+psipz^2) - ...
	    ubpz/(psipr^2+psipz^2)*(psipr*psiprw+psipz*psipzw);

	  % Vector from 1 to n
	  drn1 = rn-r1;
	  dzn1 = zn-z1;
	  dn1s = drn1^2+dzn1^2;
	  dn1 = sqrt(dn1s);

	  % Response of vector 1 to n
	  drn1w = rnw-r1w;
	  drn1b = rnb-r1b;
	  dzn1w = znw-z1w;
	  dzn1b = znb-z1b;

	  % Unit vector pointing from 1 to n
	  u1n = [drn1, dzn1]/dn1;

	  % Response of u1n
	  dumw = (drn1*drn1w+dzn1*dzn1w)/dn1s;
	  dumb = (drn1*drn1b+dzn1*dzn1b)/dn1s;
	  u1nr = u1n(1);
	  u1nrw = drn1w/dn1 - u1nr*dumw;
	  u1nrb = drn1b/dn1 - u1nr*dumb;
	  u1nz = u1n(2);
	  u1nzw = dzn1w/dn1 - u1nz*dumw;
	  u1nzb = dzn1b/dn1 - u1nz*dumb;

	  % Vector from n to 2
	  dr2n = r2-rn;
	  dz2n = z2-zn;
	  d2ns = dr2n^2+dz2n^2;
	  d2n = sqrt(d2ns);

	  % Response of vector n to 2
	  dr2nw = r2w-rnw;
	  dr2nb = r2b-rnb;
	  dz2nw = z2w-znw;
	  dz2nb = z2b-znb;

	  % Unit vector pointing from n to 2
	  un2 = [dr2n, dz2n]/d2n;

	  % Response of un2
	  dumw = (dr2n*dr2nw+dz2n*dz2nw)/d2ns;
	  dumb = (dr2n*dr2nb+dz2n*dz2nb)/d2ns;
	  un2r = un2(1);
	  un2rw = dr2nw/d2n - un2r*dumw;
	  un2rb = dr2nb/d2n - un2r*dumb;
	  un2z = un2(2);
	  un2zw = dz2nw/d2n - un2z*dumw;
	  un2zb = dz2nb/d2n - un2z*dumb;

	  % Unit vector pointing along boundary at r1, z1
	  d2 = psi1r^2+psi1z^2;
	  d = sqrt(d2);
	  ub1 = s*[psi1z, -psi1r]/d;

	  % Response of ub1
	  dumw = (psi1r*psi1rw+psi1z*psi1zw)/d2;
	  dumb = (psi1r*psi1rb+psi1z*psi1zb)/d2;
	  ub1r = ub1(1);
	  ub1rw = s*psi1zw/d - ub1r*dumw;
	  ub1rb = s*psi1zb/d - ub1r*dumb;
	  ub1z = ub1(2);
	  ub1zw = -s*psi1rw/d - ub1z*dumw;
	  ub1zb = -s*psi1rb/d - ub1z*dumb;

	  % Unit vector pointing along boundary at r2, z2
	  d2 = psi2r^2+psi2z^2;
	  d = sqrt(d2);
	  ub2 = s*[psi2z, -psi2r]/d;

	  % Response of ub2
	  dumw = (psi2r*psi2rw+psi2z*psi2zw)/d2;
	  dumb = (psi2r*psi2rb+psi2z*psi2zb)/d2;
	  ub2r = ub2(1);
	  ub2rw = s*psi2zw/d - ub2r*dumw;
	  ub2rb = s*psi2zb/d - ub2r*dumb;
	  ub2z = ub2(2);
	  ub2zw = -s*psi2rw/d - ub2z*dumw;
	  ub2zb = -s*psi2rb/d - ub2z*dumb;

	  % Unit vector pointing along boundary at rn, zn is gradually broken up
	  % (would otherwise break suddenly when rn,zn becomes x-point)
	  ubn1 = rhn/dnx*ubp + (1-rhn/dnx)*u1n; % vector on the side facing point 1
	  ubn2 = rhn/dnx*ubp + (1-rhn/dnx)*un2; % vector on the side facing point 2

	  % Response of unb1
	  ubn1r = ubn1(1);
	  ubn1rw = rhnw/dnx*ubpr+rhn/dnx*ubprw - rhnw/dnx*u1nr+(1-rhn/dnx)*u1nrw;
	  ubn1rb = rhnb/dnx*ubpr               - rhnb/dnx*u1nr+(1-rhn/dnx)*u1nrb;
	  ubn1z = ubn1(2);
	  ubn1zw = rhnw/dnx*ubpz+rhn/dnx*ubpzw - rhnw/dnx*u1nz+(1-rhn/dnx)*u1nzw;
	  ubn1zb = rhnb/dnx*ubpz               - rhnb/dnx*u1nz+(1-rhn/dnx)*u1nzb;

	  % Response of unb2
	  ubn2r = ubn2(1);
	  ubn2rw = rhnw/dnx*ubpr+rhn/dnx*ubprw - rhnw/dnx*un2r+(1-rhn/dnx)*un2rw;
	  ubn2rb = rhnb/dnx*ubpr               - rhnb/dnx*un2r+(1-rhn/dnx)*un2rb;
	  ubn2z = ubn2(2);
	  ubn2zw = rhnw/dnx*ubpz+rhn/dnx*ubpzw - rhnw/dnx*un2z+(1-rhn/dnx)*un2zw;
	  ubn2zb = rhnb/dnx*ubpz               - rhnb/dnx*un2z+(1-rhn/dnx)*un2zb;

	  % Original boundary points within the nearx region are deleted
	  % New boundary points including rn,zn are added

	  % Flag points to be deleted
	  V = R1 < dnx; % R1 is distances between R0,Z0 and this x-point rx,zx

	  % Number of points to be deleted
	  nremove = sum(double(V));

	  % Remove flagged points and update n
	  N(1:end-nremove) = N(~V);
	  NX(1:end-nremove) = NX(~V);
	  R0(1:end-nremove) = R0(~V);
	  R0b(1:end-nremove) = R0b(~V);
	  R0w(1:end-nremove,:) = R0w(~V,:);
	  Z0(1:end-nremove) = Z0(~V);
	  Z0b(1:end-nremove) = Z0b(~V);
	  Z0w(1:end-nremove,:) = Z0w(~V,:);
	  IR(1:end-nremove) = IR(~V);
	  IZ(1:end-nremove) = IZ(~V);
	  Br0(1:end-nremove) = Br0(~V);
	  Br0b(1:end-nremove) = Br0b(~V);
	  Br0w(1:end-nremove,:) = Br0w(~V,:);
	  Bz0(1:end-nremove) = Bz0(~V);
	  Bz0b(1:end-nremove) = Bz0b(~V);
	  Bz0w(1:end-nremove,:) = Bz0w(~V,:);
	  Br1(1:end-nremove) = Br1(~V);
	  Br1b(1:end-nremove) = Br1b(~V);
	  Br1w(1:end-nremove,:) = Br1w(~V,:);
	  Bz1(1:end-nremove) = Bz1(~V);
	  Bz1b(1:end-nremove) = Bz1b(~V);
	  Bz1w(1:end-nremove,:) = Bz1w(~V,:);

	  n = n-nremove;

	  % If the bdef point was deleted then it must be added back
	  insert_bdef = V(1);

	  % Note that if the bdef point is within the nearx region but isn't the x-point
	  % then it will no longer exist in R0, Z0. In its place will be rn, zn.
	  % This is by design. The bdef point will still exist as rbdef, zbdef and 
	  % it is from those coordinates that psibry should always be calculated

	  % insert will be the index where new points are inserted
	  if insert_bdef
	    insert = n+1;
	  else
	    if nremove == 0
	      if [r2-r1,z2-z1]*[R0(K(1))-rn;Z0(K(1))-zn] < 0
		insert = K(1)+1; % K(1) is before entering nearx region
	      else
		insert = K(1); % K(1) is after entering nearx region
	      end
	    else
	      insert = min(K(1:nremove));
	    end
	  end

	  % Boundary interpolation in the two intervals in the nearx region
	  % is by the same method that will be used for all boundary points
	  % r = r0x(j) + r1x(j)*x + r2x(j)*x^2 + r3x(j)*x^3, with 0<=x<=1
	  % z = z0x(j) + z1x(j)*x + z2x(j)*x^2 + z3x(j)*x^3, with 0<=x<=1

	  drx = [drn1; dr2n];
	  drxb = [drn1b; dr2nb];
	  drxw = [drn1w(:)'; dr2nw(:)'];

	  dzx = [dzn1; dz2n];
	  dzxb = [dzn1b; dz2nb];
	  dzxw = [dzn1w(:)'; dz2nw(:)'];

	  br0 = [ub1r; ubn2r];
	  br0b = [ub1rb; ubn2rb];
	  br0w = [ub1rw(:)'; ubn2rw(:)'];

	  bz0 = [ub1z; ubn2z];
	  bz0b = [ub1zb; ubn2zb];
	  bz0w = [ub1zw(:)'; ubn2zw(:)'];

	  br1 = [ubn1r; ub2r];
	  br1b = [ubn1rb; ub2rb];
	  br1w = [ubn1rw(:)'; ub2rw(:)'];

	  bz1 = [ubn1z; ub2z];
	  bz1b = [ubn1zb; ub2zb];
	  bz1w = [ubn1zw(:)'; ub2zw(:)'];

	  nom = bz0.*drx-br0.*dzx;
	  nomb = bz0b.*drx+bz0.*drxb-br0b.*dzx-br0.*dzxb;
	  nomw = drx*T.*bz0w+bz0*T.*drxw-dzx*T.*br0w-br0*T.*dzxw;

	  den = br0.*drx+bz0.*dzx;
	  denb = br0b.*drx+br0.*drxb+bz0b.*dzx+bz0.*dzxb;
	  denw = drx*T.*br0w+br0*T.*drxw+dzx*T.*bz0w+bz0*T.*dzxw;

	  t0x = nom./den;
	  t0xb = nomb./den-t0x./den.*denb;
	  t0xw = 1./den*T.*nomw-t0x./den*T.*denw;

	  nom = bz1.*drx-br1.*dzx;
	  nomb = bz1b.*drx+bz1.*drxb-br1b.*dzx-br1.*dzxb;
	  nomw = drx*T.*bz1w+bz1*T.*drxw-dzx*T.*br1w-br1*T.*dzxw;

	  den = br1.*drx+bz1.*dzx;
	  denb = br1b.*drx+br1.*drxb+bz1b.*dzx+bz1.*dzxb;
	  denw = drx*T.*br1w+br1*T.*drxw+dzx*T.*bz1w+bz1*T.*dzxw;

	  t1x = nom./den;
	  t1xb = nomb./den-t1x./den.*denb;
	  t1xw = 1./den*T.*nomw-t1x./den*T.*denw;

	  r0x = [r1; rn];
	  r0xb = [r1b; rnb];
	  r0xw = [r1w(:)'; rnw(:)'];

	  r1x = drx-t0x.*dzx;
	  r1xb = drxb-t0xb.*dzx-t0x.*dzxb;
	  r1xw = drxw-dzx*T.*t0xw-t0x*T.*dzxw;

	  r2x = (2*t0x+t1x).*dzx;
	  r2xb = (2*t0xb+t1xb).*dzx + (2*t0x+t1x).*dzxb;
	  r2xw = dzx*T.*(2*t0xw+t1xw) + (2*t0x+t1x)*T.*dzxw;	  

	  r3x = -(t0x+t1x).*dzx;
	  r3xb = -(t0xb+t1xb).*dzx - (t0x+t1x).*dzxb;
	  r3xw = -dzx*T.*(t0xw+t1xw) - (t0x+t1x)*T.*dzxw;

	  z0x = [z1; zn];
	  z0xb = [z1b; znb];
	  z0xw = [z1w(:)'; znw(:)'];

	  z1x = dzx+t0x.*drx;
	  z1xb = dzxb+t0xb.*drx + t0x.*drxb;
	  z1xw = dzxw+drx*T.*t0xw+t0x*T.*drxw;

	  z2x = -(2*t0x+t1x).*drx;
	  z2xb = -(2*t0xb+t1xb).*drx - (2*t0x+t1x).*drxb;
	  z2xw = -drx*T.*(2*t0xw+t1xw) - (2*t0x+t1x)*T.*drxw;

	  z3x = (t0x+t1x).*drx;
	  z3xb = (t0xb+t1xb).*drx + (t0x+t1x).*drxb;
	  z3xw = drx*T.*(t0xw+t1xw) + (t0x+t1x)*T.*drxw;

	  rs = [r0x*xs.^0 + r1x*xs.^1 + r2x*xs.^2 + r3x*xs.^3]';
	  irs = floor(2*rs)/2;
	  zs = [z0x*xs.^0 + z1x*xs.^1 + z2x*xs.^2 + z3x*xs.^3]';
	  izs = floor(2*zs)/2;

	  % Find intersections between grid lines and boundary
	  rhv = [nan nan]; % r of interpolated boundary points
	  zhv = [nan nan]; % z of interpolated boundary points
	  rhvx = [nan nan];% dr/dx
	  zhvx = [nan nan];% dz/dx
	  rhvb = [0 0];
	  zhvb = [0 0];
	  rhvxb = [0 0];
	  zhvxb = [0 0];
	  rhvw = zeros(2,36);
	  zhvw = zeros(2,36);
	  rhvxw = zeros(2,36);
	  zhvxw = zeros(2,36);
	  for j = 1:2 % j = 1 for interval from point 1 to n, j = 2 for interval n to 2
	    for k = 1:nr-1 % Check all xs-intervals for a grid line intersection
	      xhv = [nan nan]; % x where boundary intersects horizontal and vertical line
	      if izs(k,j) ~= izs(k+1,j); % intersected horizontal grid line
		zhv(1) = round(zs(k,j)+zs(k+1,j))/2; % z for grid line (multiple of 0.5)
		x = (xs(k)+xs(k+1))/2; % Initial value for x in interpolation formula
		for m = 1:32 % Iterate up to 32 times
	          z = z0x(j) + z1x(j)*x + z2x(j)*x^2 + z3x(j)*x^3;
	          z_x = z1x(j) + 2*z2x(j)*x + 3*z3x(j)*x^2;
		  dx = max(-0.1,min(0.1,(zhv(1)-z)/z_x));
		  x = max(0,min(1,x+dx));
		  if abs(dx) < 1e-14
		    break % now z = zhv(1) and hence x is what it should be
		  end
		end
		xhv(1) = x;
		rhv(1) = r0x(j) + r1x(j)*x + r2x(j)*x^2 + r3x(j)*x^3;
		zhvx(1) = z1x(j) + 2*z2x(j)*x + 3*z3x(j)*x^2;
		rhvx(1) = r1x(j) + 2*r2x(j)*x + 3*r3x(j)*x^2;
		% zhvw(1) = 0, find xb, xw
		% z0xb(j) + z1xb(j)*x + z2xb(j)*x^2 + z3xb(j)*x^3 + ...
		% z1x(j)*xb + 2*z2x(j)*x*xb + 3*z3x(j)*x^2*xb = 0
		xb = -(z0xb(j) + z1xb(j)*x + z2xb(j)*x^2 + z3xb(j)*x^3)/zhvx(1);
		xw = -(z0xw(j,:) + z1xw(j,:)*x + z2xw(j,:)*x^2 + z3xw(j,:)*x^3)/zhvx(1);
		zhvb(1) = 0;
		rhvb(1) = r0xb(j) + r1xb(j)*x + r2xb(j)*x^2 + r3xb(j)*x^3 + rhvx(1)*xb;
		zhvw(1,:) = 0;
		rhvw(1,:) = r0xw(j,:) + r1xw(j,:)*x + r2xw(j,:)*x^2 + r3xw(j,:)*x^3 + rhvx(1)*xw;
		zhvxb(1) = z1xb(j) + 2*z2xb(j)*x + 3*z3xb(j)*x^2 + (2*z2x(j)+6*z3x(j)*x)*xb;
		rhvxb(1) = r1xb(j) + 2*r2xb(j)*x + 3*r3xb(j)*x^2 + (2*r2x(j)+6*r3x(j)*x)*xb;
		zhvxw(1,:) = z1xw(j,:) + 2*z2xw(j,:)*x + 3*z3xw(j,:)*x^2 + (2*z2x(j)+6*z3x(j)*x)*xw;
		rhvxw(1,:) = r1xw(j,:) + 2*r2xw(j,:)*x + 3*r3xw(j,:)*x^2 + (2*r2x(j)+6*r3x(j)*x)*xw;
	      end
	      if irs(k,j) ~= irs(k+1,j); % intersected vertical grid line
		rhv(2) = round(rs(k,j)+rs(k+1,j))/2; % r for grid line (multiple of 0.5)
		x = (xs(k)+xs(k+1))/2; % Initial value for x in interpolation formula
		for m = 1:32 % Iterate up to 32 times
	          r = r0x(j) + r1x(j)*x + r2x(j)*x^2 + r3x(j)*x^3;
	          r_x = r1x(j) + 2*r2x(j)*x + 3*r3x(j)*x^2;
		  dx = max(-0.1,min(0.1,(rhv(2)-r)/r_x));
		  x = max(0,min(1,x+dx));
		  if abs(dx) < 1e-14
		    break % now r = rhv(2) and hence x is what it should be
		  end
		end
		xhv(2) = x;
		zhv(2) = z0x(j) + z1x(j)*x + z2x(j)*x^2 + z3x(j)*x^3;
		rhvx(2) = r1x(j) + 2*r2x(j)*x + 3*r3x(j)*x^2;
		zhvx(2) = z1x(j) + 2*z2x(j)*x + 3*z3x(j)*x^2;
		% rhvw(2) = 0, find xb, xw
		% r0xb(j) + r1xb(j)*x + r2xb(j)*x^2 + r3xb(j)*x^3 + ...
		% r1x(j)*xb + 2*r2x(j)*x*xb + 3*r3x(j)*x^2*xb = 0
		xb = -(r0xb(j) + r1xb(j)*x + r2xb(j)*x^2 + r3xb(j)*x^3)/rhvx(2);
		xw = -(r0xw(j,:) + r1xw(j,:)*x + r2xw(j,:)*x^2 + r3xw(j,:)*x^3)/rhvx(2);
		rhvb(2) = 0;
		zhvb(2) = z0xb(j) + z1xb(j)*x + z2xb(j)*x^2 + z3xb(j)*x^3 + zhvx(2)*xb;
		rhvw(2,:) = 0;
		zhvw(2,:) = z0xw(j,:) + z1xw(j,:)*x + z2xw(j,:)*x^2 + z3xw(j,:)*x^3 + zhvx(2)*xw;
		rhvxb(2) = r1xb(j) + 2*r2xb(j)*x + 3*r3xb(j)*x^2 + (2*r2x(j)+6*r3x(j)*x)*xb;
		zhvxb(2) = z1xb(j) + 2*z2xb(j)*x + 3*z3xb(j)*x^2 + (2*z2x(j)+6*z3x(j)*x)*xb;
		rhvxw(2,:) = r1xw(j,:) + 2*r2xw(j,:)*x + 3*r3xw(j,:)*x^2 + (2*r2x(j)+6*r3x(j)*x)*xw;
		zhvxw(2,:) = z1xw(j,:) + 2*z2xw(j,:)*x + 3*z3xw(j,:)*x^2 + (2*z2x(j)+6*z3x(j)*x)*xw;
	      end
	      if j == 1 & k == 1 & all(isnan(xhv))
		% Add point r1, z1 since it won't be close to another point
		xhv(1) = 0;
		rhv(1) = r1;
		zhv(1) = z1;
		rhvb(1) = r1b;
		zhvb(1) = z1b;
		rhvw(1,:) = r1w(:)';
		zhvw(1,:) = z1w(:)';
		rhvx(1) = r1x(j);
		zhvx(1) = z1x(j);
		rhvxb(1) = r1xb(j);
		zhvxb(1) = z1xb(j);
		rhvxw(1,:) = r1xw(j,:);
		zhvxw(1,:) = z1xw(j,:);
	      end
	      if j == 2 & k == nr-1 & all(isnan(xhv))
		% Add point r2, z2 since it won't be close to another point
		xhv(1) = 1;
		rhv(1) = r2;
		zhv(1) = z2;
		rhvb(1) = r2b;
		zhvb(1) = z2b;
		rhvw(1,:) = r2w(:)';
		zhvw(1,:) = z2w(:)';
		rhvx(1) = r1x(j)+2*r2x(j)+3*r3x(j);
		zhvx(1) = z1x(j)+2*z2x(j)+3*z3x(j);
		rhvxb(1) = r1xb(j)+2*r2xb(j)+3*r3xb(j);
		zhvxb(1) = z1xb(j)+2*z2xb(j)+3*z3xb(j);
		rhvxw(1,:) = r1xw(j,:)+2*r2xw(j,:)+3*r3xw(j,:);
		zhvxw(1,:) = z1xw(j,:)+2*z2xw(j,:)+3*z3xw(j,:);
	      end
	      [xhv,ihv] = sort(xhv); % Insert points in right order
	      rhv = rhv(ihv);
	      zhv = zhv(ihv);
	      rhvb = rhvb(ihv);
	      zhvb = zhvb(ihv);
	      rhvw = rhvw(ihv,:);
	      zhvw = zhvw(ihv,:);
	      rhvx = rhvx(ihv);
	      zhvx = zhvx(ihv);
	      rhvxb = rhvxb(ihv);
	      zhvxb = zhvxb(ihv);
	      rhvxw = rhvxw(ihv,:);
	      zhvxw = zhvxw(ihv,:);
	      for m = 1:2
		if ~isnan(xhv(m))

		  N(insert+1:end) = N(insert:end-1);	      
		  N(insert) = 0;

		  NX(insert+1:end) = NX(insert:end-1);	      
		  NX(insert) = 1;

		  R0(insert+1:end) = R0(insert:end-1);
		  R0(insert) = rhv(m);
		  R0b(insert+1:end) = R0b(insert:end-1);
		  R0b(insert) = rhvb(m);
		  R0w(insert+1:end,:) = R0w(insert:end-1,:);
		  R0w(insert,:) = rhvw(m,:);

		  Z0(insert+1:end) = Z0(insert:end-1);
		  Z0(insert) = zhv(m);
		  Z0b(insert+1:end) = Z0b(insert:end-1);
		  Z0b(insert) = zhvb(m);
		  Z0w(insert+1:end,:) = Z0w(insert:end-1,:);
		  Z0w(insert,:) = zhvw(m,:);

		  IR(insert+1:end) = IR(insert:end-1);	      
		  IR(insert) = irnearx;

		  IZ(insert+1:end) = IZ(insert:end-1);	      
		  IZ(insert) = iznearx;

		  Br0(insert+1:end) = Br0(insert:end-1);
		  Br0(insert) = rhvx(m);
		  Br0b(insert+1:end) = Br0b(insert:end-1);
		  Br0b(insert) = rhvxb(m);
		  Br0w(insert+1:end,:) = Br0w(insert:end-1,:);
		  Br0w(insert,:) = rhvxw(m,:);

		  Bz0(insert+1:end) = Bz0(insert:end-1);
		  Bz0(insert) = zhvx(m);
		  Bz0b(insert+1:end) = Bz0b(insert:end-1);
		  Bz0b(insert) = zhvxb(m);
		  Bz0w(insert+1:end,:) = Bz0w(insert:end-1,:);
		  Bz0w(insert,:) = zhvxw(m,:);

		  Br1(insert+1:end) = Br1(insert:end-1);
		  Br1(insert) = rhvx(m);
		  Br1b(insert+1:end) = Br1b(insert:end-1);
		  Br1b(insert) = rhvxb(m);
		  Br1w(insert+1:end,:) = Br1w(insert:end-1,:);
		  Br1w(insert,:) = rhvxw(m,:);

		  Bz1(insert+1:end) = Bz1(insert:end-1);
		  Bz1(insert) = zhvx(m);
		  Bz1b(insert+1:end) = Bz1b(insert:end-1);
		  Bz1b(insert) = zhvxb(m);
		  Bz1w(insert+1:end,:) = Bz1w(insert:end-1,:);
		  Bz1w(insert,:) = zhvxw(m,:);

		  insert = insert+1;
		  n = n+1;
		end
	      end	      
	    end
	    if j == 1 % Insert point n

	      N(insert+1:end) = N(insert:end-1);	      
	      N(insert) = 1;

	      NX(insert+1:end) = NX(insert:end-1);	      
	      NX(insert) = 1;

	      R0(insert+1:end) = R0(insert:end-1);
	      R0(insert) = rn;
	      R0b(insert+1:end) = R0b(insert:end-1);
	      R0b(insert) = rnb;
	      R0w(insert+1:end,:) = R0w(insert:end-1,:);
	      R0w(insert,:) = rnw;

	      Z0(insert+1:end) = Z0(insert:end-1);
	      Z0(insert) = zn;
	      Z0b(insert+1:end) = Z0b(insert:end-1);
	      Z0b(insert) = znb;
	      Z0w(insert+1:end,:) = Z0w(insert:end-1,:);
	      Z0w(insert,:) = znw;

	      IR(insert+1:end) = IR(insert:end-1);	      
	      IR(insert) = irnearx;

	      IZ(insert+1:end) = IZ(insert:end-1);	      
	      IZ(insert) = iznearx;

	      Br0(insert+1:end) = Br0(insert:end-1);
	      Br0(insert) = r1x(2);
	      Br0b(insert+1:end) = Br0b(insert:end-1);
	      Br0b(insert) = r1xb(2);
	      Br0w(insert+1:end,:) = Br0w(insert:end-1,:);
	      Br0w(insert,:) = r1xw(2,:);

	      Bz0(insert+1:end) = Bz0(insert:end-1);
	      Bz0(insert) = z1x(2);
	      Bz0b(insert+1:end) = Bz0b(insert:end-1);
	      Bz0b(insert) = z1xb(2);
	      Bz0w(insert+1:end,:) = Bz0w(insert:end-1,:);
	      Bz0w(insert,:) = z1xw(2,:);

	      Br1(insert+1:end) = Br1(insert:end-1);
	      Br1(insert) = r1x(1)+2*r2x(1)+3*r3x(1);
	      Br1b(insert+1:end) = Br1b(insert:end-1);
	      Br1b(insert) = r1xb(1)+2*r2xb(1)+3*r3xb(1);
	      Br1w(insert+1:end,:) = Br1w(insert:end-1,:);
	      Br1w(insert,:) = r1xw(1,:)+2*r2xw(1,:)+3*r3xw(1,:);

	      Bz1(insert+1:end) = Bz1(insert:end-1);
	      Bz1(insert) = z1x(1)+2*z2x(1)+3*z3x(1);
	      Bz1b(insert+1:end) = Bz1b(insert:end-1);
	      Bz1b(insert) = z1xb(1)+2*z2xb(1)+3*z3xb(1);
	      Bz1w(insert+1:end,:) = Bz1w(insert:end-1,:);
	      Bz1w(insert,:) = z1xw(1,:)+2*z2xw(1,:)+3*z3xw(1,:);

	      insert = insert+1;
	      n = n+1;

	      if insert_bdef
		insert = 1;

		N(insert+1:end) = N(insert:end-1);	      
		N(insert) = 1;

		NX(insert+1:end) = NX(insert:end-1);	      
		NX(insert) = 1;

		R0(insert+1:end) = R0(insert:end-1);
		R0(insert) = rn;
		R0b(insert+1:end) = R0b(insert:end-1);
		R0b(insert) = rnb;
		R0w(insert+1:end,:) = R0w(insert:end-1,:);
		R0w(insert,:) = rnw;

		Z0(insert+1:end) = Z0(insert:end-1);
		Z0(insert) = zn;
		Z0b(insert+1:end) = Z0b(insert:end-1);
		Z0b(insert) = znb;
		Z0w(insert+1:end,:) = Z0w(insert:end-1,:);
		Z0w(insert,:) = znw;

		IR(insert+1:end) = IR(insert:end-1);	      
		IR(insert) = irnearx;

		IZ(insert+1:end) = IZ(insert:end-1);	      
		IZ(insert) = iznearx;

		Br0(insert+1:end) = Br0(insert:end-1);
		Br0(insert) = r1x(2);
		Br0b(insert+1:end) = Br0b(insert:end-1);
		Br0b(insert) = r1xb(2);
		Br0w(insert+1:end,:) = Br0w(insert:end-1,:);
		Br0w(insert,:) = r1xw(2,:);

		Bz0(insert+1:end) = Bz0(insert:end-1);
		Bz0(insert) = z1x(2);
		Bz0b(insert+1:end) = Bz0b(insert:end-1);
		Bz0b(insert) = z1xb(2);
		Bz0w(insert+1:end,:) = Bz0w(insert:end-1,:);
		Bz0w(insert,:) = z1xw(2,:);

		Br1(insert+1:end) = Br1(insert:end-1);
		Br1(insert) = r1x(1)+2*r2x(1)+3*r3x(1);
		Br1b(insert+1:end) = Br1b(insert:end-1);
		Br1b(insert) = r1xb(1)+2*r2xb(1)+3*r3xb(1);
		Br1w(insert+1:end,:) = Br1w(insert:end-1,:);
		Br1w(insert,:) = r1xw(1,:)+2*r2xw(1,:)+3*r3xw(1,:);

		Bz1(insert+1:end) = Bz1(insert:end-1);
		Bz1(insert) = z1x(1)+2*z2x(1)+3*z3x(1);
		Bz1b(insert+1:end) = Bz1b(insert:end-1);
		Bz1b(insert) = z1xb(1)+2*z2xb(1)+3*z3xb(1);
		Bz1w(insert+1:end,:) = Bz1w(insert:end-1,:);
		Bz1w(insert,:) = z1xw(1,:)+2*z2xw(1,:)+3*z3xw(1,:);

		insert = insert+1;
		n = n+1;

	      end
	    end % End of for k = 1:nr-1 (scanning xs values)
	  end % End of for j = 1:2 (inserting nearx points for nulls.r(i), nulls.z(i))	  
	end % End of if (psip-psimag)/(psibry-psimag) < 1
      end % End of if R3(1) < 2*dnx
    end % End of if nearxregion(i)
  end % End of for i = 1:nulls.count

  %                       Done inserting NEARX points                         %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  K = IZ+nz*(IR-1);
  I36 = [K-(2*nz+2) K-(2*nz+1) K-2*nz   K-(2*nz-1) K-(2*nz-2) K-(2*nz-3) ...
	 K-(nz+2)   K-(nz+1)   K-nz     K-(nz-1)   K-(nz-2)   K-(nz-3) ...
	 K-2        K-1        K        K+1        K+2        K+3 ...
	 K+(nz-2)   K+(nz-1)   K+nz     K+(nz+1)   K+(nz+2)   K+(nz+3) ...
	 K+(2*nz-2) K+(2*nz-1) K+(2*nz) K+(2*nz+1) K+(2*nz+2) K+(2*nz+3) ...
	 K+(3*nz-2) K+(3*nz-1) K+(3*nz) K+(3*nz+1) K+(3*nz+2) K+(3*nz+3)];  

  % igg is true for grid points covered by plasma
  rranges = nan(nz,6);
  K = max(2,min(nz-2,floor(Z0))); % K = floor(Z0)
  zisint = Z0 == K;
  for i = 1:n
    if zisint(i)
      for j = 1:5
	if isnan(rranges(K(i),j))
	  rranges(K(i),j) = R0(i);
	  break
	end
	if R0(i) < rranges(K(i),j)
	  rranges(K(i),j:6) = [R0(i) rranges(K(i),j:5)];
	  break
	end
	if j == 5
	  rranges(K(i),6) = R0(i);
	end
      end
    end
  end
  for i = 1:nz
    for j = 2:2:6
      if isnan(rranges(i,j))
	break
      end
      j1 = ceil(rranges(i,j-1));
      j2 = floor(rranges(i,j));
      igg(i,j1:j2) = 1;
    end
  end

  V = 2*R0 == floor(2*R0) & ~NX;
  H = 2*Z0 == floor(2*Z0) & ~NX;

  % Derivative of flux Yr, Yz
  J = max(2,min(nr-2,floor(R0))); % J = floor(R0)
  I = max(2,min(nz-2,floor(Z0))); % I = floor(Z0)
  K = I+nz*(J-1);
  II = [K-(nz+1)   K-nz     K-(nz-1)   K-(nz-2) ...
	K-1        K        K+1        K+2      ...
	K+(nz-1)   K+nz     K+(nz+1)   K+(nz+2) ...
	K+(2*nz-1) K+(2*nz) K+(2*nz+1) K+(2*nz+2)];  
  TR = R0-J;
  Wr0 = [l TR TR.^2 TR.^3]*mx;
  Wr1 = [O l 2*TR 3*TR.^2]*mx;
  Wr2 = [O O 2*l 6*TR]*mx;
  TZ = Z0-I;
  Wz0 = [l TZ TZ.^2 TZ.^3]*mx;
  Wz1 = [O l 2*TZ 3*TZ.^2]*mx;
  Wz2 = [O O 2*l 6*TZ]*mx;
  W = ...
  [Wz0(:,1).*Wr0(:,1) Wz0(:,2).*Wr0(:,1) Wz0(:,3).*Wr0(:,1) Wz0(:,4).*Wr0(:,1) ...
   Wz0(:,1).*Wr0(:,2) Wz0(:,2).*Wr0(:,2) Wz0(:,3).*Wr0(:,2) Wz0(:,4).*Wr0(:,2) ...
   Wz0(:,1).*Wr0(:,3) Wz0(:,2).*Wr0(:,3) Wz0(:,3).*Wr0(:,3) Wz0(:,4).*Wr0(:,3) ...
   Wz0(:,1).*Wr0(:,4) Wz0(:,2).*Wr0(:,4) Wz0(:,3).*Wr0(:,4) Wz0(:,4).*Wr0(:,4)];
  Wr = ...
  [Wz0(:,1).*Wr1(:,1) Wz0(:,2).*Wr1(:,1) Wz0(:,3).*Wr1(:,1) Wz0(:,4).*Wr1(:,1) ...
   Wz0(:,1).*Wr1(:,2) Wz0(:,2).*Wr1(:,2) Wz0(:,3).*Wr1(:,2) Wz0(:,4).*Wr1(:,2) ...
   Wz0(:,1).*Wr1(:,3) Wz0(:,2).*Wr1(:,3) Wz0(:,3).*Wr1(:,3) Wz0(:,4).*Wr1(:,3) ...
   Wz0(:,1).*Wr1(:,4) Wz0(:,2).*Wr1(:,4) Wz0(:,3).*Wr1(:,4) Wz0(:,4).*Wr1(:,4)];
  Wz = ...
  [Wz1(:,1).*Wr0(:,1) Wz1(:,2).*Wr0(:,1) Wz1(:,3).*Wr0(:,1) Wz1(:,4).*Wr0(:,1) ...
   Wz1(:,1).*Wr0(:,2) Wz1(:,2).*Wr0(:,2) Wz1(:,3).*Wr0(:,2) Wz1(:,4).*Wr0(:,2) ...
   Wz1(:,1).*Wr0(:,3) Wz1(:,2).*Wr0(:,3) Wz1(:,3).*Wr0(:,3) Wz1(:,4).*Wr0(:,3) ...
   Wz1(:,1).*Wr0(:,4) Wz1(:,2).*Wr0(:,4) Wz1(:,3).*Wr0(:,4) Wz1(:,4).*Wr0(:,4)];
  Wrr = ...
  [Wz0(:,1).*Wr2(:,1) Wz0(:,2).*Wr2(:,1) Wz0(:,3).*Wr2(:,1) Wz0(:,4).*Wr2(:,1) ...
   Wz0(:,1).*Wr2(:,2) Wz0(:,2).*Wr2(:,2) Wz0(:,3).*Wr2(:,2) Wz0(:,4).*Wr2(:,2) ...
   Wz0(:,1).*Wr2(:,3) Wz0(:,2).*Wr2(:,3) Wz0(:,3).*Wr2(:,3) Wz0(:,4).*Wr2(:,3) ...
   Wz0(:,1).*Wr2(:,4) Wz0(:,2).*Wr2(:,4) Wz0(:,3).*Wr2(:,4) Wz0(:,4).*Wr2(:,4)];
  Wrz = ...
  [Wz1(:,1).*Wr1(:,1) Wz1(:,2).*Wr1(:,1) Wz1(:,3).*Wr1(:,1) Wz1(:,4).*Wr1(:,1) ...
   Wz1(:,1).*Wr1(:,2) Wz1(:,2).*Wr1(:,2) Wz1(:,3).*Wr1(:,2) Wz1(:,4).*Wr1(:,2) ...
   Wz1(:,1).*Wr1(:,3) Wz1(:,2).*Wr1(:,3) Wz1(:,3).*Wr1(:,3) Wz1(:,4).*Wr1(:,3) ...
   Wz1(:,1).*Wr1(:,4) Wz1(:,2).*Wr1(:,4) Wz1(:,3).*Wr1(:,4) Wz1(:,4).*Wr1(:,4)];
  Wzz = ...
  [Wz2(:,1).*Wr0(:,1) Wz2(:,2).*Wr0(:,1) Wz2(:,3).*Wr0(:,1) Wz2(:,4).*Wr0(:,1) ...
   Wz2(:,1).*Wr0(:,2) Wz2(:,2).*Wr0(:,2) Wz2(:,3).*Wr0(:,2) Wz2(:,4).*Wr0(:,2) ...
   Wz2(:,1).*Wr0(:,3) Wz2(:,2).*Wr0(:,3) Wz2(:,3).*Wr0(:,3) Wz2(:,4).*Wr0(:,3) ...
   Wz2(:,1).*Wr0(:,4) Wz2(:,2).*Wr0(:,4) Wz2(:,3).*Wr0(:,4) Wz2(:,4).*Wr0(:,4)];
  Yr(1:n) = sum(Wr(1:n,:)'.*psizr(II(1:n,:))');
  Yz(1:n) = sum(Wz(1:n,:)'.*psizr(II(1:n,:))');
  Yrr(1:n) = sum(Wrr(1:n,:)'.*psizr(II(1:n,:))');
  Yrz(1:n) = sum(Wrz(1:n,:)'.*psizr(II(1:n,:))');
  Yzz(1:n) = sum(Wzz(1:n,:)'.*psizr(II(1:n,:))');

  % Response for R0, Z0
  R0w(H,ic) = -s./Yr(H)*S.*W(H,:);
  R0b(H) = s./Yr(H);
  Z0w(V,ic) = -s./Yz(V)*S.*W(V,:);
  Z0b(V) = s./Yz(V);
  if ~NX(1)
    R0w(1,ic) = bd.drdpsi(1,:);
    R0w(n,ic) = bd.drdpsi(1,:);
    R0b([1 n]) = 0;
    Z0w(1,ic) = bd.dzdpsi(1,:);
    Z0w(n,ic) = bd.dzdpsi(1,:);
    Z0b([1 n]) = 0;
  end

  % Response of Yr, Yz at stationary R0, Z0
  Yrw = zeros(nbmax,36);
  Yzw = zeros(nbmax,36);
  Yrw(:,ic) = Wr;
  Yzw(:,ic) = Wz;

  % I, J can be shifted by 1 compared to IZ, IR
  for i = 1:n
    if ~NX(i)
      if J(i) < IR(i)
	R0w(i,:) = [R0w(i,7:36) zeros(1,6)];
	Z0w(i,:) = [Z0w(i,7:36) zeros(1,6)];
      end
      if J(i) > IR(i)
	R0w(i,:) = [zeros(1,6) R0w(i,1:30)];
	Z0w(i,:) = [zeros(1,6) Z0w(i,1:30)];
      end
      if I(i) < IZ(i)
	R0w(i,:) = [R0w(i,2:6) 0 R0w(i,8:12) 0 R0w(i,14:18) 0 R0w(i,20:24) 0 R0w(i,26:30) 0 R0w(i,32:36) 0];
	Z0w(i,:) = [Z0w(i,2:6) 0 Z0w(i,8:12) 0 Z0w(i,14:18) 0 Z0w(i,20:24) 0 Z0w(i,26:30) 0 Z0w(i,32:36) 0];
      end
      if I(i) > IZ(i)
	R0w(i,:) = [0 R0w(i,1:5) 0 R0w(i,7:11) 0 R0w(i,13:17) 0 R0w(i,19:23) 0 R0w(i,25:29) 0 R0w(i,31:35)];
	Z0w(i,:) = [0 Z0w(i,1:5) 0 Z0w(i,7:11) 0 Z0w(i,13:17) 0 Z0w(i,19:23) 0 Z0w(i,25:29) 0 Z0w(i,31:35)];
      end
    end
    if J(i) < IR(i)
      Yrw(i,:) = [Yrw(i,7:36) zeros(1,6)];
      Yzw(i,:) = [Yzw(i,7:36) zeros(1,6)];
    end
    if J(i) > IR(i)
      Yrw(i,:) = [zeros(1,6) Yrw(i,1:30)];
      Yzw(i,:) = [zeros(1,6) Yzw(i,1:30)];
    end
    if I(i) < IZ(i)
      Yrw(i,:) = [Yrw(i,2:6) 0 Yrw(i,8:12) 0 Yrw(i,14:18) 0 Yrw(i,20:24) 0 Yrw(i,26:30) 0 Yrw(i,32:36) 0];
      Yzw(i,:) = [Yzw(i,2:6) 0 Yzw(i,8:12) 0 Yzw(i,14:18) 0 Yzw(i,20:24) 0 Yzw(i,26:30) 0 Yzw(i,32:36) 0];
    end
    if I(i) > IZ(i)
      Yrw(i,:) = [0 Yrw(i,1:5) 0 Yrw(i,7:11) 0 Yrw(i,13:17) 0 Yrw(i,19:23) 0 Yrw(i,25:29) 0 Yrw(i,31:35)];
      Yzw(i,:) = [0 Yzw(i,1:5) 0 Yzw(i,7:11) 0 Yzw(i,13:17) 0 Yzw(i,19:23) 0 Yzw(i,25:29) 0 Yzw(i,31:35)];
    end
  end

  % Response of Yr, Yz including displacement of R0, Z0
  Yrw = Yrw + Yrr*T.*R0w + Yrz*T.*Z0w;
  Yzw = Yzw + Yrz*T.*R0w + Yzz*T.*Z0w;
  Yrb = Yrr.*R0b + Yrz.*Z0b;
  Yzb = Yrz.*R0b + Yzz.*Z0b;

  % Square of gradient
  Yg2 = Yr.^2+Yz.^2;

  % Br, Bz are vectors that point along boundary
  Br0(~NX) = -Yz(~NX);
  Bz0(~NX) = +Yr(~NX);
  Br1(~NX) = -Yz(~NX);
  Bz1(~NX) = +Yr(~NX);
  Br0b(~NX) = -Yzb(~NX);
  Bz0b(~NX) = +Yrb(~NX);
  Br1b(~NX) = -Yzb(~NX);
  Bz1b(~NX) = +Yrb(~NX);
  Br0w(~NX,:) = -Yzw(~NX,:);
  Bz0w(~NX,:) = +Yrw(~NX,:);
  Br1w(~NX,:) = -Yzw(~NX,:);
  Bz1w(~NX,:) = +Yrw(~NX,:);
  
  % Interpolation between boundary points with 3:rd degree polynomials
  % decreases errors by 1000 times compared to linear interpolation
  % Interpolated values between j and j+1 are:
  % r(j) + r1(j)*x + r2(j)*x^2 + r3(j)*x^3
  % z(j) + z1(j)*x + z2(j)*x^2 + z3(j)*x^3
  % where x is a parameter that goes from 0 to 1
  % The important reason for using this interpolation method is that it
  % eliminates discontinuities in the rate at which surface area changes.
  % These otherwise happen when the boundary moves across grid lines

  % R02, Z02, Br1, Bz1
  R02 = [R0(2:end); nan];
  Z02 = [Z0(2:end); nan];
  R02b = [R0b(2:end); 0];
  Z02b = [Z0b(2:end); 0];
  R02w = [R0w(2:end,:); zeros(1,36)];
  Z02w = [Z0w(2:end,:); zeros(1,36)];
  Br1 = [Br1(2:end); nan];
  Bz1 = [Bz1(2:end); nan];
  Br1b = [Br1b(2:end); 0];
  Bz1b = [Bz1b(2:end); 0];
  Br1w = [Br1w(2:end,:); zeros(1,36)];
  Bz1w = [Bz1w(2:end,:); zeros(1,36)];
  % shift the weights so they align with IR, IZ
  for i = 1:n-1
    if IR(i+1) < IR(i)
      R02w(i,:) = [R02w(i,7:36) zeros(1,6)];
      Z02w(i,:) = [Z02w(i,7:36) zeros(1,6)];
      Br1w(i,:) = [Br1w(i,7:36) zeros(1,6)];
      Bz1w(i,:) = [Bz1w(i,7:36) zeros(1,6)];
    end
    if IR(i+1) > IR(i)
      R02w(i,:) = [zeros(1,6) R02w(i,1:30)];
      Z02w(i,:) = [zeros(1,6) Z02w(i,1:30)];
      Br1w(i,:) = [zeros(1,6) Br1w(i,1:30)];
      Bz1w(i,:) = [zeros(1,6) Bz1w(i,1:30)];
    end
    if IZ(i+1) < IZ(i)
      R02w(i,:) = [R02w(i,2:6) 0 R02w(i,8:12) 0 R02w(i,14:18) 0 R02w(i,20:24) 0 R02w(i,26:30) 0 R02w(i,32:36) 0];
      Z02w(i,:) = [Z02w(i,2:6) 0 Z02w(i,8:12) 0 Z02w(i,14:18) 0 Z02w(i,20:24) 0 Z02w(i,26:30) 0 Z02w(i,32:36) 0];
      Br1w(i,:) = [Br1w(i,2:6) 0 Br1w(i,8:12) 0 Br1w(i,14:18) 0 Br1w(i,20:24) 0 Br1w(i,26:30) 0 Br1w(i,32:36) 0];
      Bz1w(i,:) = [Bz1w(i,2:6) 0 Bz1w(i,8:12) 0 Bz1w(i,14:18) 0 Bz1w(i,20:24) 0 Bz1w(i,26:30) 0 Bz1w(i,32:36) 0];
    end
    if IZ(i+1) > IZ(i)
      R02w(i,:) = [0 R02w(i,1:5) 0 R02w(i,7:11) 0 R02w(i,13:17) 0 R02w(i,19:23) 0 R02w(i,25:29) 0 R02w(i,31:35)];
      Z02w(i,:) = [0 Z02w(i,1:5) 0 Z02w(i,7:11) 0 Z02w(i,13:17) 0 Z02w(i,19:23) 0 Z02w(i,25:29) 0 Z02w(i,31:35)];
      Br1w(i,:) = [0 Br1w(i,1:5) 0 Br1w(i,7:11) 0 Br1w(i,13:17) 0 Br1w(i,19:23) 0 Br1w(i,25:29) 0 Br1w(i,31:35)];
      Bz1w(i,:) = [0 Bz1w(i,1:5) 0 Bz1w(i,7:11) 0 Bz1w(i,13:17) 0 Bz1w(i,19:23) 0 Bz1w(i,25:29) 0 Bz1w(i,31:35)];
    end
  end

  DR = R02-R0;
  DZ = Z02-Z0;
  DRb = R02b-R0b;
  DZb = Z02b-Z0b;
  DRw = R02w-R0w;
  DZw = Z02w-Z0w;

  NOM = Bz0.*DR - Br0.*DZ;
  NOMb = Bz0b.*DR+Bz0.*DRb - Br0b.*DZ-Br0.*DZb;
  NOMw = DR*T.*Bz0w+Bz0*T.*DRw - DZ*T.*Br0w-Br0*T.*DZw;

  DEN = Bz0.*DZ + Br0.*DR;
  DENb = Bz0b.*DZ+Bz0.*DZb + Br0b.*DR+Br0.*DRb;
  DENw = DZ*T.*Bz0w+Bz0*T.*DZw + DR*T.*Br0w+Br0*T.*DRw;

  T0 = NOM./DEN;
  T0b = NOMb./DEN - T0./DEN.*DENb;
  T0w = 1./DEN*T.*NOMw - T0./DEN*T.*DENw;

  NOM = Bz1.*DR - Br1.*DZ;
  NOMb = Bz1b.*DR+Bz1.*DRb - Br1b.*DZ-Br1.*DZb;
  NOMw = DR*T.*Bz1w+Bz1*T.*DRw - DZ*T.*Br1w-Br1*T.*DZw;

  DEN = Bz1.*DZ + Br1.*DR;
  DENb = Bz1b.*DZ+Bz1.*DZb + Br1b.*DR+Br1.*DRb;
  DENw = DZ*T.*Bz1w+Bz1*T.*DZw + DR*T.*Br1w+Br1*T.*DRw;

  T1 = NOM./DEN;
  T1b = NOMb./DEN - T1./DEN.*DENb;
  T1w = 1./DEN*T.*NOMw - T1./DEN*T.*DENw;

  if lim == 0 & dnx <= 0
    ybrr = wbrr*psizr(iib);
    ybrz = wbrz*psizr(iib);
    ybzz = wbzz*psizr(iib);
    ybrrr = wbrrr*psizr(iib);
    ybrrz = wbrrz*psizr(iib);
    ybrzz = wbrzz*psizr(iib);
    ybzzz = wbzzz*psizr(iib);
    ybrrw = ybrrr*R0w(1,:)+ybrrz*Z0w(1,:);
    ybrzw = ybrrz*R0w(1,:)+ybrzz*Z0w(1,:);
    ybzzw = ybrzz*R0w(1,:)+ybzzz*Z0w(1,:);
    ybrrw(ic) = ybrrw(ic) + wbrr;
    ybrzw(ic) = ybrzw(ic) + wbrz;
    ybzzw(ic) = ybzzw(ic) + wbzz;

    t = -ybrz/ybzz; 
    tw = -ybrzw/ybzz-t/ybzz*ybzzw;
    s2 = t^2-ybrr/ybzz;
    sw = 2*t*tw - ybrrw/ybzz + ybrr/ybzz^2*ybzzw;
    if s2 > 0
      s1 = sqrt(s2);
      sw = 1/2/s1*sw;
    else
      t = 1;
      s1 = 0; % Should never happen
      sw = 0*T;
    end
    u1 = [1; t-s1]/sqrt(1+(t-s1)^2);
    u2 = [1; t+s1]/sqrt(1+(t+s1)^2);
    if abs([DR(1) DZ(1)]*u1) > abs([DR(1) DZ(1)]*u2)
      dzdr = t+s1*[-1; 1];
      dzdrw = [tw-sw; tw+sw];
    else
      dzdr = t+s1*[1; -1];
      dzdrw = [tw+sw; tw-sw];
    end

    T0(1) = (dzdr(1)*DR(1)-DZ(1))/(dzdr(1)*DZ(1)+DR(1));

    T0b(1) = (dzdr(1)*DRb(1)-DZb(1))/(dzdr(1)*DZ(1)+DR(1)) - ...
      T0(1)/(dzdr(1)*DZ(1)+DR(1))*(dzdr(1)*DZb(1)+DRb(1));

    T0w(1,:) = (dzdrw(1,:)*DR(1)+dzdr(1)*DRw(1,:)-DZw(1,:))/(dzdr(1)*DZ(1)+DR(1))-...
      T0(1)/(dzdr(1)*DZ(1)+DR(1))*(dzdrw(1,:)*DZ(1)+dzdr(1)*DZw(1,:)+DRw(1,:));

    T1(n-1) = (dzdr(2)*DR(n-1)-DZ(n-1))/(dzdr(2)*DZ(n-1)+DR(n-1));

    T1b(n-1) = (dzdr(2)*DRb(n-1)-DZb(n-1))/(dzdr(2)*DZ(n-1)+DR(n-1)) - ...
      T1(n-1)/(dzdr(2)*DZ(n-1)+DR(n-1))*(dzdr(2)*DZb(n-1)+DRb(n-1));

    T1w(n-1,:) = (dzdrw(2,:)*DR(n-1)+dzdr(2)*DRw(n-1,:)-DZw(n-1,:))/(dzdr(2)*DZ(n-1)+DR(n-1))-...
      T1(n-1)/(dzdr(2)*DZ(n-1)+DR(n-1))*(dzdrw(2,:)*DZ(n-1)+dzdr(2)*DZw(n-1,:)+DRw(n-1,:));
  end

  % Coefficients for interpolation [floating index]
  R1 = DR-T0.*DZ;
  R2 = (2*T0+T1).*DZ;
  R3 = -(T0+T1).*DZ;
  Z1 = DZ+T0.*DR;
  Z2 = -(2*T0+T1).*DR;
  Z3 = (T0+T1).*DR;

  R1b = DRb-T0b.*DZ-T0.*DZb;
  R2b = (2*T0b+T1b).*DZ+(2*T0+T1).*DZb;
  R3b = -(T0b+T1b).*DZ-(T0+T1).*DZb;
  Z1b = DZb+T0b.*DR+T0.*DRb;
  Z2b = -(2*T0b+T1b).*DR-(2*T0+T1).*DRb;
  Z3b = (T0b+T1b).*DR+(T0+T1).*DRb;

  R1w = DRw-DZ*T.*T0w-T0*T.*DZw;
  R2w = DZ*T.*(2*T0w+T1w)+(2*T0+T1)*T.*DZw;
  R3w = -DZ*T.*(T0w+T1w)-(T0+T1)*T.*DZw;
  Z1w = DZw+DR*T.*T0w+T0*T.*DRw;
  Z2w = -DR*T.*(2*T0w+T1w)-(2*T0+T1)*T.*DRw;
  Z3w = DR*T.*(T0w+T1w)+(T0+T1)*T.*DRw;

  % IC, JC index cells cut by boundary in interpolated intervals
  IC = round(Z0+0.5*Z1+0.25*Z2+0.125*Z3);
  JC = round(R0+0.5*R1+0.25*R2+0.125*R3);
  
  % Flags V, H will now include nearx points
  V = 2*R0 == floor(2*R0);
  H = 2*Z0 == floor(2*Z0);
  
  % Remember R0, Z0 in index units
  IR0 = R0;
  IZ0 = Z0;

  % Change units from floating index to physics
  rmaxis = (rmaxis-1)*dr+rg(1);
  zmaxis = (zmaxis-1)*dz+zg(1);
  rbdef = (rbdef-1)*dr+rg(1);
  zbdef = (zbdef-1)*dz+zg(1);
  DR = DR*dr;
  R0 = (R0-1)*dr+rg(1);
  R1 = R1*dr;
  R2 = R2*dr;
  R3 = R3*dr;
  DZ = DZ*dz;
  Z0 = (Z0-1)*dz+zg(1);
  Z1 = Z1*dz;
  Z2 = Z2*dz;
  Z3 = Z3*dz;
  DRb = DRb*dr;
  R0b = R0b*dr;
  R1b = R1b*dr;
  R2b = R2b*dr;
  R3b = R3b*dr;
  DZb = DZb*dz;
  Z0b = Z0b*dz;
  Z1b = Z1b*dz;
  Z2b = Z2b*dz;
  Z3b = Z3b*dz;
  DRw = DRw*dr;
  R0w = R0w*dr;
  R1w = R1w*dr;
  R2w = R2w*dr;
  R3w = R3w*dr;
  DZw = DZw*dz;
  Z0w = Z0w*dz;
  Z1w = Z1w*dz;
  Z2w = Z2w*dz;
  Z3w = Z3w*dz;

  % dC = distances *along boundary* = integral(dR^2+dZ^2)
  % dRdx = R1+2*R2*x+3*R3*x^2, dZdx = Z1+2*Z2*x+3*Z3*x^2
  % dC = Integral[sqrt(dRdx^2+dZdx^2)*dx] from x = 0 to x = 1
  % Must be solved numerically by chopping into smaller pieces
  Rs = R1*xs+R2*xs.^2+R3*xs.^3;
  Zs = Z1*xs+Z2*xs.^2+Z3*xs.^3;
  dRs = diff(Rs');
  dZs = diff(Zs');
  dC = sum(sqrt(dRs.^2+dZs.^2))';
  dCb = zeros(nbmax,1);
  dCw = zeros(nbmax,36);
  for i = 1:n-1
    R((i-1)*(nr-1)+[1:nr]) = R0(i) + Rs(i,:);
    Z((i-1)*(nr-1)+[1:nr]) = Z0(i) + Zs(i,:);
    dRw = diff((xs'*R1w(i,:)+xs'.^2*R2w(i,:)+xs'.^3*R3w(i,:)));
    dZw = diff((xs'*Z1w(i,:)+xs'.^2*Z2w(i,:)+xs'.^3*Z3w(i,:)));
    dCw(i,:) = sum((dRs(:,i)*T.*dRw+dZs(:,i)*T.*dZw)./(sqrt(dRs(:,i).^2+dZs(:,i).^2)*T));
    dRb = diff((xs'*R1b(i)+xs'.^2*R2b(i)+xs'.^3*R3b(i)));
    dZb = diff((xs'*Z1b(i)+xs'.^2*Z2b(i)+xs'.^3*Z3b(i)));
    dCb(i) = sum((dRs(:,i).*dRb+dZs(:,i).*dZb)./(sqrt(dRs(:,i).^2+dZs(:,i).^2)));
  end

  % dA0 = integral(r(x)*dzdx) from x = 0 to 1
  dA0 =   R0 .* DZ + ...
      1/2*R1 .* Z1 + ...
      2/3*R1 .* Z2 + ...
      3/4*R1 .* Z3 + ...
      1/3*R2 .* Z1 + ...
      2/4*R2 .* Z2 + ...
      3/5*R2 .* Z3 + ...
      1/4*R3 .* Z1 + ...
      2/5*R3 .* Z2 + ...
      3/6*R3 .* Z3;
  dA0b = (R0 .* DZb+DZ .* R0b) + ...
     1/2*(R1 .* Z1b+Z1 .* R1b) + ...
     2/3*(R1 .* Z2b+Z2 .* R1b) + ...
     3/4*(R1 .* Z3b+Z3 .* R1b) + ...
     1/3*(R2 .* Z1b+Z1 .* R2b) + ...
     2/4*(R2 .* Z2b+Z2 .* R2b) + ...
     3/5*(R2 .* Z3b+Z3 .* R2b) + ...
     1/4*(R3 .* Z1b+Z1 .* R3b) + ...
     2/5*(R3 .* Z2b+Z2 .* R3b) + ...
     3/6*(R3 .* Z3b+Z3 .* R3b);
  dA0w = (R0*T.*DZw+DZ*T.*R0w) + ...
     1/2*(R1*T.*Z1w+Z1*T.*R1w) + ...
     2/3*(R1*T.*Z2w+Z2*T.*R1w) + ...
     3/4*(R1*T.*Z3w+Z3*T.*R1w) + ...
     1/3*(R2*T.*Z1w+Z1*T.*R2w) + ...
     2/4*(R2*T.*Z2w+Z2*T.*R2w) + ...
     3/5*(R2*T.*Z3w+Z3*T.*R2w) + ...
     1/4*(R3*T.*Z1w+Z1*T.*R3w) + ...
     2/5*(R3*T.*Z2w+Z2*T.*R3w) + ...
     3/6*(R3*T.*Z3w+Z3*T.*R3w);

  % dA1 = integral(x*r(x)*dzdx) from x = 0 to 1
  dA1 =   1/2*R0 .* Z1 + ...
          2/3*R0 .* Z2 + ...
	  3/4*R0 .* Z3 + ...
	  1/3*R1 .* Z1 + ...
	  2/4*R1 .* Z2 + ...
	  3/5*R1 .* Z3 + ...
          1/4*R2 .* Z1 + ...
	  2/5*R2 .* Z2 + ...
	  3/6*R2 .* Z3 + ...
	  1/5*R3 .* Z1 + ...
	  2/6*R3 .* Z2 + ...
	  3/7*R3 .* Z3;
  dA1b = 1/2*(R0 .* Z1b+Z1 .* R0b) + ...
	 2/3*(R0 .* Z2b+Z2 .* R0b) + ...
	 3/4*(R0 .* Z3b+Z3 .* R0b) + ...
	 1/3*(R1 .* Z1b+Z1 .* R1b) + ...
	 2/4*(R1 .* Z2b+Z2 .* R1b) + ...
	 3/5*(R1 .* Z3b+Z3 .* R1b) + ...
	 1/4*(R2 .* Z1b+Z1 .* R2b) + ...
	 2/5*(R2 .* Z2b+Z2 .* R2b) + ...
	 3/6*(R2 .* Z3b+Z3 .* R2b) + ...
	 1/5*(R3 .* Z1b+Z1 .* R3b) + ...
	 2/6*(R3 .* Z2b+Z2 .* R3b) + ...
	 3/7*(R3 .* Z3b+Z3 .* R3b);
  dA1w = 1/2*(R0*T.*Z1w+Z1*T.*R0w) + ...
	 2/3*(R0*T.*Z2w+Z2*T.*R0w) + ...
	 3/4*(R0*T.*Z3w+Z3*T.*R0w) + ...
	 1/3*(R1*T.*Z1w+Z1*T.*R1w) + ...
	 2/4*(R1*T.*Z2w+Z2*T.*R1w) + ...
	 3/5*(R1*T.*Z3w+Z3*T.*R1w) + ...
	 1/4*(R2*T.*Z1w+Z1*T.*R2w) + ...
	 2/5*(R2*T.*Z2w+Z2*T.*R2w) + ...
	 3/6*(R2*T.*Z3w+Z3*T.*R2w) + ...
	 1/5*(R3*T.*Z1w+Z1*T.*R3w) + ...
	 2/6*(R3*T.*Z2w+Z2*T.*R3w) + ...
	 3/7*(R3*T.*Z3w+Z3*T.*R3w);

  % dA2 = integral(x^2*r(x)*dzdx) from x = 0 to 1
  dA2 =   1/3*R0 .* Z1 + ...
          2/4*R0 .* Z2 + ...
	  3/5*R0 .* Z3 + ...
	  1/4*R1 .* Z1 + ...
	  2/5*R1 .* Z2 + ...
	  3/6*R1 .* Z3 + ...
          1/5*R2 .* Z1 + ...
	  2/6*R2 .* Z2 + ...
	  3/7*R2 .* Z3 + ...
	  1/6*R3 .* Z1 + ...
	  2/7*R3 .* Z2 + ...
	  3/8*R3 .* Z3;
  dA2b = 1/3*(R0 .* Z1b+Z1 .* R0b) + ...
	 2/4*(R0 .* Z2b+Z2 .* R0b) + ...
	 3/5*(R0 .* Z3b+Z3 .* R0b) + ...
	 1/4*(R1 .* Z1b+Z1 .* R1b) + ...
	 2/5*(R1 .* Z2b+Z2 .* R1b) + ...
	 3/6*(R1 .* Z3b+Z3 .* R1b) + ...
	 1/5*(R2 .* Z1b+Z1 .* R2b) + ...
	 2/6*(R2 .* Z2b+Z2 .* R2b) + ...
	 3/7*(R2 .* Z3b+Z3 .* R2b) + ...
	 1/6*(R3 .* Z1b+Z1 .* R3b) + ...
	 2/7*(R3 .* Z2b+Z2 .* R3b) + ...
	 3/8*(R3 .* Z3b+Z3 .* R3b);
  dA2w = 1/3*(R0*T.*Z1w+Z1*T.*R0w) + ...
	 2/4*(R0*T.*Z2w+Z2*T.*R0w) + ...
	 3/5*(R0*T.*Z3w+Z3*T.*R0w) + ...
	 1/4*(R1*T.*Z1w+Z1*T.*R1w) + ...
	 2/5*(R1*T.*Z2w+Z2*T.*R1w) + ...
	 3/6*(R1*T.*Z3w+Z3*T.*R1w) + ...
	 1/5*(R2*T.*Z1w+Z1*T.*R2w) + ...
	 2/6*(R2*T.*Z2w+Z2*T.*R2w) + ...
	 3/7*(R2*T.*Z3w+Z3*T.*R2w) + ...
	 1/6*(R3*T.*Z1w+Z1*T.*R3w) + ...
	 2/7*(R3*T.*Z2w+Z2*T.*R3w) + ...
	 3/8*(R3*T.*Z3w+Z3*T.*R3w);

  % dA3 = integral(x^3*r(x)*dzdx) from x = 0 to 1
  dA3 =   1/4*R0 .* Z1 + ...
          2/5*R0 .* Z2 + ...
	  3/6*R0 .* Z3 + ...
	  1/5*R1 .* Z1 + ...
	  2/6*R1 .* Z2 + ...
	  3/7*R1 .* Z3 + ...
          1/6*R2 .* Z1 + ...
	  2/7*R2 .* Z2 + ...
	  3/8*R2 .* Z3 + ...
	  1/7*R3 .* Z1 + ...
	  2/8*R3 .* Z2 + ...
	  3/9*R3 .* Z3;
  dA3b = 1/4*(R0 .* Z1b+Z1 .* R0b) + ...
	 2/5*(R0 .* Z2b+Z2 .* R0b) + ...
	 3/6*(R0 .* Z3b+Z3 .* R0b) + ...
	 1/5*(R1 .* Z1b+Z1 .* R1b) + ...
	 2/6*(R1 .* Z2b+Z2 .* R1b) + ...
	 3/7*(R1 .* Z3b+Z3 .* R1b) + ...
	 1/6*(R2 .* Z1b+Z1 .* R2b) + ...
	 2/7*(R2 .* Z2b+Z2 .* R2b) + ...
	 3/8*(R2 .* Z3b+Z3 .* R2b) + ...
	 1/7*(R3 .* Z1b+Z1 .* R3b) + ...
	 2/8*(R3 .* Z2b+Z2 .* R3b) + ...
	 3/9*(R3 .* Z3b+Z3 .* R3b);
  dA3w = 1/4*(R0*T.*Z1w+Z1*T.*R0w) + ...
	 2/5*(R0*T.*Z2w+Z2*T.*R0w) + ...
	 3/6*(R0*T.*Z3w+Z3*T.*R0w) + ...
	 1/5*(R1*T.*Z1w+Z1*T.*R1w) + ...
	 2/6*(R1*T.*Z2w+Z2*T.*R1w) + ...
	 3/7*(R1*T.*Z3w+Z3*T.*R1w) + ...
	 1/6*(R2*T.*Z1w+Z1*T.*R2w) + ...
	 2/7*(R2*T.*Z2w+Z2*T.*R2w) + ...
	 3/8*(R2*T.*Z3w+Z3*T.*R2w) + ...
	 1/7*(R3*T.*Z1w+Z1*T.*R3w) + ...
	 2/8*(R3*T.*Z2w+Z2*T.*R3w) + ...
	 3/9*(R3*T.*Z3w+Z3*T.*R3w);


  % dRA = integral(r(x)^2/2*dzdx) from x = 0 to 1
  dRA = (R0.*dA0 + R1.*dA1 + R2.*dA2 + R3.*dA3)/2;
  dRAb = 1/2*((R0 .* dA0b+dA0 .* R0b) + ...
              (R1 .* dA1b+dA1 .* R1b) + ...
              (R2 .* dA2b+dA2 .* R2b) + ...
              (R3 .* dA3b+dA3 .* R3b));
  dRAw = 1/2*((R0*T.*dA0w+dA0*T.*R0w) + ...
              (R1*T.*dA1w+dA1*T.*R1w) + ...
              (R2*T.*dA2w+dA2*T.*R2w) + ...
              (R3*T.*dA3w+dA3*T.*R3w));

  % dZA = integral(r(x)*z(x)*dzdx) from x = 0 to 1
  dZA = (Z0.*dA0 + Z1.*dA1 + Z2.*dA2 + Z3.*dA3);
  dZAb = (Z0 .* dA0b+dA0 .* Z0b) + ...
	 (Z1 .* dA1b+dA1 .* Z1b) + ...
	 (Z2 .* dA2b+dA2 .* Z2b) + ...
	 (Z3 .* dA3b+dA3 .* Z3b);
  dZAw = (Z0*T.*dA0w+dA0*T.*Z0w) + ...
	 (Z1*T.*dA1w+dA1*T.*Z1w) + ...
	 (Z2*T.*dA2w+dA2*T.*Z2w) + ...
	 (Z3*T.*dA3w+dA3*T.*Z3w);

  % dAR = integral(log(r(x))*dzdx) from x = 0 to 1, where log(r) is Taylor expanded
  dAR = (log(R0)-3/2).*DZ+2*dA0./R0-dRA./R0.^2;
  dARb = ((log(R0)-3/2) .* DZb+DZ./R0 .* R0b) + ...
	 (2./R0 .* dA0b-2*dA0./R0.^2 .* R0b) - ...
	 (1./R0.^2 .* dRAb-2*dRA./R0.^3 .* R0b);
  dARw = ((log(R0)-3/2)*T.*DZw+DZ./R0*T.*R0w) + ...
	 (2./R0*T.*dA0w-2*dA0./R0.^2*T.*R0w) - ...
	 (1./R0.^2*T.*dRAw-2*dRA./R0.^3*T.*R0w);

  % Major radius of inner edge of cells indexed by JC
  rim = (JC-1.5)*dr+rg(1);

  % Integral(1*dA) from r=0 to rim
  dA0e = rim.*DZ;
  dA0eb = rim .* DZb;
  dA0ew = rim*T.*DZw;

  % Integral(R*dA) from r=0 to rim
  dRAe = rim.^2/2.*DZ;
  dRAeb = rim.^2/2 .* DZb;
  dRAew = rim.^2/2*T.*DZw;

  % Integral(Z*dA) from r=0 to rim
  dZAe = dA0e.*(Z0+DZ/2);
  dZAeb = dA0e .* (Z0b+DZb/2)+(Z0+DZ/2) .* dA0eb;
  dZAew = dA0e*T.*(Z0w+DZw/2)+(Z0+DZ/2)*T.*dA0ew;

  % Integral(dA/R) from r=0 to rim
  dARe = log(rim).*DZ;
  dAReb = log(rim) .* DZb;
  dARew = log(rim)*T.*DZw;

  % 6x6 indices into the upper outer side of 7x7 indices
  kuo = 8+[(1:6) 7+(1:6) 14+(1:6) 21+(1:6) 28+(1:6) 35+(1:6)];

  for m = 1:n-1
    i = IC(m);
    j = JC(m);
    k = i+nz*(j-1); % Line segment m:m+1 cuts through cell k
    % Contribution 1 is from Z-intervals where boundary is within cell k
     Acell(k) =  Acell(k) + dA0(m) - dA0e(m);
    RAcell(k) = RAcell(k) + dRA(m) - dRAe(m);
    ZAcell(k) = ZAcell(k) + dZA(m) - dZAe(m);
    ARcell(k) = ARcell(k) + dAR(m) - dARe(m);
    % Contribution 2 is from Z-intervals where boundary is outside cell k
     Acell(i,1:j-1) =  Acell(i,1:j-1) + dr*DZ(m);
    RAcell(i,1:j-1) = RAcell(i,1:j-1) + RA(i,1:j-1)*DZ(m)/dz;
    ZAcell(i,1:j-1) = ZAcell(i,1:j-1) + dZAe(m)*dr/rim(m);
    ARcell(i,1:j-1) = ARcell(i,1:j-1) + AR(i,1:j-1)*DZ(m)/dz;
    egg(i,j) = 1;
    % Response for contributuion 1 to psibry 
     Acellb(k) =  Acellb(k) + dA0b(m) - dA0eb(m);
    RAcellb(k) = RAcellb(k) + dRAb(m) - dRAeb(m);
    ZAcellb(k) = ZAcellb(k) + dZAb(m) - dZAeb(m);
    ARcellb(k) = ARcellb(k) + dARb(m) - dAReb(m);
    % Response for contributuion 1 to flux change at point m 
    kkk = kuo + max(-8,min(0,IZ(m)-i + (IR(m)-j)*7));
    f36 = kkk>0 & kkk<=ngg;
     Acellw(k,kkk(f36)) =  Acellw(k,kkk(f36)) + dA0w(m,f36) - dA0ew(m,f36);
    RAcellw(k,kkk(f36)) = RAcellw(k,kkk(f36)) + dRAw(m,f36) - dRAew(m,f36);
    ZAcellw(k,kkk(f36)) = ZAcellw(k,kkk(f36)) + dZAw(m,f36) - dZAew(m,f36);
    ARcellw(k,kkk(f36)) = ARcellw(k,kkk(f36)) + dARw(m,f36) - dARew(m,f36);
    % Responses for contributuion 2
    if IR0(m) == j+0.5
      % Line segment m:m+1 begins on outer edge of cell k
       Acellb(k) =  Acellb(k) + dr*Z0b(m);
      RAcellb(k) = RAcellb(k) + RA(k)/dz*Z0b(m);
      ZAcellb(k) = ZAcellb(k) + Z0(m)*dr*Z0b(m);
      ARcellb(k) = ARcellb(k) + AR(k)/dz*Z0b(m);
       Acellw(k,kkk(f36)) =  Acellw(k,kkk(f36)) + dr*Z0w(m,f36);
      RAcellw(k,kkk(f36)) = RAcellw(k,kkk(f36)) + RA(k)/dz*Z0w(m,f36);
      ZAcellw(k,kkk(f36)) = ZAcellw(k,kkk(f36)) + Z0(m)*dr*Z0w(m,f36);
      ARcellw(k,kkk(f36)) = ARcellw(k,kkk(f36)) + AR(k)/dz*Z0w(m,f36);
    end
    if IR0(m+1) == j+0.5
      % Line segment m:m+1 ends on outer edge of cell k
      kkk = kuo + max(-8,min(0,IZ(m+1)-i + (IR(m+1)-j)*7));
      f36 = kkk>0 & kkk<=ngg;
       Acellb(k) =  Acellb(k) - dr*Z0b(m+1);
      RAcellb(k) = RAcellb(k) - RA(k)/dz*Z0b(m+1);
      ZAcellb(k) = ZAcellb(k) - Z0(m+1)*dr*Z0b(m+1);
      ARcellb(k) = ARcellb(k) - AR(k)/dz*Z0b(m+1);
       Acellw(k,kkk(f36)) =  Acellw(k,kkk(f36)) - dr*Z0w(m+1,f36);
      RAcellw(k,kkk(f36)) = RAcellw(k,kkk(f36)) - RA(k)/dz*Z0w(m+1,f36);
      ZAcellw(k,kkk(f36)) = ZAcellw(k,kkk(f36)) - Z0(m+1)*dr*Z0w(m+1,f36);
      ARcellw(k,kkk(f36)) = ARcellw(k,kkk(f36)) - AR(k)/dz*Z0w(m+1,f36);
    end
  end
  ivac = ~igg&~egg;
   Acell(ivac) = 0;
  RAcell(ivac) = 0;
  ZAcell(ivac) = 0;
  ARcell(ivac) = 0;

  Ac = Acell(:);

  rgc(egg) = RAcell(egg)./Acell(egg);
  zgc(egg) = ZAcell(egg)./Acell(egg);
  rgcb(egg) = 1./Ac(egg).*RAcellb(egg) - rgc(egg)./Ac(egg).* Acellb(egg);
  zgcb(egg) = 1./Ac(egg).*ZAcellb(egg) - zgc(egg)./Ac(egg).* Acellb(egg);
  rgcw(egg,:) = 1./Ac(egg)*F.*RAcellw(egg,:) - ...
               rgc(egg)./Ac(egg)*F.* Acellw(egg,:);
  zgcw(egg,:) = 1./Ac(egg)*F.*ZAcellw(egg,:) - ...
               zgc(egg)./Ac(egg)*F.* Acellw(egg,:);

  for i = 3:nz-2 % Calculate psizrc
    for j = 3:nr-2 % Calculate psizrc
      k = i+nz*(j-1);
      if egg(k)
	kk = [17 18 19 20 24 25 26 27 31 32 33 34 38 39 40 41]; % the 4x4 in 7x7
	ir = j;
	tr = (rgc(i,j)-rg(j))/dr;
	if tr < 0
	  ir = ir-1;
	  tr = tr+1;
	  kk = kk-7;
	end
	wr0 = [1 tr tr^2 tr^3]*mx;
	wr1 = [0 1 2*tr 3*tr^2]*mx;
	wr2 = [0 0 2 6*tr]*mx;
	iz = i;
	tz = (zgc(i,j)-zg(i))/dz;
	if tz < 0
	  iz = iz-1;
	  tz = tz+1;
	  kk = kk-1;
	end
	wz0 = [1 tz tz^2 tz^3]*mx;
	wz1 = [0 1 2*tz 3*tz^2]*mx;
	wz2 = [0 0 2 6*tz]*mx;
	p = wz0*psizr(iz-1:iz+2,ir-1:ir+2)*wr0';
	pr = wz0*psizr(iz-1:iz+2,ir-1:ir+2)*wr1';
	pz = wz1*psizr(iz-1:iz+2,ir-1:ir+2)*wr0';
	prr = wz0*psizr(iz-1:iz+2,ir-1:ir+2)*wr2';
	prz = wz1*psizr(iz-1:iz+2,ir-1:ir+2)*wr1';
	pzz = wz2*psizr(iz-1:iz+2,ir-1:ir+2)*wr0';

	psizrc(i,j) = p;
	psizrcr(i,j) = pr/dr;
	psizrcz(i,j) = pz/dz;
	psizrcrr(i,j) = prr/dr/dr;
	psizrcrz(i,j) = prz/dr/dz;
	psizrczz(i,j) = pzz/dz/dz;

	psizrcb(k) = psizrcr(i,j)*rgcb(k)+psizrcz(i,j)*zgcb(k);
	psizrcw(k,:) = psizrcr(i,j)*rgcw(k,:)+psizrcz(i,j)*zgcw(k,:);
	psizrcw(k,kk) = psizrcw(k,kk)+reshape(wz0'*wr0,1,16);

	psizrcrb(k) = psizrcrr(i,j)*rgcb(k)+psizrcrz(i,j)*zgcb(k);
	psizrcrw(k,:) = psizrcrr(i,j)*rgcw(k,:)+psizrcrz(i,j)*zgcw(k,:);
	psizrcrw(k,kk) = psizrcrw(k,kk)+reshape(wz0'*wr1,1,16)/dr;

	psizrczb(k) = psizrcrz(i,j)*rgcb(k)+psizrczz(i,j)*zgcb(k);
	psizrczw(k,:) = psizrcrz(i,j)*rgcw(k,:)+psizrczz(i,j)*zgcw(k,:);
	psizrczw(k,kk) = psizrczw(k,kk)+reshape(wz1'*wr0,1,16)/dz;
      end
    end
  end
  
  % Shape parameters
  nn = (n-1)*(nr-1)+1; % Number of boundary points in R, Z
  k8 = ones(8,1);
  [R8(1),k8(1)] = max(R(1:nn-1));
  [Z8(3),k8(3)] = max(Z(1:nn-1));
  [R8(5),k8(5)] = min(R(1:nn-1));
  [Z8(7),k8(7)] = min(Z(1:nn-1));
  i8 = floor(1+(k8-1)/(nr-1)); % region
  j8 = k8 - (i8-1)*(nr-1);     % index in xs
  % Find exact values
  for k = [1 3 5 7]
    i = i8(k);
    j = 1+mod(i-2,n-1);
    x = xs(j8(k));
    if k == 1 | k == 5
      dum = R1(i)*(R1(j)+2*R2(j)+3*R3(j));
    else
      dum = Z1(i)*(Z1(j)+2*Z2(j)+3*Z3(j));
    end
    if j8(k) ~= 1 | dum > 0 % No break in derivative
      for j = 1:4
        if k == 1 | k == 5
	  dx = min(0.5/nr,max(-0.5/nr,-(R1(i)+2*R2(i)*x+3*R3(i)*x^2)/(2*R2(i)+6*R3(i)*x)));
	else
	  dx = min(0.5/nr,max(-0.5/nr,-(Z1(i)+2*Z2(i)*x+3*Z3(i)*x^2)/(2*Z2(i)+6*Z3(i)*x)));
	end
	x = x+dx;
	if x < 0
          x = x+1;
	  i = 1+mod(i-2,n-1);
	elseif x > 1
          x = x-1;
	  i = 1+mod(i,n-1);
	end
      end
      % Response to psizr
      if k == 1 | k == 5 % D[(R1(i)+2*R2(i)*x+3*R3(i)*x^2)] = 0
	xw = -(R1w(i,:)+2*R2w(i,:)*x+3*R3w(i,:)*x^2)/(2*R2(i)+6*R3(i)*x);
	xb = -(R1b(i  )+2*R2b(i  )*x+3*R3b(i  )*x^2)/(2*R2(i)+6*R3(i)*x);
      else % D[(Z1(i)+2*Z2(i)*x+3*Z3(i)*x^2)] = 0
	xw = -(Z1w(i,:)+2*Z2w(i,:)*x+3*Z3w(i,:)*x^2)/(2*Z2(i)+6*Z3(i)*x);
	xb = -(Z1b(i  )+2*Z2b(i  )*x+3*Z3b(i  )*x^2)/(2*Z2(i)+6*Z3(i)*x);
      end
    else % The extreme point is at a tip
      xw = zeros(1,36);
      xb = 0;
    end
    R8(k) = R0(i) + R1(i)*x + R2(i)*x^2 + R3(i)*x^3;
    Z8(k) = Z0(i) + Z1(i)*x + Z2(i)*x^2 + Z3(i)*x^3;
    R8p(k,I36(i,:)) = R0w(i,:)+R1w(i,:)*x+R2w(i,:)*x^2+R3w(i,:)*x^3 + (R1(i)+2*R2(i)*x+3*R3(i)*x^2)*xw;
    Z8p(k,I36(i,:)) = Z0w(i,:)+Z1w(i,:)*x+Z2w(i,:)*x^2+Z3w(i,:)*x^3 + (Z1(i)+2*Z2(i)*x+3*Z3(i)*x^2)*xw;
    R8p(k,iib) = R8p(k,iib) + (R0b(i)+R1b(i)*x+R2b(i)*x^2+R3b(i)*x^3 + (R1(i)+2*R2(i)*x+3*R3(i)*x^2)*xb)*wb;
    Z8p(k,iib) = Z8p(k,iib) + (Z0b(i)+Z1b(i)*x+Z2b(i)*x^2+Z3b(i)*x^3 + (Z1(i)+2*Z2(i)*x+3*Z3(i)*x^2)*xb)*wb;
  end
  % Squareness point k+1 is on the diagonal of square defined by k and k+2
  for k = [2 4 6 8]
    i = k-1;
    j = 1+mod(k,8);
    ratio = (Z8(j)-Z8(i))/(R8(i)-R8(j));
    m = i8(i);
    ratio1 = nan;
    for l = 1:n
      ratio2 = (Z0(m)-Z8(i))/(R0(m)-R8(j));
      if (ratio1-ratio)*(ratio2-ratio) <= 0
        break
      end
      m = 1+mod(m,n-1);
      ratio1 = ratio2;
    end
    % Approximate position of squareness point k
    R8(k) = R0(m);
    Z8(k) = Z0(m);
    x = 0; % Find exact x in interval m that makes ratio2 = ratio for point k
    for l = 1:7
      drdx = R1(m)+2*R2(m)*x+3*R3(m)*x^2;
      dzdx = Z1(m)+2*Z2(m)*x+3*Z3(m)*x^2;
      dratio2dx = dzdx/(R8(k)-R8(j)) - ratio2/(R8(k)-R8(j))*drdx;
      dx = min(0.5,max(-0.5,(ratio-ratio2)/dratio2dx));
      x = x+dx;
      if x < 0
        x = x+1;
	m = 1+mod(m-2,n-1);
      elseif x > 1
        x = x-1;
	m = 1+mod(m,n-1);
      end      
      R8(k) = R0(m) + R1(m)*x + R2(m)*x^2 + R3(m)*x^3;
      Z8(k) = Z0(m) + Z1(m)*x + Z2(m)*x^2 + Z3(m)*x^3;
      ratio2 = (Z8(k)-Z8(i))/(R8(k)-R8(j));
    end
    % Response to psizr
    ratiop = (Z8p(j,:)-Z8p(i,:))/(R8(i)-R8(j)) - ratio/(R8(i)-R8(j))*(R8p(i,:)-R8p(j,:));
    % Response of point k excluding the response of x
    R8p(k,I36(m,:)) = R0w(m,:) + R1w(m,:)*x + R2w(m,:)*x^2 + R3w(m,:)*x^3;
    R8p(k,iib) = R8p(k,iib) + (R0b(m) + R1b(m)*x + R2b(m)*x^2 + R3b(m)*x^3)*wb;
    Z8p(k,I36(m,:)) = Z0w(m,:) + Z1w(m,:)*x + Z2w(m,:)*x^2 + Z3w(m,:)*x^3;
    Z8p(k,iib) = Z8p(k,iib) + (Z0b(m) + Z1b(m)*x + Z2b(m)*x^2 + Z3b(m)*x^3)*wb;
    % The response of x makes ratiop-ratio2p = 0
    r = R1(m) + 2*R2(m)*x + 3*R3(m)*x^2;
    z = Z1(m) + 2*Z2(m)*x + 3*Z3(m)*x^2;
    % ratio2p = (Z8p(k,:)+z*xp-Z8p(i,:))/(R8(k)-R8(j)) - ratio2/(R8(k)-R8(j))*(R8p(k,:)+r*xp-R8p(j,:));
    % ratio2p excluding point k's response to x
    dum1 = 1/(R8(k)-R8(j));
    dum2 = ratio2/(R8(k)-R8(j));
    ratio2p = dum1*(Z8p(k,:)-Z8p(i,:)) - dum2*(R8p(k,:)-R8p(j,:));
    % ratiop - ratio2p - dum1*z*xp + dum2*r*xp = 0
    xp = (ratiop-ratio2p)/(dum1*z-dum2*r);
    R8p(k,:) = R8p(k,:) + r*xp;
    Z8p(k,:) = Z8p(k,:) + z*xp;
  end
  
  % Position
  rsurf = (R8(1)+R8(5))/2;
  zsurf = (Z8(3)+Z8(7))/2;
  rsurfp = (R8p(1,:)+R8p(5,:))/2;
  zsurfp = (Z8p(3,:)+Z8p(7,:))/2;

  % Size
  aminor = (R8(1)-R8(5))/2;
  bminor = (Z8(3)-Z8(7))/2;
  aminorp = (R8p(1,:)-R8p(5,:))/2;
  bminorp = (Z8p(3,:)-Z8p(7,:))/2;
  
  % Elongation
  elong = bminor/aminor;
  elongp = bminorp/aminor - elong/aminor*aminorp;
  
  % Triangularity
  triu = (rsurf-R8(3))/aminor;
  tril = (rsurf-R8(7))/aminor;
  triup = (rsurfp-R8p(3,:))/aminor - triu/aminor*aminorp;
  trilp = (rsurfp-R8p(7,:))/aminor - tril/aminor*aminorp;

  % Squareness
  Sq = zeros(4,1);
  Sqp = zeros(4,ngg);
  for k = 1:4
    
    % Calculating squareness of point R8(2*k), Z8(2*k)
    
    % Points i & j define a square, point 2*k is on a diagonal in that square
    if k == 1 | k == 3
      i = 2*k-1;
      j = 1+mod(2*k,8);
    else
      j = 2*k-1;
      i = 1+mod(2*k,8);
    end
    
    % A is the corner of the square that is inside the plasma
    rA = R8(j);
    zA = Z8(i);
    rAp = R8p(j,:);
    zAp = Z8p(i,:);

    % D is the corner of the square that is outside the plasma
    rD = R8(i);
    zD = Z8(j);
    rDp = R8p(i,:);
    zDp = Z8p(j,:);
    
    % The other corners of the square are R8(i),Z8(i) and R8(j),Z8(j) that are on the boundary
    
    % B is the boundary point on the diagonal from A to D
    rB = R8(2*k);
    zB = Z8(2*k);
    rBp = R8p(2*k,:);
    zBp = Z8p(2*k,:);
    
    % C is the point on the diagonal where an ellipse would be
    rC = rA+sqrt(0.5)*(rD-rA);
    zC = zA+sqrt(0.5)*(zD-zA);
    rCp = rAp+sqrt(0.5)*(rDp-rAp);
    zCp = zAp+sqrt(0.5)*(zDp-zAp);
    
    % Distance from A to boundary
    AB = sqrt((rB-rA).^2+(zB-zA).^2);
    ABp = (rB-rA)/AB*(rBp-rAp) + (zB-zA)/AB*(zBp-zAp);

    % Distance from A to ellipse
    AC = sqrt((rC-rA).^2+(zC-zA).^2);
    ACp = (rC-rA)/AC*(rCp-rAp) + (zC-zA)/AC*(zCp-zAp);

    % Distance from ellipse to corner
    CD = sqrt((rD-rC).^2+(zD-zC).^2);
    CDp = (rD-rC)/CD*(rDp-rCp) + (zD-zC)/CD*(zDp-zCp);
    
    % Squareness = the boundary's deviation from ellipse divided by CD
    Sq(k) = (AB-AC)/CD;
    Sqp(k,:) = (ABp-ACp)/CD - Sq(k)/CD*CDp;
  end
  squo = Sq(1);
  squi = Sq(2);
  sqli = Sq(3);
  sqlo = Sq(4);
  squop = Sqp(1,:);
  squip = Sqp(2,:);
  sqlip = Sqp(3,:);
  sqlop = Sqp(4,:);
  
  % drsep = radial distance between upper divertor separatrix and lower
  r = (R8(1)-rg(1))/dr+1;
  z = (Z8(1)-zg(1))/dz+1;
  ir = floor(r);
  iz = floor(z);
  ii = iz+(ir-1)*nz+i16';
  tr = r-ir;
  tz = z-iz;
  wz0 = [1 tz tz^2 tz^3]*mx;
  wz1 = [0 1 2*tz 3*tz^2]*mx;
  wr0 = [1 tr tr^2 tr^3]*mx;
  wr1 = [0 1 2*tr 3*tr^2]*mx;
  wr2 = [0 0 2 6*tr]*mx;
  p = wz0*psizr(iz-1:iz+2,ir-1:ir+2)*wr0';
  pr = wz0*psizr(iz-1:iz+2,ir-1:ir+2)*wr1'/dr;
  prr = wz0*psizr(iz-1:iz+2,ir-1:ir+2)*wr2'/dr^2;
  prz = wz1*psizr(iz-1:iz+2,ir-1:ir+2)*wr1'/dr/dz;
  k = 1;
  p2 = -inf;
  for i = 2:bd.n
    if s*bd.psi(i) > p2 & abs(bd.z(i)-bd.z(1)) > bminor/dz
      k = i;
      p2 = s*bd.psi(i);
    end
  end
  % k now points to competing bdef point
  p2 = bd.psi(k);
  drsep = sign(bd.z(k)-bd.z(1))*s*(p2-psibry)/pr;
  p2p = zeros(1,ngg);
  p2p(bd.ii(k,:)) = bd.w(k,:);
  prp = prr*R8p(1,:) + prz*Z8p(1,:);
  prp(ii) = prp(ii) + reshape(wz0'*wr1,1,16)/dr;
  psibryp = zeros(1,ngg);
  psibryp(iib) = wb;
  drsepp = sign(bd.z(k)-bd.z(1))*s*((p2p-psibryp)/pr - (p2-psibry)/pr^2*prp);
  
  % Integrals
  Cl = sum(dC(1:n-1));
  Vtot = sum(dRA(1:n-1))*(2*pi);
  Atot = sum(dA0(1:n-1));
  Ltot = sum(dAR(1:n-1));
  Clp(iib) = sum(dCb(1:n-1))*wb;
  Vtotp(iib) = 2*pi*sum(dRAb(1:n-1))*wb;
  Atotp(iib) = sum(dA0b(1:n-1))*wb;
  Ltotp(iib) = sum(dARb(1:n-1))*wb;
  for i = 1:n-1
    Clp(I36(i,:)) = Clp(I36(i,:)) + dCw(i,:);
    Vtotp(I36(i,:)) = Vtotp(I36(i,:)) + 2*pi*dRAw(i,:);
    Atotp(I36(i,:)) = Atotp(I36(i,:)) + dA0w(i,:);
    Ltotp(I36(i,:)) = Ltotp(I36(i,:)) + dARw(i,:);
  end

end

% Archive result

de.I36 = 'Indices used with weights of size [nbmax,36]';
re.I36 = I36;

de.R0 = 'dR0(i) = R0.w(i,:)*dpsizr(I36(i,:)'') + R0.b(i)*dpsibry';
re.R0.w = R0w;
re.R0.b = R0b;

de.R1 = 'dR1(i) = R1.w(i,:)*dpsizr(I36(i,:)'') + R1.b(i)*dpsibry';
re.R1.w = R1w;
re.R1.b = R1b;

de.R2 = 'dR2(i) = R2.w(i,:)*dpsizr(I36(i,:)'') + R2.b(i)*dpsibry';
re.R2.w = R2w;
re.R2.b = R2b;

de.R3 = 'dR3(i) = R3.w(i,:)*dpsizr(I36(i,:)'') + R3.b(i)*dpsibry';
re.R3.w = R3w;
re.R3.b = R3b;

de.Z0 = 'dZ0(i) = Z0.w(i,:)*dpsizr(I36(i,:)'') + Z0.b(i)*dpsibry';
re.Z0.w = Z0w;
re.Z0.b = Z0b;

de.Z1 = 'dZ1(i) = Z1.w(i,:)*dpsizr(I36(i,:)'') + Z1.b(i)*dpsibry';
re.Z1.w = Z1w;
re.Z1.b = Z1b;

de.Z2 = 'dZ2(i) = Z2.w(i,:)*dpsizr(I36(i,:)'') + Z2.b(i)*dpsibry';
re.Z2.w = Z2w;
re.Z2.b = Z2b;

de.Z3 = 'dZ3(i) = Z3.w(i,:)*dpsizr(I36(i,:)'') + Z3.b(i)*dpsibry';
re.Z3.w = Z3w;
re.Z3.b = Z3b;

de.dC = 'ddC(i) = dC.w(i,:)*dpsizr(I36(i,:)'') + dC.b(i)*dpsibry';
re.dC.w = dCw;
re.dC.b = dCb;

de.Yr = 'dYr(i) = Yr.w(i,:)*dpsizr(I36(i,:)'') + Yr.b(i)*dpsibry';
re.Yr.w = Yrw;
re.Yr.b = Yrb;

de.Yz = 'dYz(i) = Yz.w(i,:)*dpsizr(I36(i,:)'') + Yz.b(i)*dpsibry';
re.Yz.w = Yzw;
re.Yz.b = Yzb;

de.DR = 'dDR(i) = DR.w(i,:)*dpsizr(I36(i,:)'') + DR.b(i)*dpsibry';
re.DR.w = DRw;
re.DR.b = DRb;

de.DZ = 'dDZ(i) = DZ.w(i,:)*dpsizr(I36(i,:)'') + DZ.b(i)*dpsibry';
re.DZ.w = DZw;
re.DZ.b = DZb;    

de.T0 = 'dT0(i) = T0.w(i,:)*dpsizr(I36(i,:)'') + T0.b(i)*dpsibry';
re.T0.w = T0w;
re.T0.b = T0b;

de.T1 = 'dT1(i) = T1.w(i,:)*dpsizr(I36(i,:)'') + T1.b(i)*dpsibry';
re.T1.w = T1w;
re.T1.b = T1b;

de.dA0 = 'ddA0(i) = dA0.w(i,:)*dpsizr(I36(i,:)'') + dA0.b(i)*dpsibry';
re.dA0.w = dA0w; 
re.dA0.b = dA0b;

de.dRA = 'ddRA(i) = dRA.w(i,:)*dpsizr(I36(i,:)'') + dRA.b(i)*dpsibry';
re.dRA.w = dRAw; 
re.dRA.b = dRAb; 

de.dZA = 'ddZA(i) = dZA.w(i,:)*dpsizr(I36(i,:)'') + dZA.b(i)*dpsibry';
re.dZA.w = dZAw; 
re.dZA.b = dZAb; 

de.dAR = 'ddAR(i) = dAR.w(i,:)*dpsizr(I36(i,:)'') + dAR.b(i)*dpsibry';
re.dAR.w = dARw; 
re.dAR.b = dARb;

de.dA0e = 'ddA0e(i) = dA0e.w(i,:)*dpsizr(I36(i,:)'') + dA0e.b(i)*dpsibry';
re.dA0e.w = dA0ew; 
re.dA0e.b = dA0eb;

de.dRAe = 'ddRAe(i) = dRAe.w(i,:)*dpsizr(I36(i,:)'') + dRAe.b(i)*dpsibry';
re.dRAe.w = dRAew; 
re.dRAe.b = dRAeb;

de.dZAe = 'ddZAe(i) = dZAe.w(i,:)*dpsizr(I36(i,:)'') + dZAe.b(i)*dpsibry';
re.dZAe.w = dZAew; 
re.dZAe.b = dZAeb;

de.dARe = 'ddARe(i) = dARe.w(i,:)*dpsizr(I36(i,:)'') + dARe.b(i)*dpsibry';
re.dARe.w = dARew; 
re.dARe.b = dAReb;

de.di7x7 = 'Index differences to 7x7 points around a grid point';
re.di7x7 = di7x7;

de.Acell = 'dAcell(k) = Acell.w(k,:)*dpsizr(k+di7x7) + Acell.b(k)*dpsibry';
re.Acell.w = Acellw;
re.Acell.b = Acellb;

de.RAcell = 'dRAcell(k) = RAcell.w(k,:)*dpsizr(k+di7x7) + RAcell.b(k)*dpsibry';
re.RAcell.w = RAcellw;
re.RAcell.b = RAcellb;

de.ZAcell = 'dZAcell(k) = ZAcell.w(k,:)*dpsizr(k+di7x7) + ZAcell.b(k)*dpsibry';
re.ZAcell.w = ZAcellw;
re.ZAcell.b = ZAcellb;

de.ARcell = 'dARcell(k) = ARcell.w(k,:)*dpsizr(k+di7x7) + ARcell.b(k)*dpsibry';
re.ARcell.w = ARcellw;
re.ARcell.b = ARcellb;

de.rgc = 'drgc(k) = rgc.w(k,:)*dpsizr(k+di7x7) + rgc.b(k)*dpsibry';
re.rgc.w = rgcw;
re.rgc.b = rgcb;

de.zgc = 'dzgc(k) = zgc.w(k,:)*dpsizr(k+di7x7) + zgc.b(k)*dpsibry';
re.zgc.w = zgcw;
re.zgc.b = zgcb;

de.psizrc = 'dpsizrc(k) = psizrc.w(k,:)*dpsizr(k+di7x7) + psizrc.b(k)*dpsibry';
re.psizrc.w = psizrcw;
re.psizrc.b = psizrcb;

de.psizrcr = 'dpsizrcr(k) = psizrcr.w(k,:)*dpsizr(k+di7x7) + psizrcr.b(k)*dpsibry';
re.psizrcr.w = psizrcrw;
re.psizrcr.b = psizrcrb;

de.psizrcz = 'dpsizrcz(k) = psizrcz.w(k,:)*dpsizr(k+di7x7) + psizrcz.b(k)*dpsibry';
re.psizrcz.w = psizrczw;
re.psizrcz.b = psizrczb;

de.R8p = 'dR8 = R8p*dpsizr(:)';
re.R8p = R8p;

de.Z8p = 'dZ8 = Z8p*dpsizr(:)';
re.Z8p = Z8p;

de.rsurfp = 'drsurf = rsurfp*dpsizr(:)';
re.rsurfp = rsurfp;

de.zsurfp = 'dzsurf = zsurfp*dpsizr(:)';
re.zsurfp = zsurfp;

de.aminorp = 'daminor = aminorp*dpsizr(:)';
re.aminorp = aminorp;

de.bminorp = 'dbminor = bminorp*dpsizr(:)';
re.bminorp = bminorp;

de.elongp = 'delong = elongp*dpsizr(:)';
re.elongp = elongp;

de.triup = 'dtriu = triup*dpsizr(:)';
re.triup = triup;

de.trilp = 'dtril = trilp*dpsizr(:)';
re.trilp = trilp;

de.squop = 'dsquo = squop*dpsizr(:)';
re.squop = squop;

de.squip = 'dsqui = squip*dpsizr(:)';
re.squip = squip;

de.sqlip = 'dsqli = sqlip*dpsizr(:)';
re.sqlip = sqlip;

de.sqlop = 'dsqlo = sqlop*dpsizr(:)';
re.sqlop = sqlop;

de.drsepp = 'drsep = drsepp*dpsizr(:)';
re.drsepp = drsepp;

de.Clp = 'dCl = Clp*dpsizr(:)';
re.Clp = Clp;

de.Vtotp = 'dVtot = Vtotp*dpsizr(:)';
re.Vtotp = Vtotp;

de.Atotp = 'dAtot = Atotp*dpsizr(:)';
re.Atotp = Atotp;

de.Ltotp = 'dLtot = Ltotp*dpsizr(:)';
re.Ltotp = Ltotp;

de.iia = 'Indices used to calculate quantities at magnetic axis';
re.iia = iia;

de.drmaxisdpsi = 'drmaxis = drmaxisdpsi*dpsizr(iia)';
re.drmaxisdpsi = drmaxisdpsi;

de.dzmaxisdpsi = 'dzmaxis = dzmaxisdpsi*dpsizr(iia)';
re.dzmaxisdpsi = dzmaxisdpsi;

de.wa = 'dpsimag = wa*dpsizr(iia)';
re.wa = wa;

de.warr = 'dpsimagrr = warr*dpsizr(iia), how d2(flux)/dr2 changes at the axis';
re.warr = warr;

de.warz = 'dpsimagrz = warz*dpsizr(iia), how d2(flux)/drdz changes when psizr changes';
re.warz = warz;

de.wazz = 'dpsimagzz = wazz*dpsizr(iia), how d2(flux)/dz2 changes when psizr changes';
re.wazz = wazz;

de.iib = 'Indices used to calculate quantities at boundary-defining point';
re.iib = iib;

de.drbdefdpsi = 'drbdef = drbdefdpsi*dpsizr(iib)';
re.drbdefdpsi = drbdefdpsi;

de.dzbdefdpsi = 'dzbdef = dzbdefdpsi*dpsizr(iib)';
re.dzbdefdpsi = dzbdefdpsi;

de.wb = 'Weights used to calculate flux at boundary-defining point';
re.wb = wb;

de.wbrr = 'Weights used to calculate d2(flux)/dr2 at boundary-defining point';
re.wbrr = wbrr;

de.wbrz = 'Weights used to calculate d2(flux)/drdz at boundary-defining point';
re.wbrz = wbrz;

de.wbzz = 'Weights used to calculate d2(flux)/dz2 at boundary-defining point';
re.wbzz = wbzz;

re.info = de;

ds.plasma = 'True if plasma exists';
b.plasma = plasma;

ds.nn = 'R(1:nn), Z(1:nn) are high density boundary points, remaining are nans';
b.nn = nn;

ds.R = 'R of boundary with nr-2 interpolated points between each traced point';
b.R = R;

ds.Z = 'Z of boundary with nr-2 interpolated points between each traced point';
b.Z = Z;

ds.n = 'R0(1:n), Z0(1:n) are boundary points, remaining are nans';
b.n = n;

ds.IC = '1:st index of cell with boundary segment';
b.IC = IC;

ds.JC = '2:nd index of cell with boundary segment';
b.JC = JC;

ds.R0 = 'R of n boundary points and then nans';
b.R0 = R0;

ds.R1 = 'To find boundary between points in R0:';
b.R1 = R1;

ds.R2 = 'R = R0 + R1*x^1 + R2*x^2 + R3*x^3;';
b.R2 = R2;

ds.R3 = 'for x in range 0 to 1';
b.R3 = R3;

ds.Z0 = 'Z of n boundary points and then nans';
b.Z0 = Z0;

ds.Z1 = 'To find boundary between points in Z0:';
b.Z1 = Z1;

ds.Z2 = 'Z = Z0 + Z1*x^1 + Z2*x^2 + Z3*x^3;';
b.Z2 = Z2;

ds.Z3 = 'for x in range 0 to 1';
b.Z3 = Z3;

ds.dC = 'boundary length to next point';
b.dC = dC;

ds.Yr = 'dr*dY/dR';
b.Yr = Yr;

ds.Yz = 'dz*dY/dZ';
b.Yz = Yz;

ds.DR = 'diff(R0)';
b.DR = DR;

ds.DZ = 'diff(Z0)';
b.DZ = DZ;

ds.T0 = 'Was used to calculate R0, R1, R2, R3, Z0, Z1, Z2, Z3';
b.T0 = T0;

ds.T1 = 'Was used to calculate R0, R1, R2, R3, Z0, Z1, Z2, Z3';
b.T1 = T1;

ds.dA0 = 'Contributions to area of whole plasma from boundary segments';
b.dA0 = dA0; 

ds.dRA = 'Contributions to Integral[R*dA] over plasma from boundary segments';
b.dRA = dRA; 

ds.dZA = 'Contributions to Integral[Z*dA] over plasma from boundary segments';
b.dZA = dZA; 

ds.dAR = 'Contributions to Integral[dA/R] over plasma from boundary segments';
b.dAR = dAR; 

ds.dA0e = 'Integral[dA] for R from 0 to inner edge of cell for boundary segments';
b.dA0e = dA0e; 

ds.dRAe = 'Integral[R*dA] for R from 0 to inner edge of cell for boundary segments';
b.dRAe = dRAe; 

ds.dZAe = 'Integral[Z*dA] for R from 0 to inner edge of cell for boundary segments';
b.dZAe = dZAe; 

ds.dARe = 'Integral[dA/R] for R from 1 to inner edge of cell for boundary segments';
b.dARe = dARe; 

ds.N = 'True for points closest to x-points within dnx regions';
b.N = N;

ds.NX = 'True for points within a distance dnx of x-points';
b.NX = NX;

ds.H = 'True when Z0 on horizontal grid line in 2x2 denser grid';
b.H = H;

ds.V = 'True when R0 on vertical grid line in 2x2 denser grid';
b.V = V;

ds.Acell = 'Integral(dA) over plasma within cell';
b.Acell = Acell;

ds.RAcell = 'Integral(R*dA) over plasma within cell';
b.RAcell = RAcell;

ds.ZAcell = 'Integral(Z*dA) over plasma within cell';
b.ZAcell = ZAcell;

ds.ARcell = 'Integral(dA/R) over plasma within cell';
b.ARcell = ARcell;

ds.rgc = 'R for center of current within cell';
b.rgc = rgc;

ds.zgc = 'Z for center of current within cell';
b.zgc = zgc;

ds.psizrc = 'psizr at rgc, zgc';
b.psizrc = psizrc;

ds.psizrcr = 'd(psizr)/dr at rgc, zgc';
b.psizrcr = psizrcr;

ds.psizrcz = 'd(psizr)/dz at rgc, zgc';
b.psizrcz = psizrcz;

ds.igg = 'True for grid points covered by plasma';
b.igg = igg;

ds.egg = 'True for grid cells cut by boundary';
b.egg = egg;

ds.R8 = 'R for points at Rmax,~45,Zmax,~135,Rmin,~225,Zmin,~315';
b.R8 = R8;

ds.Z8 = 'Z for points at Rmax,~45,Zmax,~135,Rmin,~225,Zmin,~315';
b.Z8 = Z8;

ds.rsurf = 'R for geometric center';
b.rsurf = rsurf;

ds.zsurf = 'Z for geometric center';
b.zsurf = zsurf;

ds.aminor = 'Half the radial width of the plasma';
b.aminor = aminor;

ds.bminor = 'Half the vertical height of the plasma';
b.bminor = bminor;

ds.elong = 'Elongation = bminor/aminor';
b.elong = elong;

ds.triu = 'Upper triangularity';
b.triu = triu;

ds.tril = 'Lower triangularity';
b.tril = tril;

ds.squo = 'Upper outer squareness';
b.squo = squo;

ds.squi = 'Upper inner squareness';
b.squi = squi;

ds.sqli = 'Lower inner squareness';
b.sqli = sqli;

ds.sqlo = 'Lower outer squareness';
b.sqlo = sqlo;

ds.drsep = 'flux difference to competing bdef/outboard midplane gradient';
b.drsep = drsep;

ds.Cl = 'Contour length = length of the boundary';
b.Cl = Cl;

ds.Vtot = 'Integral(2*pi*R*dA) = Total volume';
b.Vtot = Vtot;

ds.Atot = 'Integral(dA) = Total area of cross section';
b.Atot = Atot;

ds.Ltot = 'Integral(dA/R)';
b.Ltot = Ltot;

ds.rmaxis = 'R of magnetic axis';
b.rmaxis = rmaxis;

ds.zmaxis = 'Z of magnetic axis';
b.zmaxis = zmaxis;

ds.emaxis = 'Elongation at magnetic axis = DZ/DR for ellipse';
b.emaxis = emaxis;

ds.kmaxis = 'Elongation = max(radius)/min(radius) for ellipse';
b.kmaxis = kmaxis;

ds.psimag = 'Flux at magnetic axis';
b.psimag = psimag;

ds.rbdef = 'R of boundary-defining point';
b.rbdef = rbdef;

ds.zbdef = 'Z of boundary-defining point';
b.zbdef = zbdef;

ds.psibry = 'Flux along boundary';
b.psibry = psibry;

ds.nulls = 'nulls in the flux';
b.nulls = nulls;

ds.touches = 'possible touch points';
b.touches = touches;

ds.bd = 'Present and prospective boundary-defining points';
b.bd = bd;

ds.resp = 'Responses to psizr and psibry';
b.resp = re;

b.info = ds;
