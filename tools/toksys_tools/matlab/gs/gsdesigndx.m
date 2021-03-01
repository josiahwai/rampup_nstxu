%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gsdesigndx
%
%  PURPOSE: Find the change of x needed to approach targets
%
%  INPUTS:  targets, weights, limits, locks - gsdesign variables
%
%  OUTPUTS: dx, the change of x (=[ic;iv;sf;sp])
%           de, fraction of flux error to remove
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	
%
%  WRITTEN BY:  Anders Welander ON 2017-07-19
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x-points
rx = (b.nulls.r(b.nulls.type=='X')-1)*dr+rg(1);
zx = (b.nulls.z(b.nulls.type=='X')-1)*dz+zg(1);
drxdpsi = b.nulls.drdy(b.nulls.type=='X',:)*dr;
dzxdpsi = b.nulls.dzdy(b.nulls.type=='X',:)*dz;
iix = b.nulls.ii(b.nulls.type=='X',:);

psizr = e.psizr;
dpsizrds = [r.dpsizrdx r.dpsizrde];

pcurrt = e.pcurrt;
dpcurrtds = [r.dpcurrtdx r.dpcurrtde];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quantities that can be targeted, locked, or limited %

% state vector
dxds = [eye(nx) zeros(nx,1)];

% coil circuit currents = x(ix.ic)
ci = x(c.ix.ic);
dcids = [eye(nic) zeros(nic,ns-nic)];

% vessel circuit currents = x(ix.iv)
vi = x(c.ix.iv);
dvids = [zeros(niv,nic) eye(niv) zeros(niv,nih+nsf+nsp+1)];

% halo current = x(ix.ih)
ih = x(c.ix.ih);
dihds = [zeros(nih,nic+niv) eye(nih) zeros(nih,nsf+nsp+1)];

% fpol parameters = x(ix.sf)
sf = x(c.ix.sf);
dsfds = [zeros(nsf,nic+niv+nih) eye(nsf) zeros(nsf,nsp+1)];

% pres parameters = x(ix.sp)
sp = x(c.ix.sp);
dspds = [zeros(nsp,nic+niv+nih+nsf) eye(nsp) zeros(nsp,1)];

% TokSys coil current
ic = c.Pcc*ci;
dicds = [c.Pcc zeros(nc,niv+nih+nsf+nsp+1)];

% TokSys vessel current
iv = c.Pvc*vi;
divds = [zeros(nv,nic) c.Pvc zeros(nv,nih+nsf+nsp+1)];

% Flux on axis
psimag = e.psimag;
dpsimagds = [r.dpsimagdx r.dpsimagde];

% Flux on boundary
psibry = e.psibry;
dpsibryds = [r.dpsibrydx r.dpsibryde];

% Plasma flux
psipla = e.psipla;
dpsiplads = [r.dpsipladx r.dpsiplade];

% R of magnetic axis
rmaxis = e.rmaxis;
drmaxisds = [r.drmaxisdx r.drmaxisde];

% Z of magnetic axis
zmaxis = e.zmaxis;
dzmaxisds = [r.dzmaxisdx r.dzmaxisde];

% R of current centroid
rcur = e.rcur;
drcurds = [r.drcurdx r.drcurde];

% Z of current centroid
zcur = e.zcur;
dzcurds = [r.dzcurdx r.dzcurde];

% Second order current moment
y2 = e.y2;
dy2ds = [r.dy2dx r.dy2de];

% Normalized second order current moment
y2n = e.y2n;
dy2nds = [r.dy2ndx r.dy2nde];

% R of 8 points used for shape parameters
R8 = b.R8;
dR8ds = b.resp.R8p*dpsizrds;

% Z of 8 points used for shape parameters
Z8 = b.Z8;
dZ8ds = b.resp.Z8p*dpsizrds;

% R of geometric center
rsurf = e.rsurf;
drsurfds = [r.drsurfdx r.drsurfde];

% Z of geometric center
zsurf = e.zsurf;
dzsurfds = [r.dzsurfdx r.dzsurfde];

% Half of maximum radial width
aminor = e.aminor;
daminords = [r.daminordx r.daminorde];

% Half of maximum vertical height
bminor = e.bminor;
dbminords = [r.dbminordx r.dbminorde];

% Elongation
elong = e.elong;
delongds = [r.delongdx r.delongde];

% Z of outermost - Z of innermost point
dzoi = Z8(1)-Z8(5);
ddzoids = dZ8ds(1,:)-dZ8ds(5,:);

% Upper triangularity
triu = e.triu;
dtriuds = [r.dtriudx r.dtriude];

% Lower triangularity
tril = e.tril;
dtrilds = [r.dtrildx r.dtrilde];

% Upper outer squareness
squo = e.squo;
dsquods = [r.dsquodx r.dsquode];

% Upper inner squareness
squi = e.squi;
dsquids = [r.dsquidx r.dsquide];

% Lower inner squareness
sqli = e.sqli;
dsqlids = [r.dsqlidx r.dsqlide];

% Lower outer squareness
sqlo = e.sqlo;
dsqlods = [r.dsqlodx r.dsqlode];

% drsep = distance between upper and lower separatrix
drsep = b.drsep;
ddrsepds = b.resp.drsepp*dpsizrds;

% Contour length = length of the boundary
Cl = e.Cl;
dClds = [r.dCldx r.dClde];

% Integral(2*pi*R*dA) = Total volume
Vtot = e.Vtot;
dVtotds = [r.dVtotdx r.dVtotde];

% Integral(dA) = Total area of cross section
Atot = e.Atot;
dAtotds = [r.dAtotdx r.dAtotde];

% Integral(dA/R)
Ltot = e.Ltot;
dLtotds = [r.dLtotdx r.dLtotde];

% Total plasma current
cpasma = e.cpasma;
dcpasmads = [r.dcpasmadx r.dcpasmade];

% normalized inductance
li = e.li;
dlids = [r.dlidx r.dlide];

% normalized inductance for elliptical plasma

% normalized inductance for circular plasma

% Higher order profile parameter

% Total thermal energy
Wth = e.Wth;
dWthds = [r.dWthdx r.dWthde];

% Normalized average pressure
betap = e.betap;
dbetapds = [r.dbetapdx r.dbetapde];

% betan
betan = e.betan;
dbetands = [r.dbetandx r.dbetande];

% q on axis

% q at normalized poloidal flux = 0.95

% fpol profile

% parallel current profile

% pressure profile

% Radial force on coils

% Radial force on coils

% Flux loop signals
fl = c.mlc*ci + c.mlv*vi + c.mpl'*pcurrt(:);
dflds = [c.mlc c.mlv zeros(c.nfl,c.nih+c.nsf+c.nsp+1)] + c.mpl'*dpcurrtds;

% Loop voltage signals

% Magnetic probe signals
bp = c.gbc*ci + c.gbv*vi + c.gpb'*pcurrt(:);
dbpds = [c.gbc c.gbv zeros(c.nbp,c.nih+c.nsf+c.nsp+1)] + c.gpb'*dpcurrtds;

% Rogowski signals
rog = c.rldata*[ci;vi;cpasma];
drogds = c.rldata*[eye(c.nic+c.niv,ns); dcpasmads];

% Flux error
fluxerror = max(abs(e.psizr_err(:)/(psimag-psibry)));
dfluxerrords = [zeros(1,nx) 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For determining range of linearity %

% fltest measures change of psibar at test points
fltest = [d.fl.test d.fl.teste];

% bdtest measures when the boundary-defining point will jump
bdtest = zeros(c.nbdtest,ns);
for i = 2:d.bd.n
  % Don't use d.bd.test because flux error removal is handled differently here
  dpsibd = d.bd.psi(i)-psibry;
  bdtest(i,:) = -(d.bd.w(i,:)*dpsizrds(d.bd.ii(i,:),:)-dpsibryds)/dpsibd;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build error vector and constraints %

% Index of item in error vector
iev = 0;

% Error vector = measured - target
ev = [];

% Error vector response to s
devds = [];

% Target names
strev = '';

% ds will hold changes to be made
ds = zeros(ns,1);

% Flag locked states
ilocked = false(ns,1);

% Lock change of flux error fraction to -1 (remove all flux error)
ds(end) = -1;
ilocked(end) = 1;

% Locks on ci and vi can be applied directly on ds
for i = 1:nic
  if ~isnan(locks.ci(i))
    ds(c.ix.ic(i)) = locks.ci(i) - ci(i);
    ilocked(i) = true;  
  end
end
for i = 1:niv
  if ~isnan(locks.vi(i))
    ds(c.ix.iv(i)) = locks.vi(i) - vi(i);
    ilocked(c.ix.iv(i)) = true;  
  end
end

% Remaining states will be used to match targets
iopen = ~ilocked;

% Number of primal variables
n = sum(iopen);

% Equality constraints matrix
Ae = zeros(0,n);

% Equality constraints r.h.s.
be = [];

% Inequility constraints matrix
Ai = zeros(0,n);

% Inequility constraints r.h.s.
bi = [];

% Will become true if any point anywhere is out of range
outofrange = false;

% Separatrix
if isfield(targets,'rsep') & isfield(targets,'zsep')
  for i = 1:min(length(targets.rsep),length(targets.zsep))
    if isfield(weights,'sep') & length(weights.sep) >= i
      weight = weights.sep(i);
    else
      weight = 1;
    end
    R = targets.rsep(i);
    Z = targets.zsep(i);
    if ~isnan(R) & ~isnan(Z) & ~isnan(weight) & R > 0
      ir = (R-rg(1))/dr+1;
      iz = (Z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration == 1
          disp(['warning in gsdesign: target rsep = ' num2str(R), ...
            ', zsep = ' num2str(Z) ' is out of range, and will be ignored'])
	end
      else
	iev = iev+1;
	w = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	ii = kz+nz*(kr-1)+i16';
	psi = w*psizr(ii);
	dpsids = w*dpsizrds(ii,:);
	ev(iev,1) = weight*(psi - psibry);
	devds(iev,:) = weight*(dpsids - dpsibryds);
	str = ['sep' num2str(i)];
	strev(iev,1:length(str)) = str;
      end
    end  
  end
end
if isfield(locks,'rsep') & isfield(locks,'zsep')
  n = min(size(locks.rsep,1),size(locks.zsep,1));
  for i = 1:n
    R = locks.rsep(i);
    Z = locks.zsep(i);
    if ~isnan(R) & ~isnan(Z) & R > 0
      ir = (R-rg(1))/dr+1;
      iz = (Z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration == 1
          disp(['warning in gsdesign: locked rsep = ' num2str(R), ...
            ', zsep = ' num2str(Z) ' is out of range, and will be ignored'])
	end
      else
        wsep = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
        iisep = kz+nz*(kr-1)+i16';
        psisep = wsep*psizr(iisep);
        dpsisepds = wsep*dpsizrds(iisep,:);
        Ae(end+1,:) = dpsisepds(iopen)-dpsibryds(iopen);
        be(end+1,1) = psibry+dpsibryds*ds-psisep-dpsisepds*ds;  
      end
    end
  end
end
if isfield(limits,'rsep') & isfield(limits,'zsep')
  n = min(size(limits.rsep,1),size(limits.zsep,1));
  for i = 1:n
    rs = linspace(limits.rsep(i,1),limits.rsep(i,2));
    zs = linspace(limits.zsep(i,1),limits.zsep(i,2));
    [is, ws] = gs_interp2(rg,zg,psizr,rs,zs,'WEIGHTS');
    ps = sum(ws'.*psizr(is)');
    [~, i1] = min(ps); % Upper limit for flux at i1 is psibry
    Ai(end+1,:) = +(ws(i1,:)*dpsizrdx(is(i1,:),:)-dpsibrydx)*dxdxcirc(:,iopen);
    bi(end+1,1) = +(psibry-ps(i1));
    [~, i2] = max(ps); % Lower limit for flux at i2 is psibry
    Ai(end+1,:) = -(ws(i2,:)*dpsizrdx(is(i2,:),:)-dpsibrydx)*dxdxcirc(:,iopen);
    bi(end+1,1) = -(psibry-ps(i2));
  end
end

% x-points
if isfield(targets,'rx') & isfield(targets,'zx')
  for i = 1:min(length(targets.rx),length(targets.zx))
    if isfield(weights,'x') & numel(weights.x) >= i
      weight = weights.x(i);
    else
      weight = 1;
    end
    R = targets.rx(i);
    Z = targets.zx(i);
    if ~isnan(R) & ~isnan(Z) & ~isnan(weight)
      ir = (R-rg(1))/dr+1;
      iz = (Z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration == 1
          disp(['warning in gsdesign: target rx = ' num2str(R), ...
            ', zx = ' num2str(Z) ' is out of range, and will be ignored'])
	end
      else
	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wr = reshape(wz0*wr1, 1, 16);
	wz = reshape(wz1*wr0, 1, 16);
	ii = kz+nz*(kr-1)+i16';
	psir = wr*psizr(ii);
	psiz = wz*psizr(ii);
	dpsirds = wr*dpsizrds(ii,:);
	dpsizds = wz*dpsizrds(ii,:);
	iev = iev+1;
	ev(iev,1) = weight*psir;
	devds(iev,:) = weight*dpsirds;
	str = ['\partial\psi_x_' num2str(i) '/\partialR'];
	strev(iev,1:length(str)) = str;
	iev = iev+1;
	ev(iev,1) = weight*psiz;
	devds(iev,:) = weight*dpsizds;
	str = ['\partial\psi_x_' num2str(i) '/\partialZ'];
	strev(iev,1:length(str)) = str;
      end
    end
  end
end
if isfield(locks,'rx') & isfield(locks,'zx')
  for i = 1:min(length(locks.rx),length(locks.zx))
    R = locks.rx(i);
    Z = locks.zx(i);
    if ~isnan(R) & ~isnan(Z) & R > 0
      ir = (R-rg(1))/dr+1;
      iz = (Z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration == 1
          disp(['warning in gsdesign: LOCKED rx = ' num2str(r), ...
            ', zx = ' num2str(z) ' is out of range, and will be ignored'])
	end
      else
 	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wr = reshape(wz0*wr1, 1, 16);
	wz = reshape(wz1*wr0, 1, 16);
	ii = kz+nz*(kr-1)+i16';
	dpsirds = wr*dpsizrds(ii,:);
	dpsizds = wz*dpsizrds(ii,:);
	psir = wr*psizr(ii) + dpsirds*ds;
	psiz = wz*psizr(ii) + dpsizds*ds;
        Ae(end+1,:) = dpsirds(:,iopen);
	be(end+1,1) = -psir;
        Ae(end+1,:) = dpsizds(:,iopen);
	be(end+1,1) = -psiz;
      end
    end
  end
end
m = numel(rx);
if m > 0
  % Polygons where x must be
  if isfield(limits,'rx') & isfield(limits,'zx')
    for i = 1:min(size(limits.rx,1),size(limits.zx,1))
      ilimx = ~isnan(limits.rx(i,:));
      rs = limits.rx(i,ilimx);
      zs = limits.zx(i,ilimx);
      ks = inpolygon(rx,zx,rs,zs);
      rsm = mean(rs(1:end-1));
      zsm = mean(zs(1:end-1));
      [~,km] = sort((rx-rsm).^2+(zx-zsm).^2);
      if any(ks)
	xPointInsidePolygon = true;
	k = km(min(find(ks(km)))); % index to x-point nearest rsm,zsm inside rs,zs
      else
	xPointInsidePolygon = false;
	k = km(1); % index to x-point nearest rsm,zsm
      end
      dRds = drxdpsi(k,:)*dpsizrds(iix(k,:),:);
      dZds = dzxdpsi(k,:)*dpsizrds(iix(k,:),:);
      dR = dRds*ds;
      dZ = dZds*ds;
      f = min([1 abs(dr/dR) abs(dz/dZ)]); % A proxy for f that will be calculated
      % Position of x-point when locked ds taken into account
      R = rx(k) + f*dR;
      Z = zx(k) + f*dZ;      
      xB = [rs(1:end-1)' rs(2:end)'];
      yB = [zs(1:end-1)' zs(2:end)'];
      if xPointInsidePolygon
	% Limit x-point displacement in 16 directions
	for j = 1:16
	  pr = cos(j*pi/8);
	  pz = sin(j*pi/8);
	  xA = R+[0 pr];
	  yA = Z+[0 pz];
	  [~, ~, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
	  if any(flag > 3 & tA > 0)
            Ai(end+1,:) = pr*dRds(:,iopen)+pz*dZds(:,iopen);
            bi(end+1,1) = max(tA(flag > 3 & tA > 0));
	  else
            Ai(end+1,:) = pr*dRds(:,iopen)+pz*dZds(:,iopen);
            bi(end+1,1) = 0;
	  end
	end
      else % x-point outside polygon limits, find closest distance
	xA = R+diff(zs)'*[0 1];
	yA = Z-diff(rs)'*[0 1];
	[rr, zz, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
	rt = [rr(flag>3)' rs];
	zt = [zz(flag>3)' zs];
	[dmin,k] = min(sqrt((rt-R).^2+(zt-Z).^2));
	if dmin > 0
	  pr = (rt(k)-R)/dmin;
	  pz = (zt(k)-Z)/dmin;
	  xA = R+[0 pr];
	  yA = Z+[0 pz];
	  [~, ~, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
	  if sum(flag > 3 & tA > 0) > 1
            Ai(end+1,:) = +pr*dRds(:,iopen)+pz*dZds(:,iopen);
            bi(end+1,1) = +max(tA(flag > 3 & tA > 0));
            Ai(end+1,:) = -pr*dRds(:,iopen)-pz*dZds(:,iopen);
            bi(end+1,1) = -min(tA(flag > 3 & tA > 0));
	    Ae(end+1,:) = +pz*dRds(:,iopen)-pr*dZds(:,iopen);
	    be(end+1,1) = 0;
	  end
	end
      end
    end
  end
  
  % Polygons where x must not be
  if isfield(limits,'rnox') & isfield(limits,'znox')
    for i = 1:min(size(limits.rnox,1),size(limits.znox,1))
      
      % Store polygon in rs, zs
      ilimx = ~isnan(limits.rnox(i,:));
      rs = limits.rnox(i,ilimx);
      zs = limits.znox(i,ilimx);
      
      % Polygon center
      rsc = (max(rs)+min(rs))/2;
      zsc = (max(zs)+min(zs))/2;
      
      % Flag x-points that are inside polygon rs, zs
      ks = inpolygon(rx,zx,rs,zs);
      
      % Sort x-points w.r.t. distance to center of the polygon
      [~,km] = sort((rx-rsc).^2+(zx-zsc).^2);
      
      % Find index (k) to nearest x-point, preferably inside polygon
      x_in_polygon = any(ks);
      if x_in_polygon
	k = km(min(find(ks(km)))); % index to x-point nearest rsc,zsc inside rs,zs
	s = +1;
      else
	k = km(1); % index to x-point nearest rsc,zsc
	s = -1;
      end
      
      % Response of nearest x-point
      dRds = drxdpsi(k,:)*dpsizrds(iix(k,:),:);
      dZds = dzxdpsi(k,:)*dpsizrds(iix(k,:),:);

      % Shift of nearest x-point due to the locked part of ds
      dR = dRds*ds;
      dZ = dZds*ds;
      
      % A proxy for f that will be calculated
      f = min([1 abs(dr/dR) abs(dz/dZ)]);
      
      % Position of x-point when locked ds taken into account
      R = rx(k) + f*dR;
      Z = zx(k) + f*dZ;
      
      % Lines that run from the x-point in directions perpendicular to polygon sides
      xA = R+diff(zs)'*[0 1];
      yA = Z-diff(rs)'*[0 1];
      
      % Lines that run along polygon sides
      xB = [rs(1:end-1)' rs(2:end)'];
      yB = [zs(1:end-1)' zs(2:end)'];
      
      % Points where lines cross
      [rr, zz, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
      
      % All points that can be the closest point
      rt = [rr(flag>3)' rs];
      zt = [zz(flag>3)' zs];
      
      % Find closest distance to edge of polygon
      [dmin,k] = min(sqrt((rt-R).^2+(zt-Z).^2));
      
      % Find unity vector toward polygon
      if dmin == 0 % Smack on polygon edge, point to center
        pr = (rsc-R)/sqrt((rsc-R)^2+(zsc-Z)^2);
	pz = (zsc-Z)/sqrt((rsc-R)^2+(zsc-Z)^2);
      else % Point along vector defined by nearest point
	pr = s*(R-rt(k))/dmin;
	pz = s*(Z-zt(k))/dmin;
      end
      
      % Limit on displacement of x along [pr,pz]
      Ai(end+1,:) = pr*dRds(:,iopen)+pz*dZds(:,iopen);
      bi(end+1,1) = -s*dmin;
      
      % Algorithm might fail if more than one x-point tries to invade the region

    end
  end
end

% Snowflakes
if isfield(targets,'rsnf') & isfield(targets,'zsnf')
  for i = 1:min(length(targets.rsnf),length(targets.zsnf))
    if isfield(weights,'snf') & i <= length(weights.snf)
      weight = weights.snf(i);
    else
      weight = 1;
    end
    R = targets.rsnf(i);
    Z = targets.zsnf(i);
    if ~isnan(R) & ~isnan(Z) & ~isnan(weight)
      ir = (R-rg(1))/dr+1;
      iz = (Z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration == 1
          disp(['warning in gsdesign: target rsnf = ' num2str(R), ...
            ', zsnf = ' num2str(Z) ' is out of range, and will be ignored'])
	end
      else
	iisnf = kz+nz*(kr-1)+i16';
	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wr2 = [0 0 2 6*tr]*mx/dr^2;
	wz2 = mx'*[0 0 2 6*tz]'/dz^2;
	wsnf = reshape(wz0*wr0, 1, 16);
	wsnfr = reshape(wz0*wr1, 1, 16);
	wsnfz = reshape(wz1*wr0, 1, 16);
	wsnfrr = reshape(wz0*wr2, 1, 16);
	wsnfrz = reshape(wz1*wr1, 1, 16);
	wsnfzz = reshape(wz2*wr0, 1, 16);
	% tsnf is the angle (theta) toward outer x-point
	if isfield(targets,'tsnf') & i <= length(targets.tsnf)
	  t = targets.tsnf(i); 
	else % Pick an angle aimed straight away from the plasma
          t = 40*sign(zmaxis-Z);
	end
	% nsnf is number of rays going out from the snowflake
	if isfield(targets,'nsnf') & i <= length(targets.nsnf)
	  nsnf = targets.nsnf(i);
	else
          nsnf = 6;
	end
	% dt is angle in degrees between snowflake rays
	dt = 360/nsnf;
	% rhosnf is a snowflake "radius"
	if isfield(targets,'rhosnf') & i <= length(targets.rhosnf)
	  rhosnf = targets.rhosnf(i);
	else
          rhosnf = (dr+dz)/2;
	end
	iev = iev+1;
	ev(iev,1) = weight*10*wsnfr*psizr(iisnf);
	devds(iev,:) = weight*10*wsnfr*dpsizrds(iisnf,:);
	strev(iev,1:10) = 'snowflakes';
	iev = iev+1;
	ev(iev,1) = weight*10*wsnfz*psizr(iisnf);
	devds(iev,:) = weight*10*wsnfz*dpsizrds(iisnf,:);
	strev(iev,1:10) = 'snowflakes';
        % Control outer x-point
	R = targets.rsnf(i)+rhosnf*cos(t*pi/180)/2;
	Z = targets.zsnf(i)+rhosnf*sin(t*pi/180)/2;
	ir = (R-rg(1))/dr+1;
	iz = (Z-zg(1))/dz+1;
	kr = max(2,min(nr-2,floor(ir)));
	kz = max(2,min(nz-2,floor(iz)));
	tr = ir-kr;
	tz = iz-kz;
	iip = kz+nz*(kr-1)+i16';
	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wp = reshape(wz0*wr0, 1, 16);
	wpr = reshape(wz0*wr1, 1, 16);
	wpz = reshape(wz1*wr0, 1, 16);
	iev = iev+1;
	ev(iev,1) = weight/10*wpr*psizr(iip);
	devds(iev,:) = weight/10*wpr*dpsizrds(iip,:);
	strev(iev,1:10) = 'snowflakes';
	iev = iev+1;
	ev(iev,1) = weight/10*wpz*psizr(iip);
	devds(iev,:) = weight/10*wpz*dpsizrds(iip,:);
	strev(iev,1:10) = 'snowflakes';
	for it = 1:nsnf % Control flux to form the snowflake
	  R = targets.rsnf(i)+rhosnf*cos((t+(it-0.5)*dt)*pi/180);
	  Z = targets.zsnf(i)+rhosnf*sin((t+(it-0.5)*dt)*pi/180);
	  ir = (R-rg(1))/dr+1;
	  iz = (Z-zg(1))/dz+1;
	  kr = max(2,min(nr-2,floor(ir)));
	  kz = max(2,min(nz-2,floor(iz)));
	  tr = ir-kr;
	  tz = iz-kz;
          iip = kz+nz*(kr-1)+i16';
          wp = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx, 1, 16);
	  iev = iev+1;
	  ev(iev,1) = weight*(wp*psizr(iip)-wsnf*psizr(iisnf));
	  devds(iev,:) = weight*(wp*dpsizrds(iip,:)-wsnf*dpsizrds(iisnf,:));
	  strev(iev,1:10) = 'snowflakes';
	end
      end
    end
  end
end
if isfield(locks,'rsnf') & isfield(locks,'zsnf')
  for i = 1:min(length(locks.rsnf),length(locks.zsnf))
    R = locks.rsnf(i);
    Z = locks.zsnf(i);
    if ~isnan(R) & ~isnan(Z) & R > 0
      ir = (R-rg(1))/dr+1;
      iz = (Z-zg(1))/dz+1;
      kr = floor(ir);
      kz = floor(iz);
      tr = ir-kr;
      tz = iz-kz;
      if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
	outofrange = true;
	if iteration == 1
          disp(['warning in gsdesign: LOCKED rsnf = ' num2str(R), ...
            ', zsnf = ' num2str(Z) ' is out of range, and will be ignored'])
	end
      else
	iisnf = kz+nz*(kr-1)+i16';
	wr0 = [1 tr tr^2 tr^3]*mx;
	wz0 = mx'*[1 tz tz^2 tz^3]';
	wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	wr2 = [0 0 2 6*tr]*mx/dr^2;
	wz2 = mx'*[0 0 2 6*tz]'/dz^2;
	wsnf = reshape(wz0*wr0, 1, 16);
	wsnfr = reshape(wz0*wr1, 1, 16);
	wsnfz = reshape(wz1*wr0, 1, 16);
	wsnfrr = reshape(wz0*wr2, 1, 16);
	wsnfrz = reshape(wz1*wr1, 1, 16);
	wsnfzz = reshape(wz2*wr0, 1, 16);
	dpsisnfrds = wsnfr*dpsizrds(iisnf,:);
	dpsisnfzds = wsnfz*dpsizrds(iisnf,:);
	psisnfr = wsnfr*psizr(iisnf) + dpsisnfrds*ds;
	psisnfz = wsnfz*psizr(iisnf) + dpsisnfzds*ds;
        Ae(end+1,:) = dpsisnfrds(:,iopen);
	be(end+1,1) = -psisnfr;
        Ae(end+1,:) = dpsisnfzds(:,iopen);
	be(end+1,1) = -psisnfz;
	% tsnf is the angle (theta) toward outer x-point
	if isfield(locks,'tsnf') & i <= length(locks.tsnf)
	  t = locks.tsnf(i); 
	else % Pick an angle aimed straight away from the plasma
          t = 40*sign(zmaxis-z);
	end
	% nsnf is number of rays going out from the snowflake
	if isfield(locks,'nsnf') & i <= length(locks.nsnf)
	  nsnf = locks.nsnf(i);
	else
          nsnf = 6;
	end
	% dt is angle in degrees between snowflake rays
	dt = 360/nsnf;
	% rhosnf is a snowflake "radius"
	if isfield(locks,'rhosnf') & i <= length(locks.rhosnf)
	  rhosnf = locks.rhosnf(i);
	else
          rhosnf = (dr+dz);
	end
	for it = 1:nsnf/2 % Control flux to form the snowflake
	  R = locks.rsnf(i)+rhosnf*cos((t+(it-0.5)*dt)*pi/180);
	  Z = locks.zsnf(i)+rhosnf*sin((t+(it-0.5)*dt)*pi/180);
	  ir = (R-rg(1))/dr+1;
	  iz = (Z-zg(1))/dz+1;
	  kr = max(2,min(nr-2,floor(ir)));
	  kz = max(2,min(nz-2,floor(iz)));
	  tr = ir-kr;
	  tz = iz-kz;
          iip = kz+nz*(kr-1)+i16';
          wp = reshape(mx'*[1 tz tz^2 tz^3]'*[1 tr tr^2 tr^3]*mx, 1, 16);
          psip = wp*psizr(iip);
          dpsipds = wp*dpsizrds(iip,:);
          Ae(end+1,:) = dpsipds(:,iopen)-dpsibryds(:,iopen);
          if fluxerror < 0.01
            be(end+1,1) = psibry-dpsibryds*ds-psip+dpsipds*ds;  
          else
	    be(end+1,1) = psibry-psip;
	  end
	end
      end
    end
  end
end

bdef.target.exists = false;
if isfield(targets,'rbdef') & isfield(targets,'zbdef') & ...
   ~isempty(targets.rbdef)  & ~isempty(targets.zbdef)
  if min(length(targets.rbdef),length(targets.zbdef)) > 1 & iteration == 1
    disp('warning in gsdesign: Only first target for rbdef, zbdef used')
  end
  if isfield(weights,'bdef') & ~isempty(weights.bdef)
    weight = weights.bdef(1);
  else
    weight = 1;
  end
  R = targets.rbdef(1);
  Z = targets.zbdef(1);
  if ~isnan(R) & ~isnan(Z) & ~isnan(weight)
    ir = (R-rg(1))/dr+1;
    iz = (Z-zg(1))/dz+1;
    kr = floor(ir);
    kz = floor(iz);
    tr = ir-kr;
    tz = iz-kz;
    xA = R+diff(c.zl')*[0 1];
    yA = Z-diff(c.rl')*[0 1];
    xB = [c.rl(1:c.nl-1)' c.rl(2:c.nl)'];
    yB = [c.zl(1:c.nl-1)' c.zl(2:c.nl)'];
    [X, Y, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
    dA = tA.*c.dl(1:c.nl-1);
    dlmin = min(abs(dA(flag > 1)));
    if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
      outofrange = true;
      if iteration == 1
        disp(['warning in gsdesign: target point rbdef = ' num2str(R), ...
          ', zbdef = ' num2str(Z) ' is out of range and will be ignored'])
      end
    else
      iev = iev+1;
      wsep = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      iisep = kz+nz*(kr-1)+i16';
      psisep = wsep*psizr(iisep);
      dpsisepds = wsep*dpsizrds(iisep,:);
      ev(iev,1) = weight*(psisep - psibry);
      devds(iev,:) = weight*(dpsisepds - dpsibryds);
      str = '\psi_b_d_e_f';
      strev(iev,1:length(str)) = str;
      % Archive data needed to make flux limits for other parts of boundary
      bdef.target.exists = true; % exists and inside limiter or very close
      bdef.target.dlmin = dlmin; % minimum distance to limiter
      bdef.target.R = R;
      bdef.target.Z = Z;
      bdef.target.w = wsep;
      bdef.target.ii = iisep;
      bdef.target.psi = psisep;
      bdef.target.dpsids = dpsisepds;            
    end
  end
end
if isfield(locks,'rbdef') & isfield(locks,'zbdef') & ...
  ~isempty(locks.rbdef) & ~isempty(locks.zbdef)
  if min(length(locks.rbdef),length(locks.zbdef)) > 1 & iteration == 1
    disp('warning in gsdesign: Only first value for locked rbdef, zbdef used')
  end
  R = locks.rbdef(1);
  Z = locks.zbdef(1);
  if ~isnan(R) & ~isnan(Z)
    ir = (R-rg(1))/dr+1;
    iz = (Z-zg(1))/dz+1;
    kr = floor(ir);
    kz = floor(iz);
    tr = ir-kr;
    tz = iz-kz;
    xA = R+diff(c.zl')*[0 1];
    yA = Z-diff(c.rl')*[0 1];
    xB = [c.rl(1:c.nl-1)' c.rl(2:c.nl)'];
    yB = [c.zl(1:c.nl-1)' c.zl(2:c.nl)'];
    [X, Y, flag, tA, tB] = line_intersection(xA, yA, xB, yB);
    dA = tA.*c.dl(1:c.nl-1);
    dlmin = min(abs(dA(flag > 1)));
    if kr < 2 | kr > nr-2 | kz < 2 | kz > nz-2
      outofrange = true;
      if iteration == 1
        disp(['warning in gsdesign: LOCKED rbdef = ' num2str(r), ...
          ', zbdef = ' num2str(z) ' is out of range, and will be ignored'])
      end
    else
      wsep = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      iisep = kz+nz*(kr-1)+i16';
      psisep = wsep*psizr(iisep);
      dpsisepds = wsep*dpsizrds(iisep,:);
      Ae(end+1,:) = dpsisepds(:,iopen)-dpsibryds(:,iopen);
      if fluxerror < 0.01
        be(end+1,1) = psibry+dpsibryds*ds-psisep-dpsisepds*ds;
      else
        be(end+1,1) = psibry-psisep;
      end
      bdef.target.exists = true; % exists and inside limiter or very close
      bdef.target.dlmin = dlmin; % minimum distance to limiter
      bdef.target.R = R;
      bdef.target.Z = Z;
      bdef.target.w = wsep;
      bdef.target.ii = iisep;
      bdef.target.psi = psisep;
      bdef.target.dpsids = dpsisepds;            
    end
  end
end
wrongway.R = [];
wrongway.Z = [];
if bdef.target.exists

  psiab = psimag - psibry;
  
  % Store target for boundary-defining point in rt, zt
  rt = bdef.target.R;
  zt = bdef.target.Z;
    
  % Decide if boundary needs to be pulled to bdef target rt, zt
  [dum, i] = min(((b.R-rt)/dr).^2+((b.Z-zt)/dz).^2);
  target_outside_plasma = [rmaxis-rt,zmaxis-zt]*[b.R(i)-rt;b.Z(i)-zt] > 0;
  far_to_boundary = (b.R(i)-rt)^2/dr^2+(b.Z(i)-zt)^2/dz^2 > 1;
  pullbnd = target_outside_plasma & far_to_boundary;

  % drsep and bdef_dpsibar limits are applied only when bdef is close to target
  far_to_bdef = (b.rbdef-rt)^2/dr^2+(b.zbdef-zt)^2/dz^2 > 9;
  if far_to_bdef
    limits.drsep = [-inf +inf];
    bdef_dpsibar = 0.01;
  else
    limits.drsep = limits_drsep;
    if isfield(limits,'bdef_dpsibar')
      bdef_dpsibar = limits.bdef_dpsibar;
    elseif isfield(limits,'drsep') & ~isempty(limits.drsep)
      bdef_dpsibar = 0;
    else
      bdef_dpsibar = 0.01;
    end
  end
             
  if pullbnd
    % Ensure monotonic flux between boundary and bdef target to pull boundary to rt, zt
    m = max(9,1+round(sqrt(dum)*2));
    rs = rt + (b.R(i)-rt)*(1:m)/m;
    zs = zt + (b.Z(i)-zt)*(1:m)/m;
    [is, ws] = gs_interp2(rg,zg,psizr,rs,zs,'WEIGHTS');
    for i = 1:m
      dpsids = ws(i,:)*dpsizrds(is(i,:),:);
      psidiff = ws(i,:)*psizr(is(i,:)') - bdef.target.psi;
      psicorr = (dpsids - bdef.target.dpsids)*ds;
      Ai(end+1,:) = -sign(psiab)*(dpsids(:,iopen)-bdef.target.dpsids(:,iopen));
      bi(end+1,1) = sign(psiab)*(psidiff+psicorr);
    end
  end
  
  % Limit flux change at critical points to get correct bdef point
  for i = 2:b.bd.n
  
    % Critical point (that could become boundary-defining)
    rn = (b.bd.r(i)-1)*dr+rg(1);
    zn = (b.bd.z(i)-1)*dz+zg(1);
    
    % Distance to boundary from critical point
    db = sqrt(min((b.R-rn).^2+(b.Z-zn).^2));
    
    % Distance to bdef target from critical point
    dt = sqrt((rt-rn)^2+(zt-zn)^2);
    
    % If critical point close to boundary and far from rt, zt then block it
    if db < 3*dr & dt > 9*dr
      dpsids = b.bd.w(i,:)*dpsizrds(b.bd.ii(i,:),:);
      psidiff = b.bd.psi(i) - bdef.target.psi;
      psicorr = (dpsids - bdef.target.dpsids)*ds;
      Ai(end+1,:) = sign(psiab)*(dpsids(:,iopen)-bdef.target.dpsids(:,iopen));
      bi(end+1,1) = -sign(psiab)*(psidiff+psicorr+psiab*bdef_dpsibar);
      wrongway.R(end+1) = rn;
      wrongway.Z(end+1) = zn;
    end
        
  end
end

if isfield(targets,'psic')
  if ~isfield(targets,'dpsicdic')
    targets.dpsicdic = 0*targets.psic;
  end
  if ~isfield(targets,'ic0')
    targets.ic0 = 0*targets.psic;
  end
  psic = mcc*ic + mcv*iv + mpc'*pcurrt(:);
  dpsicds = [mcc mcv zeros(nc,ns-nc-nv)] + mpc'*dpcurrtds;
  for i = 1:min([length(targets.psic), ...
                 length(targets.dpsicdic), ...
                 length(targets.ic0), nc])
    if isfield(weights,'psic') & length(weights.psic) >= i
      weight = weights.psic(i);
    else
      weight = 1;
    end
    target = targets.psic(i) + targets.dpsicdic(i)*(ic(i)-targets.ic0(i));
    if ~isnan(target) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(psic(i) - target);
      devds(iev,:) = weight*dpsicds(i,:);
      devds(iev,i) = devds(iev,i)-weight*targets.dpsicdic(i);
      str = ['\psic' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

if isfield(targets,'psiv')
  if ~isfield(targets,'dpsivdiv')
    targets.dpsivdiv = 0*targets.psiv;
  end
  targets.dpsivdiv = targets.dpsivdiv(:);
  if numel(targets.dpsivdiv) < nv
    targets.dpsivdiv(nv) = 0;
  end
  if numel(targets.dpsivdiv) > nv
    targets.dpsivdiv = targets.dpsivdiv(1:nv);
  end
  if ~isfield(targets,'iv0')
    targets.iv0 = 0*targets.psiv;
  end
  psiv = mcv'*ic + mvv*iv + mpv'*pcurrt(:);
  dpsivds = [mcv' mvv-diag(targets.dpsivdiv) zeros(nv,ns-nc-nv)] + mpv'*dpcurrtds;
  for i = 1:min([length(targets.psiv), ...
                 length(targets.dpsivdiv), ...
                 length(targets.iv0), nv])
    if isfield(weights,'psiv') & length(weights.psiv) >= i
      weight = weights.psiv(i);
    else
      weight = 1;
    end
    target = targets.psiv(i) + targets.dpsivdiv(i)*(iv(i)-targets.iv0(i));
    if ~isnan(target) & ~isnan(weight)
      iev = iev+1;
      ev(iev,1) = weight*(psiv(i) - target);
      devds(iev,:) = weight*dpsivds(i,:);
      devds(iev,nc+i) = devds(iev,nc+i)-weight*targets.dpsivdiv(i);
      str = ['\psiv' num2str(i)];
      strev(iev,1:length(str)) = str;
    end
  end
end

% Flux expansion
if isfield(targets,'fluxexp') & isfield(targets,'rfluxexp') & isfield(targets,'zfluxexp')
  for i = 1:min([numel(targets.rfluxexp),numel(targets.zfluxexp),numel(targets.fluxexp)])   
    
    % target and point where flux expansion is measured
    tfe = targets.fluxexp(i);
    rfe = targets.rfluxexp(i);
    zfe = targets.zfluxexp(i);

    % dfluxexp is distance in outboard midplane to another separatrix
    if isfield(targets,'dfluxexp') & length(targets.dfluxexp) >= i
      dfe = targets.dfluxexp(i);
    else
      dfe = 0;
    end

    if ~isnan(tfe) & ~isnan(rfe) &~isnan(zfe) & ~isnan(dfe)

      % Flux expansion found by calling calc_fluxexp with equilibrium structure eqfe
      eqfe.rg = rg;
      eqfe.zg = zg;
      eqfe.rbbbs = b.R;
      eqfe.zbbbs = b.Z;
      eqfe.nbbbs = sum(~isnan(b.R));
      eqfe.psizr = psizr;
      eqfe.psibry = psibry;
      [fe, rfe2, zfe2] = calc_fluxexp(eqfe, rfe, zfe, dfe);

      % The response is found numerically, mark all relevant grid points with a nan
      dfedpsizr = zeros(1,ngg);
      kr = floor((b.R8(1)-rg(1))/dr+1);
      kz = floor((b.Z8(1)-zg(1))/dz+1);
      dfedpsizr(kz+nz*(kr-1)+i16) = nan;
      kz = floor((zfe-zg(1))/dz+1);
      kr = floor((rfe-rg(1))/dr+1);
      dfedpsizr(kz+nz*(kr-1)+i16) = nan;
      if dfe > 0
        kz = floor((zfe2-zg(1))/dz+1);
        kr = floor((rfe2-rg(1))/dr+1);
        dfedpsizr(kz+nz*(kr-1)+i16) = nan;
      end

      % Call calc_fluxexp with small flux perturbations at all relevant grid points
      dpsiba = psibry-psimag;
      for j = 1:ngg
	if isnan(dfedpsizr(j))
          eqfe.psizr(j) = psizr(j)+dpsiba/1e6; % Perturb this grid point
          dum = calc_fluxexp(eqfe, rfe, zfe, dfe);
          dfedpsizr(j) = (dum-fe)/dpsiba*1e6; % Store response
	  eqfe.psizr(j) = psizr(j); % Restore this grid point           
	end
      end
      dfeds = dfedpsizr*dpsizrds;

      % Include in error vector
      if isfield(weights,'fluxexp')
	weight = weights.fluxexp;
      else
	weight = 1;
      end
      iev = iev+1;
      ev(iev,1) = weight*(fe - tfe);
      devds(iev,:) = weight*dfeds;
      str = ['\psiexp'];
      strev(iev,1:length(str)) = str;
    end
  end
end
if isfield(locks,'fluxexp') & isfield(locks,'rfluxexp') & isfield(locks,'zfluxexp')
  for i = 1:min([numel(locks.rfluxexp),numel(locks.zfluxexp),numel(locks.fluxexp)])   
    
    % lock and point where flux expansion is measured
    lfe = locks.fluxexp(i);
    rfe = locks.rfluxexp(i);
    zfe = locks.zfluxexp(i);

    % dfluxexp is distance in outboard midplane to another separatrix
    if isfield(locks,'dfluxexp') & length(locks.dfluxexp) >= i
      dfe = locks.dfluxexp(i);
    else
      dfe = 0;
    end

    if ~isnan(lfe) & ~isnan(rfe) &~isnan(zfe) & ~isnan(dfe)

      % Flux expansion found by calling calc_fluxexp with equilibrium structure eqfe
      eqfe.rg = rg;
      eqfe.zg = zg;
      eqfe.rbbbs = b.R;
      eqfe.zbbbs = b.Z;
      eqfe.nbbbs = sum(~isnan(b.R));
      eqfe.psizr = psizr;
      eqfe.psibry = psibry;
      [fe, rfe2, zfe2] = calc_fluxexp(eqfe, rfe, zfe, dfe);

      % The response is found numerically, mark all relevant grid points with a nan
      dfedpsizr = zeros(1,ngg);
      kr = floor((b.R8(1)-rg(1))/dr+1);
      kz = floor((b.Z8(1)-zg(1))/dz+1);
      dfedpsizr(kz+nz*(kr-1)+i16) = nan;
      kz = floor((zfe-zg(1))/dz+1);
      kr = floor((rfe-rg(1))/dr+1);
      dfedpsizr(kz+nz*(kr-1)+i16) = nan;
      if dfe > 0
        kz = floor((zfe2-zg(1))/dz+1);
        kr = floor((rfe2-rg(1))/dr+1);
        dfedpsizr(kz+nz*(kr-1)+i16) = nan;
      end

      % Call calc_fluxexp with small flux perturbations at all relevant grid points
      dpsiba = psibry-psimag;
      for j = 1:ngg
	if isnan(dfedpsizr(j))
          eqfe.psizr(j) = psizr(j)+dpsiba/1e6; % Perturb this grid point
          dum = calc_fluxexp(eqfe, rfe, zfe, dfe);
          dfedpsizr(j) = (dum-fe)/dpsiba*1e6; % Store response
	  eqfe.psizr(j) = psizr(j); % Restore this grid point           
	end
      end
      dfeds = dfedpsizr*dpsizrds;

      % Include in equality constraints
      Ae(end+1,:) = dfeds(iopen);
      be(end+1,1) = lfe-fe-dfeds*ds;
      
    end
  end
end

% Directly accessible objects
name = fields(direct);
for i = 1:numel(name)
  nam = char(name(i));
  v = eval(nam);
  dvds = eval(['d' nam 'ds']);
  
  % Targets
  if isfield(targets,nam)
    target = getfield(targets,nam);
    if isfield(weights,nam)
      weight = getfield(weights,nam);
      weight(numel(weight)+1:numel(target)) = 1;
    else
      weight = ones(numel(target),1);
    end
    for j = 1:min(numel(v),numel(target))
      if weight(j) ~= 0 & ~isnan(weight(j)) & ~isnan(target(j)) & ~isnan(v(j))
	iev = iev+1;
	ev(iev,1) = weight(j)*(v(j) - target(j));
	devds(iev,:) = weight(j)*dvds(j,:);
	if numel(v) > 1
	  str = [nam num2str(j)];
	else
	  str = nam;
	end
	strev(iev,1:numel(str)) = str;
      end
    end
  end
  
  % Locks
  if isfield(locks,nam)
    lock = getfield(locks,nam);
    for j = 1:min(numel(v),numel(lock))
      Ae(end+1,:) = dvds(j,iopen);
      be(end+1,1) = lock(j)-v(j)-dvds(j,:)*ds;
    end
  end

  % Limits
  if isfield(limits,nam)
    lim = getfield(limits,nam);
    for j = 1:min(numel(v),size(lim,1))
      Ai(end+1,:) = -dvds(j,iopen);
      bi(end+1,1) = -(lim(j,1)-v(j));
      Ai(end+1,:) = +dvds(j,iopen);
      bi(end+1,1) = +(lim(j,2)-v(j));
    end
  end

end

strev(strev==0) = 32;

if isempty(ev)
  iev = iev+1;
  ev(iev,1) = fluxerror;
  devds(iev,1:ns) = dfluxerrords;
  strev(iev,1:10) = '\psi_e_r_r';  
end

evc = ev + devds*ds;

% d(error vector)/d(primal variable)
devdpv = devds(:,iopen);

% Hessian for devdpv
H = 2*devdpv'*devdpv;
if ~any(H(:))
  H = eye(n);
end

% Hessian r.h.s.
h = (2*evc'*devdpv)';

if ~isnan(xlock.rbdef(1))
  % Two boundary-defining points, lock the flux of the second

  % Find index in b.bd of the second bdef point
  [~,i] = min((b.bd.r-xlock.rbdef(1)).^2+(b.bd.z-xlock.zbdef(1)).^2);
  [~,j] = min((b.bd.r-xlock.rbdef(2)).^2+(b.bd.z-xlock.zbdef(2)).^2);
  if i == 1
    i = j;
  end
  % Now i points to the second bdef point in b.bd
  
  % Remove second bdef point from bdtest
  bdtest(i,:) = 0;
  
  % Lock the flux of the second bdef point to that of the first  
  w = b.bd.w(i,:);  
  ii = b.bd.ii(i,:)';
  psi = b.bd.psi(i);
  dpsids = w*dpsizrds(ii,:);
  Ae(end+1,:) = dpsids(:,iopen)-dpsibryds(:,iopen);
  be(end+1,1) = psibry+dpsibryds*ds-psi-dpsids*ds;
  % Modifying this after ds is known and calculating new xqp was not worth pursuing

end
if isfield(locks,'psic')
  for i = 1:min(length(locks.psic),nc)
    Ae(end+1,:) = dpsicds(:,iopen);
    be(end+1,1) = locks.psic(i)-psic(i)-dpsicds(i,:)*ds;
  end
end
if isfield(locks,'psiv')
  for i = 1:min(length(locks.psiv),nv)
    Ae(end+1,:) = dpsivds(:,iopen);
    be(end+1,1) = locks.psiv(i)-psiv(i)-dpsivds(i,:)*ds;
  end
end
ie = ~isnan(be);

if isfield(limits,'psic')
  n = min(size(limits.psic,1),nc);
  for i = 1:n
    k = indic(i);
    dpsicdpv = dysdx(k,:)*dxdxcirc(:,iopen);
    Ai(end+1,:) = -dpsicdpv;
    bi(end+1,1) = -(limits.psic(i,1)-ys(k));
    Ai(end+1,:) = +dpsicdpv;
    bi(end+1,1) = +(limits.psic(i,2)-ys(k));
  end
end
if isfield(limits,'psiv')
  n = min(size(limits.psiv,1),nv);
  for i = 1:n
    k = indiv(i);
    dpsivdpv = dysdx(k,:)*dxdxcirc(:,iopen);
    Ai(end+1,:) = -dpsivdpv;
    bi(end+1,1) = -(limits.psiv(i,1)-ys(k));
    Ai(end+1,:) = +dpsivdpv;
    bi(end+1,1) = +(limits.psiv(i,2)-ys(k));
  end
end
ii = any(Ai,2) & ~isinf(bi) & ~isnan(bi);

if outofrange & iteration == 1
  disp(['Allowed r,z are: ' num2str(rg(2)) ' < r < ' num2str(rg(nr-1)), ...
    ' & ' num2str(zg(2)) ' < z < ' num2str(zg(nz-1))])
end


% Time to optimize

if license('checkout','optimization_toolbox')
  opts = optimset(...
    'Algorithm','interior-point-convex',...
    'Display','off',...
    'TolX',1e-12,...
    'TolFun',1e-12);
  [xqp,fqp,eqp] = quadprog((H+H')/2,h,Ai(ii,:),bi(ii),Ae(ie,:),be(ie),[],[],[],opts);
  if eqp == -6
    % Non-convex problem detected
    if iteration == 1
      disp('Warning gsdesign: spec does not contain enough information for unique solution')
      disp('Looking for solution close to the initial equilibrium')
    end
    % Lock all targets and minimize weighted sum of changes to the primal variables
    wp = sum(abs(dR8ds(:,iopen))+abs(dZ8ds(:,iopen)))';
    Ae = [Ae(ie,:); devds(:,iopen)];
    be = [be(ie); -evc];
    H = diag(wp);
    h = zeros(sum(iopen),1);
    [xqp,fqp,eqp] = quadprog((H+H')/2,h,Ai(ii,:),bi(ii),Ae,be,[],[],[],opts);
  end
  ds(iopen) = xqp;  
else
  for i = 1:length(be)
    Ai(end+1,:) = +Ae(i,:);
    bi(end+1,1) = +be(i,1);
    Ai(end+1,:) = -Ae(i,:);
    bi(end+1,1) = -be(i,1);
  end
  ii = any(Ai,2) & ~isinf(bi);
  rep = pdipmqpneq2(H,h,Ai(ii,:),bi(ii),100,1e-9,0.97);
  ds(iopen) = rep.x;
end

dev = devds*ds;

dpsizr = reshape(dpsizrds*ds,nz,nr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decide what fraction of ds that can be made with a close to linear response %

% flmax1 < 1 if change of normalized poloidal flux < dpsibarlimit
flmax1 = max(abs(fltest*ds)); % < 1 for linear approx

% The fraction of ds that can be made without deforming the plasma too much
f = min(1,1/flmax1);

% if it was concluded that two points should define the boundary then this may no longer be true if f < 1
if f < 1
  xlock.rbdef = nan(1,2);
  xlock.zbdef = nan(1,2);
end

% Archive bdef point [floating index units]
rbdefs(iteration) = b.bd.r(1);
zbdefs(iteration) = b.bd.z(1);

% bdmax1 < 1 if the boundary-defining point does not jump
[bdmax1, j] = max(bdtest*ds); % > 1 means bdef point jumps, bdmax1(1) = 0 always

% The fraction of ds that can be made before the boundary-defining point jumps
fbd = min(1,1/bdmax1);

% Check if f is limited by a jump of the boundary-defining point
run_gsdesigndx_again = false;
if fbd < f
  % Check if this is a second flip back to where the bdef point was in the previous iteration
  if iteration > 1 & controlled_bdef_switch(iteration-1)
    % The bdef point also switched in the previous iteration
    rp = rbdefs(iteration-1);
    zp = zbdefs(iteration-1);
    rj = b.bd.r(j);
    zj = b.bd.z(j);
    if abs(rp-rj) < 0.1 & abs(zp-zj) < 0.1
      % The bdef point in iteration+1 will switch back to what it was in iteration-1
      % The minimum in the cost function occurs when both points are on the boundary
      xlock.rbdef = b.bd.r([1 j]);
      xlock.zbdef = b.bd.z([1 j]);
      run_gsdesigndx_again = true;
    end
  end
  % Update f to what it can be with approximately linear response
  f = fbd;
  % Estimated displacements of bdef point (1) and new bdef point (j) [floating index units]
  dr1 = b.bd.drdpsi(1,:)*dpsizr(b.bd.ii(1,:)')*f;
  dz1 = b.bd.dzdpsi(1,:)*dpsizr(b.bd.ii(1,:)')*f;
  drj = b.bd.drdpsi(j,:)*dpsizr(b.bd.ii(1,:)')*f;
  dzj = b.bd.dzdpsi(j,:)*dpsizr(b.bd.ii(1,:)')*f;
  % Estimated new positions of bdef point (1) and new bdef point (j) [floating index units]
  r1 = b.bd.r(1) + dr1;
  z1 = b.bd.z(1) + dz1;
  rj = b.bd.r(j) + drj;
  zj = b.bd.z(j) + dzj;
  if max(abs([dr1,dz1,drj,dzj])) < 0.25
    % The change f*dpsizr moves both old bdef (1) and new (j) by small amount
    % Tune f until 1e-7 < (yj-y1)/(psimag-psibry) < 1e-6
    for k = 1:9
      psizrx = psizr + dpsizr*f;
      nulls = gsnulls(psizrx,4*c.nbdtest);
      touches = gstouches(c,psizrx);
      rbs = [nulls.r touches.r];
      zbs = [nulls.z touches.z];
      wbs = [nulls.w;touches.w];
      ybs = [nulls.y touches.psi];
      ibs = [nulls.ii;touches.ii];
      [~,i1] = min((rbs-r1).^2+(zbs-z1).^2);
      [~,ij] = min((rbs-rj).^2+(zbs-zj).^2);
      % Flux at points 1 and j
      y1 = ybs(i1); % flux at 1
      yj = ybs(ij); % flux at j
      if abs(rbs(i1)-r1) < 0.1 & abs(zbs(i1)-z1) < 0.1 & abs(rbs(ij)-rj) < 0.1 & abs(zbs(ij)-zj) < 0.1
	% Points 1 and j still exist and are in expected regions
	% Make (yj-y1)/(psimag-psibry) < 1e-6
	dy1df = wbs(i1,:)*dpsizr(ibs(i1,:)');
	dyjdf = wbs(ij,:)*dpsizr(ibs(ij,:)');
	% (dyjdf-dy1df)*df = 5e-7*(psimag-psibry)-(yj-y1)
	df = (5e-7*(psimag-psibry)-(yj-y1))/(dyjdf-dy1df);
	f = min(1,max(0,f+df));
      end
      dum = (yj-y1)/(psimag-psibry);
      if dum > 1e-7 & dum < 1e-6
        controlled_bdef_switch(iteration) = 1;
        break
      end
    end
    % Now f is tuned to switch bdef point to rbs(ij), zbs(ij)
  end
end

dx = f*ds(1:end-1);
de = -f*ds(end);
%save(['step' num2str(iteration)])
