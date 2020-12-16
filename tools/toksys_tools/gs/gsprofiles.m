function p = gsprofiles(e,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE: p = gsprofiles(e)
%         p = gsprofiles(e,b), where b = gsboundary(...)
%
%  PURPOSE: Calculate profiles for Grad-Shafranov equilibrium
%
%  INPUTS: e, structure with fields:
%             psizr = flux on grid [Wb]
%             rg, zg = grid coordinates [m]
%             rbbbs, zbbbs, nbbbs = boundary coordinates
%             rmaxis, zmaxis = axis coordinate
%             fpol, ffprim, pres, pprime
%             xlim, ylim = limiter coordinates (R,Z values)
%          b, output from gsboundary (for response calculation)
%
%  OUTPUTS: p, structure with profiles versus normalized flux
%              p.info contains detailed information
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  WRITTEN BY:  Anders Welander ON 2017-01-13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
  plasma = e.rmaxis > 0;
else
  plasma = b.plasma;
end

include_response = logical(1);

mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % For interpolation
mu0 = 4e-7*pi;

rg = e.rg(:)';
dr = rg(2)-rg(1);
nr = numel(rg);
zg = e.zg(:);
dz = zg(2)-zg(1);
nz = numel(zg);
psizr = e.psizr;
pres = e.pres(:)';
fpol = e.fpol(:)';
pprime = e.pprime(:)';
ffprim = e.ffprim(:)';
fprime = ffprim./fpol;
rmaxis = e.rmaxis;
zmaxis = e.zmaxis;
m = numel(fpol);
rgg = ones(nz,1)*rg;
zgg = zg*ones(1,nr);
ngg = nz*nr;
% Number of angles with contour points
if isfield(e,'nbbbs')
  n = e.nbbbs-1;
else
  n = numel(e.rbbbs)-1;
end
nm = n*m;

% Normalized flux values to contour
if isfield(e,'psibar') & numel(e.psibar) == m
  psibar = e.psibar(:)';
else
  psibar = linspace(0,1,m);
end


if plasma


  % The tilted ellipses near the axis must be treated numerically
  % Calculate small contour around the axis instead of just a point
  psibarc = psibar; % psibarc is for contours
  ismaxis = psibar == 0;
  psibaraxis = min(1e-6,min(psibar(~ismaxis))/1e4);
  psibarc(ismaxis) = psibaraxis;

  % contour points at boundary
  Rb = reshape(e.rbbbs(1:n),n,1);
  Zb = reshape(e.zbbbs(1:n),n,1);

  % The points wrap around in positive theta direction
  if sum((Rb+Rb([2:n 1])).*(Zb([2:n 1])-Zb)) < 0
    Rb(1:n) = Rb(n:-1:1);
    Zb(1:n) = Zb(n:-1:1);
  end

  % Find boundary point with weakest field and use as rbdef, zbdef
  rb = (Rb-rg(1))/dr+1;
  zb = (Zb-zg(1))/dz+1;
  jb = floor(rb);
  ib = floor(zb);
  rb = rb-jb;
  zb = zb-ib;
  kb = ib+nz*(jb-1);
  IIG = [kb-(nz+1)   kb-nz     kb-(nz-1)   kb-(nz-2) ...
	 kb-1        kb        kb+1        kb+2      ...
	 kb+(nz-1)   kb+nz     kb+(nz+1)   kb+(nz+2) ...
	 kb+(2*nz-1) kb+(2*nz) kb+(2*nz+1) kb+(2*nz+2)];
  Gr0 = [ones(n,1) rb rb.^2 rb.^3]*mx;
  Gz0 = [ones(n,1) zb zb.^2 zb.^3]*mx;
  Gr1 = [zeros(n,1) ones(n,1) 2*rb 3*rb.^2]*mx;
  Gz1 = [zeros(n,1) ones(n,1) 2*zb 3*zb.^2]*mx;
  G = ...
  [Gz0(:,1).*Gr0(:,1) Gz0(:,2).*Gr0(:,1) Gz0(:,3).*Gr0(:,1) Gz0(:,4).*Gr0(:,1) ...
   Gz0(:,1).*Gr0(:,2) Gz0(:,2).*Gr0(:,2) Gz0(:,3).*Gr0(:,2) Gz0(:,4).*Gr0(:,2) ...
   Gz0(:,1).*Gr0(:,3) Gz0(:,2).*Gr0(:,3) Gz0(:,3).*Gr0(:,3) Gz0(:,4).*Gr0(:,3) ...
   Gz0(:,1).*Gr0(:,4) Gz0(:,2).*Gr0(:,4) Gz0(:,3).*Gr0(:,4) Gz0(:,4).*Gr0(:,4)];
  Gr = ...
  [Gz0(:,1).*Gr1(:,1) Gz0(:,2).*Gr1(:,1) Gz0(:,3).*Gr1(:,1) Gz0(:,4).*Gr1(:,1) ...
   Gz0(:,1).*Gr1(:,2) Gz0(:,2).*Gr1(:,2) Gz0(:,3).*Gr1(:,2) Gz0(:,4).*Gr1(:,2) ...
   Gz0(:,1).*Gr1(:,3) Gz0(:,2).*Gr1(:,3) Gz0(:,3).*Gr1(:,3) Gz0(:,4).*Gr1(:,3) ...
   Gz0(:,1).*Gr1(:,4) Gz0(:,2).*Gr1(:,4) Gz0(:,3).*Gr1(:,4) Gz0(:,4).*Gr1(:,4)];
  Gz = ...
  [Gz1(:,1).*Gr0(:,1) Gz1(:,2).*Gr0(:,1) Gz1(:,3).*Gr0(:,1) Gz1(:,4).*Gr0(:,1) ...
   Gz1(:,1).*Gr0(:,2) Gz1(:,2).*Gr0(:,2) Gz1(:,3).*Gr0(:,2) Gz1(:,4).*Gr0(:,2) ...
   Gz1(:,1).*Gr0(:,3) Gz1(:,2).*Gr0(:,3) Gz1(:,3).*Gr0(:,3) Gz1(:,4).*Gr0(:,3) ...
   Gz1(:,1).*Gr0(:,4) Gz1(:,2).*Gr0(:,4) Gz1(:,3).*Gr0(:,4) Gz1(:,4).*Gr0(:,4)];
  yb = sum(G.*psizr(IIG),2);
  ybr = sum(Gr.*psizr(IIG),2);
  ybz = sum(Gz.*psizr(IIG),2);
  yg = sqrt(ybr.^2+ybz.^2);
  ygmax = max(yg);
  % ibdef is index to boundary-defining point
  if isfield(e,'ibdef')
    ibdef = e.ibdef;
  else
    [~,ibdef] = min(yg);
  end
  % jbdef is the index before ibdef with wrap-around
  if ibdef == 1
    jbdef = n;
  else
    jbdef = ibdef-1;
  end
  % kbdef is the index after ibdef with wrap-around
  if ibdef == n
    kbdef = 1;
  else
    kbdef = ibdef+1;
  end
  iib = IIG(ibdef,:);
  isx = yg/ygmax < 1e-4;
  rbdef = Rb(ibdef);
  zbdef = Zb(ibdef);
  psibry = yb(ibdef);
  r1 = (rbdef-rg(1))/dr+1;
  z1 = (zbdef-zg(1))/dz+1;
  j = floor(r1);
  i = floor(z1);
  tr = r1-j;
  tz = z1-i;
  wr0 = [1 tr tr^2 tr^3]*mx;
  wz0 = [1 tz tz^2 tz^3]*mx;
  wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
  wz1 = [0 1 2*tz 3*tz^2]*mx/dz;
  wr2 = [0 0 2 6*tr]*mx/dr^2;
  wz2 = [0 0 2 6*tz]*mx/dz^2;
  wr3 = [0 0 0 6]*mx/dr^3;
  wz3 = [0 0 0 6]*mx/dz^3;
  ybrr = wz0*psizr(i-1:i+2,j-1:j+2)*wr2';
  ybrz = wz1*psizr(i-1:i+2,j-1:j+2)*wr1';
  ybzz = wz2*psizr(i-1:i+2,j-1:j+2)*wr0';
  ybrrr = wz0*psizr(i-1:i+2,j-1:j+2)*wr3';
  ybrrz = wz1*psizr(i-1:i+2,j-1:j+2)*wr2';
  ybrzz = wz2*psizr(i-1:i+2,j-1:j+2)*wr1';
  ybzzz = wz3*psizr(i-1:i+2,j-1:j+2)*wr0';
  wb  =  reshape(wz0'*wr0,1,16);
  wbr =  reshape(wz0'*wr1,1,16);
  wbz =  reshape(wz1'*wr0,1,16);
  wbrr =  reshape(wz0'*wr2,1,16);
  wbrz =  reshape(wz1'*wr1,1,16);
  wbzz =  reshape(wz2'*wr0,1,16);

  % flux on axis
  for count = 1:9 % Iterate as needed to find axis more exactly
    r1 = (rmaxis-rg(1))/dr+1;
    z1 = (zmaxis-zg(1))/dz+1;
    j = floor(r1);
    i = floor(z1);
    tr = r1-j;
    tz = z1-i;
    wr0 = [1 tr tr^2 tr^3]*mx;
    wz0 = [1 tz tz^2 tz^3]*mx;
    wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
    wz1 = [0 1 2*tz 3*tz^2]*mx/dz;
    wr2 = [0 0 2 6*tr]*mx/dr^2;
    wz2 = [0 0 2 6*tz]*mx/dz^2;
    psimag = wz0*psizr(i-1:i+2,j-1:j+2)*wr0';
    yar = wz0*psizr(i-1:i+2,j-1:j+2)*wr1';
    yaz = wz1*psizr(i-1:i+2,j-1:j+2)*wr0';
    yarr = wz0*psizr(i-1:i+2,j-1:j+2)*wr2';
    yarz = wz1*psizr(i-1:i+2,j-1:j+2)*wr1';
    yazz = wz2*psizr(i-1:i+2,j-1:j+2)*wr0';
    drza = -[yarr yarz; yarz yazz]\[yar; yaz];
    s = max(abs(drza./[dr;dz]));
    if s < 1e-12
      break
    end
    t = min(0.1/s,1);
    rmaxis = rmaxis+t*drza(1);
    zmaxis = zmaxis+t*drza(2);
  end
  i16 = [-nz-1 -nz -nz+1 -nz+2 -1 0 1 2 nz-1 nz nz+1 nz+2 2*nz-1 2*nz 2*nz+1 2*nz+2];
  iia = i+nz*(j-1) + i16;
  wa  =  reshape(wz0'*wr0,1,16);

  % Find the usually small tilt of the ellipse at the axis
  v = linspace(0,pi,8);
  rv = cos(v);
  zv = sin(v);
  D2 = yarr*rv.^2/2+yarz*rv.*zv+yazz*zv.^2/2;
  [~,k] = max(D2);
  D3 = (yazz-yarr)*rv.*zv+yarz*(rv.^2-zv.^2); % d(d2)/dv should be 0
  D4 = (yazz-yarr)*(rv.^2-zv.^2)-4*yarz*rv.*zv;
  u = v(k)-D3(k)/D4(k);
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
  kmaxis = sqrt(abs(yavv/yauu));
  jmaxis = rmaxis*e.pprime(1)+e.ffprim(1)/rmaxis/mu0;
  q0 = abs(fpol(1)/mu0/jmaxis)*(1/kmaxis+kmaxis)/rmaxis^2;

  % Integrations along field lines are in positive phi direction
  ipdir = -sign(yauu);
  btdir = sign(fpol(1));
  swdir = -ipdir*btdir; % switch direction of dRs, dZs if swdir == -1

  % Normalized poloidal flux
  dba = psibry-psimag;
  psibarzr = (psizr-psimag)/dba;

  % Poloidal flux at contours
  psip = psimag + dba*psibar;

  % A vector of zeros and a vector of ones
  O = zeros(nm,1);
  l = ones(nm,1);

  % R, Z for all contour points in floating index units
  R0 = (rmaxis + (Rb-rmaxis)*psibarc - rg(1))/dr + 1;
  Z0 = (zmaxis + (Zb-zmaxis)*psibarc - zg(1))/dz + 1;
  Ra = (rmaxis - rg(1))/dr + 1;
  Za = (zmaxis - zg(1))/dz + 1;

  % Cosine and Sine components
  RH = sqrt((R0-Ra).^2+(Z0-Za).^2);
  CO = (R0(:,end)-Ra)./RH(:,end)*ones(1,m);
  SI = (Z0(:,end)-Za)./RH(:,end)*ones(1,m);

  % Target normalized flux for contour points
  YT = ones(n,1)*psibarc;

  % Define sizes
  Y = zeros(n,m);
  Yr = zeros(n,m);
  Yz = zeros(n,m);
  Yrr = zeros(n,m);
  Yrz = zeros(n,m);
  Yzz = zeros(n,m);

  J = floor(R0(:));
  I = floor(Z0(:));
  K = I+nz*(J-1);
  Iw = [K-(nz+1)   K-nz     K-(nz-1)   K-(nz-2) ...
	K-1        K        K+1        K+2      ...
	K+(nz-1)   K+nz     K+(nz+1)   K+(nz+2) ...
	K+(2*nz-1) K+(2*nz) K+(2*nz+1) K+(2*nz+2)];  
  TR = R0(:)-J;
  TZ = Z0(:)-I;
  Wr0 = [l TR TR.^2 TR.^3]*mx;
  Wz0 = [l TZ TZ.^2 TZ.^3]*mx;
  Wr1 = [O l 2*TR 3*TR.^2]*mx;
  Wz1 = [O l 2*TZ 3*TZ.^2]*mx;
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
  Y(:) = sum(W.*psibarzr(Iw),2);
  % linear interpolation
  for i = 1:n
    rh = RH(i,:);
    for j = 1:m
      k = 1;
      while k < m-1 & Y(i,k+1) < YT(i,j)
	k = k+1;
      end
      % Now the correct value for RH(i,j) is inside range rh(k:k+1)
      dY = Y(i,k+1)-Y(i,k);
      f = min(1,max(0,(Y(i,k+1)-YT(i,j))/dY));
      RH(i,j) = f*rh(k) + (1-f)*rh(k+1);
    end
  end
  R0 = Ra + RH.*CO;
  Z0 = Za + RH.*SI;

  count = 0;
  dRHmax = 1;
  while dRHmax > 1e-9 & count < 50
    count = count + 1;
    J = floor(R0(:));
    I = floor(Z0(:));
    K = I+nz*(J-1);
    Iw = [K-(nz+1)   K-nz     K-(nz-1)   K-(nz-2) ...
	  K-1        K        K+1        K+2      ...
	  K+(nz-1)   K+nz     K+(nz+1)   K+(nz+2) ...
	  K+(2*nz-1) K+(2*nz) K+(2*nz+1) K+(2*nz+2)];  
    TR = R0(:)-J;
    TZ = Z0(:)-I;
    Wr0 = [l TR TR.^2 TR.^3]*mx;
    Wz0 = [l TZ TZ.^2 TZ.^3]*mx;
    Wr1 = [O l 2*TR 3*TR.^2]*mx;
    Wz1 = [O l 2*TZ 3*TZ.^2]*mx;
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
    Y(:) = sum(W.*psibarzr(Iw),2);
    % Newton-Rhapson
    Yr(:) = sum(Wr.*psibarzr(Iw),2);
    Yz(:) = sum(Wz.*psibarzr(Iw),2);
    dRH = (YT-Y)./(CO.*Yr+SI.*Yz);
    if psibarc(end) == 1
      dRH(isx,end) = 0; % Yr, Yz are 0 or practically 0
    end
    dRH(dRH>dr/4) = dr/4;
    dRH(dRH<-dr/4) = -dr/4;
    dRH(RH+dRH<0) = 0;
    RH = RH + dRH;
    dRHmax = max(max(abs(dRH)));
    R0 = Ra + RH.*CO;
    Z0 = Za + RH.*SI;
  end
  Wr2 = [O O 2*l 6*TR]*mx;
  Wz2 = [O O 2*l 6*TZ]*mx;
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
  Yrr(:) = sum(Wrr.*psibarzr(Iw),2);
  Yrz(:) = sum(Wrz.*psibarzr(Iw),2);
  Yzz(:) = sum(Wzz.*psibarzr(Iw),2);

  % Find contour points R4, Z4 consisting of extremes in R or Z
  [~,io] = max(R0);
  [~,iu] = max(Z0);
  [~,ii] = min(R0);
  [~,il] = min(Z0);
  R4 = zeros(4,m);
  Z4 = zeros(4,m);
  for j = 1:m
    R4(:,j) = R0([io(j) iu(j) ii(j) il(j)],j);
    Z4(:,j) = Z0([io(j) iu(j) ii(j) il(j)],j);
  end
  O4 = zeros(4*m,1);
  l4 = ones(4*m,1);
  YT4 = ones(4,1)*psibarc;
  Y4 = zeros(4,m);
  Y4r = zeros(4,m);
  Y4z = zeros(4,m);
  Y4rr = zeros(4,m);
  Y4rz = zeros(4,m);
  Y4zz = zeros(4,m);
  count = 0;
  dxymax = 1;
  while dxymax > 1e-9 & count < 50
    count = count + 1;
    J4 = floor(R4(:));
    I4 = floor(Z4(:));
    K4 = I4+nz*(J4-1);
    II4 = [K4-(nz+1)   K4-nz     K4-(nz-1)   K4-(nz-2) ...
	   K4-1        K4        K4+1        K4+2      ...
	   K4+(nz-1)   K4+nz     K4+(nz+1)   K4+(nz+2) ...
	   K4+(2*nz-1) K4+(2*nz) K4+(2*nz+1) K4+(2*nz+2)];  
    ER = R4(:)-J4;
    EZ = Z4(:)-I4;
    Er0 = [l4 ER ER.^2 ER.^3]*mx;
    Ez0 = [l4 EZ EZ.^2 EZ.^3]*mx;
    Er1 = [O4 l4 2*ER 3*ER.^2]*mx;
    Ez1 = [O4 l4 2*EZ 3*EZ.^2]*mx;
    Er2 = [O4 O4 2*l4 6*ER]*mx;
    Ez2 = [O4 O4 2*l4 6*EZ]*mx;
    W4 = ...
    [Ez0(:,1).*Er0(:,1) Ez0(:,2).*Er0(:,1) Ez0(:,3).*Er0(:,1) Ez0(:,4).*Er0(:,1) ...
     Ez0(:,1).*Er0(:,2) Ez0(:,2).*Er0(:,2) Ez0(:,3).*Er0(:,2) Ez0(:,4).*Er0(:,2) ...
     Ez0(:,1).*Er0(:,3) Ez0(:,2).*Er0(:,3) Ez0(:,3).*Er0(:,3) Ez0(:,4).*Er0(:,3) ...
     Ez0(:,1).*Er0(:,4) Ez0(:,2).*Er0(:,4) Ez0(:,3).*Er0(:,4) Ez0(:,4).*Er0(:,4)];
    W4r = ...
    [Ez0(:,1).*Er1(:,1) Ez0(:,2).*Er1(:,1) Ez0(:,3).*Er1(:,1) Ez0(:,4).*Er1(:,1) ...
     Ez0(:,1).*Er1(:,2) Ez0(:,2).*Er1(:,2) Ez0(:,3).*Er1(:,2) Ez0(:,4).*Er1(:,2) ...
     Ez0(:,1).*Er1(:,3) Ez0(:,2).*Er1(:,3) Ez0(:,3).*Er1(:,3) Ez0(:,4).*Er1(:,3) ...
     Ez0(:,1).*Er1(:,4) Ez0(:,2).*Er1(:,4) Ez0(:,3).*Er1(:,4) Ez0(:,4).*Er1(:,4)];
    W4z = ...
    [Ez1(:,1).*Er0(:,1) Ez1(:,2).*Er0(:,1) Ez1(:,3).*Er0(:,1) Ez1(:,4).*Er0(:,1) ...
     Ez1(:,1).*Er0(:,2) Ez1(:,2).*Er0(:,2) Ez1(:,3).*Er0(:,2) Ez1(:,4).*Er0(:,2) ...
     Ez1(:,1).*Er0(:,3) Ez1(:,2).*Er0(:,3) Ez1(:,3).*Er0(:,3) Ez1(:,4).*Er0(:,3) ...
     Ez1(:,1).*Er0(:,4) Ez1(:,2).*Er0(:,4) Ez1(:,3).*Er0(:,4) Ez1(:,4).*Er0(:,4)];
    W4rr = ...
    [Ez0(:,1).*Er2(:,1) Ez0(:,2).*Er2(:,1) Ez0(:,3).*Er2(:,1) Ez0(:,4).*Er2(:,1) ...
     Ez0(:,1).*Er2(:,2) Ez0(:,2).*Er2(:,2) Ez0(:,3).*Er2(:,2) Ez0(:,4).*Er2(:,2) ...
     Ez0(:,1).*Er2(:,3) Ez0(:,2).*Er2(:,3) Ez0(:,3).*Er2(:,3) Ez0(:,4).*Er2(:,3) ...
     Ez0(:,1).*Er2(:,4) Ez0(:,2).*Er2(:,4) Ez0(:,3).*Er2(:,4) Ez0(:,4).*Er2(:,4)];
    W4rz = ...
    [Ez1(:,1).*Er1(:,1) Ez1(:,2).*Er1(:,1) Ez1(:,3).*Er1(:,1) Ez1(:,4).*Er1(:,1) ...
     Ez1(:,1).*Er1(:,2) Ez1(:,2).*Er1(:,2) Ez1(:,3).*Er1(:,2) Ez1(:,4).*Er1(:,2) ...
     Ez1(:,1).*Er1(:,3) Ez1(:,2).*Er1(:,3) Ez1(:,3).*Er1(:,3) Ez1(:,4).*Er1(:,3) ...
     Ez1(:,1).*Er1(:,4) Ez1(:,2).*Er1(:,4) Ez1(:,3).*Er1(:,4) Ez1(:,4).*Er1(:,4)];
    W4zz = ...
    [Ez2(:,1).*Er0(:,1) Ez2(:,2).*Er0(:,1) Ez2(:,3).*Er0(:,1) Ez2(:,4).*Er0(:,1) ...
     Ez2(:,1).*Er0(:,2) Ez2(:,2).*Er0(:,2) Ez2(:,3).*Er0(:,2) Ez2(:,4).*Er0(:,2) ...
     Ez2(:,1).*Er0(:,3) Ez2(:,2).*Er0(:,3) Ez2(:,3).*Er0(:,3) Ez2(:,4).*Er0(:,3) ...
     Ez2(:,1).*Er0(:,4) Ez2(:,2).*Er0(:,4) Ez2(:,3).*Er0(:,4) Ez2(:,4).*Er0(:,4)];

    Y4(:) = sum(W4.*psibarzr(II4),2);
    Y4r(:) = sum(W4r.*psibarzr(II4),2);
    Y4z(:) = sum(W4z.*psibarzr(II4),2);
    Y4rr(:) = sum(W4rr.*psibarzr(II4),2);
    Y4rz(:) = sum(W4rz.*psibarzr(II4),2);
    Y4zz(:) = sum(W4zz.*psibarzr(II4),2);
    % Adjust R4([1 3],:),Z4([1 3],:) by x,y a few times, with x,y found from:
    % pz + prz*x + pzz*y = 0    (d(psi)/dz = 0)
    % p + pr*x + pz*y = psibar
    % Solving for x: p*pzz-pz*pz + (pr*pzz-pz*prz)*x = psibar*pzz
    % x = (psibar*pzz+pz*pz-p*pzz)/(pr*pzz-pz*prz)
    % Solving for y: p*prz-pr*pz + (pz*prz-pr*pzz)*y = psibar*prz
    % y = -(psibar*prz+pr*pz-p*prz)/(pr*pzz-pz*prz)
    d4 = Y4r.*Y4zz-Y4z.*Y4rz;
    x4 = (YT4.*Y4zz+Y4z.*Y4z-Y4.*Y4zz)./d4;
    y4 = -(YT4.*Y4rz+Y4r.*Y4z-Y4.*Y4rz)./d4;
    % Adjust R4([2 4],:),Z4([2 4],:) by x,y a few times, with x,y found from:
    % pr + prr*x + prz*y = 0    (d(psi)/dr = 0)
    % p + pr*x + pz*y = psibar
    % Solving for x: p*prz-pr*pz + (pr*prz-pz*prr)*x = psibar*prz
    % x = (psibar*prz+pr*pz-p*prz)/(pr*prz-pz*prr)
    % Solving for y: p*prr-pr*pr + (pz*prr-pr*prz)*y = psibar*prr
    % y = -(psibar*prr+pr*pr-p*prr)/(pr*prz-pz*prr)
    d4([2 4],:) = Y4r([2 4],:).*Y4rz([2 4],:)-Y4z([2 4],:).*Y4rr([2 4],:);
    x4([2 4],:) = (YT4([2 4],:).*Y4rz([2 4],:)+Y4r([2 4],:).*Y4z([2 4],:)-...
                    Y4([2 4],:).*Y4rz([2 4],:))./d4([2 4],:);
    y4([2 4],:) = -(YT4([2 4],:).*Y4rr([2 4],:)+Y4r([2 4],:).*Y4r([2 4],:)-...
                    Y4([2 4],:).*Y4rr([2 4],:))./d4([2 4],:);
    if psibar(1) == 0
      x4(:,1) = 0;
      y4(:,1) = 0;
    end
    yg4 = sqrt(Y4r(:,end).^2+Y4z(:,end).^2);
    isx4 = yg4/ygmax < 1e-4;
    x4(isx4,end) = 0;
    y4(isx4,end) = 0;
    s4 = 0.25./sqrt(x4.^2+y4.^2);
    s4(s4>1) = 1;
    x4 = s4.*x4; % To decrease risk of convergence problem
    y4 = s4.*y4; % sqrt(x4.^2+y4.^2) <= 0.25
    R4 = R4 + x4;
    Z4 = Z4 + y4;
    dxymax = max([x4(:); y4(:)]);
  end

  % Change to physics units
  R0 = (R0-1)*dr+rg(1);
  Z0 = (Z0-1)*dz+zg(1);
  R4 = (R4-1)*dr+rg(1);
  Z4 = (Z4-1)*dz+zg(1);
  Wr = Wr/dr;
  Wz = Wz/dz;
  Yr = Yr*(dba/dr);
  Yz = Yz*(dba/dz);
  Yrr = Yrr*(dba/dr/dr);
  Yrz = Yrz*(dba/dr/dz);
  Yzz = Yzz*(dba/dz/dz);
  Y4r = Y4r*(dba/dr);
  Y4z = Y4z*(dba/dz);
  Y4rr = Y4rr*(dba/dr/dr);
  Y4rz = Y4rz*(dba/dr/dz);
  Y4zz = Y4zz*(dba/dz/dz);
  
  % Magnetic field
  Br = -Yz./(2*pi*R0);     % Radial component
  Bz = +Yr./(2*pi*R0);     % Vertical component
  Bt = ones(n,1)*fpol./R0; % Toroidal component
  Bp = sqrt(Br.^2+Bz.^2);  % Poloidal component
  B  = sqrt(Bp.^2+Bt.^2);  % Magnitude of B
  
  % Current density, j = fprime/mu0*B + [0,R*pprime,0], for [r,t,z] components
  fm = ones(n,1)*fprime/mu0;
  jr = fm.*Br;                        % Radial component
  jz = fm.*Bz;                        % Vertical component
  jt = fm.*Bt + ones(n,1)*pprime.*R0; % Toroidal component

  R = R0([1:end 1],:);
  Z = Z0([1:end 1],:);
  DR = diff(R);
  DZ = diff(Z);
  
  Yg2 = Yr.^2+Yz.^2;
  Yg = sqrt(Yg2);

  % Morph to linear interpolation when direction to next point changes > 75.5225 degrees
  DD = sqrt(DR.^2+DZ.^2);
  F = 4*(DR.*DR([2:n 1],:) + DZ.*DZ([2:n 1],:))./DD./DD([2:n 1],:);
  F(F>1) = 1; % Cubic interpolation when F = 1, direction change < 75.5225 degrees
  F(F<0) = 0; % Linear interpolation when F = 0, direction change > 90 degrees
  
  % How field lines deviate from straight lines at beginning and end of intervals  
  T0 = (Yr.*DR + Yz.*DZ)./(Yr.*DZ - Yz.*DR).*F;
  T1 = (Yr([2:end 1],:).*DR + Yz([2:end 1],:).*DZ)./...
       (Yr([2:end 1],:).*DZ - Yz([2:end 1],:).*DR).*F;

  % Coefficients for interpolation between points on the same contour
  R1 = DR-T0.*DZ;
  R2 = (2*T0+T1).*DZ;
  R3 = -(T0+T1).*DZ;
  Z1 = DZ+T0.*DR;
  Z2 = -(2*T0+T1).*DR;
  Z3 = (T0+T1).*DR;
  

  % Integrations along field lines done numerically

  dC = zeros(n,m); % distances between points along contour
  dq = zeros(n,m); % toroidal turn between contour points
  dS = zeros(n,m); % field-line distances between contour points
  jdSfp = zeros(n,m); % Integral(j*dS) / fprime
  jdSpp = zeros(n,m); % Integral(j*dS) / pprime
  jdSff = zeros(n,m); % Integral(j*dS) / ffprim
  jdSfpYr1 = zeros(n,m); % d(jdSfp)/dYr1
  jdSfpYr2 = zeros(n,m); % d(jdSfp)/dYr2
  jdSfpYz1 = zeros(n,m); % d(jdSfp)/dYz1
  jdSfpYz2 = zeros(n,m); % d(jdSfp)/dYz2
  xs = linspace(0,1,nr); % angular subdivisions for total of n*(nr-1)
  xm = (xs(1:nr-1)+xs(2:nr))/2;
  x1 = [1/(2*pi*mu0)*ones(n,1)*(1-xm)]';
  x2 = [1/(2*pi*mu0)*ones(n,1)*xm]';
  for j = 1:m
    Rs = [R0(:,j)*ones(1,nr-1)+R1(:,j)*xm+R2(:,j)*xm.^2+R3(:,j)*xm.^3]';
    % dRs, dZs point along B if ipdir == -1, otherwise opposite to B
    dRs = swdir*diff((R1(:,j)*xs+R2(:,j)*xs.^2+R3(:,j)*xs.^3)');
    dZs = swdir*diff((Z1(:,j)*xs+Z2(:,j)*xs.^2+Z3(:,j)*xs.^3)');
    Ygs = [Yg(:,j)*(1-xm) + Yg([2:n 1],j)*xm]'; % Decent approximation
    dCs = sqrt(dRs.^2+dZs.^2);
    % If Bt > 0, integrate along B, otherwise in opposite direction
    dTs = btdir*fpol(j)*2*pi*dCs./Ygs;
    dSs = sqrt(dCs.^2+dTs.^2);
    dC(:,j) = sum(dCs)';
    dq(:,j) = sum(dCs./Ygs./Rs)*abs(fpol(j));
    dS(:,j) = sum(dSs);
    % Yrs = [Yr(:,j)*(1-xm) + Yr([2:n 1],j)*xm]'; % Decent approximation
    % Yzs = [Yz(:,j)*(1-xm) + Yz([2:n 1],j)*xm]'; % Decent approximation
    % jr = fprime*Br/mu0, jz = fprime*Bz/mu0, jt = R*pprime+ffprim/R/mu0
    % Integrals of jpar along field in positive phi direction, 3 contributions
    % jdSfp(:,j) = sum((Yrs.*dZs-Yzs.*dRs)./Rs)/(2*pi*mu0);
    jdSpp(:,j) = sum(dTs.*Rs);
    jdSff(:,j) = sum(dTs./Rs)/mu0;
    jdSfpYr1(:,j) =  sum(x1.*dZs./Rs);
    jdSfpYr2(:,j) =  sum(x2.*dZs./Rs);
    jdSfpYz1(:,j) = -sum(x1.*dRs./Rs);
    jdSfpYz2(:,j) = -sum(x2.*dRs./Rs);
  end
  jdSfp = jdSfpYr1.*Yr+jdSfpYr2.*Yr([2:n 1],:)+jdSfpYz1.*Yz+jdSfpYz2.*Yz([2:n 1],:);
  
  % For calculating flux surface averages
  Bpmax = max(Bp);
  Bpx = Bp; % Like Bp but never extremely small
  for j = 1:m
    weak = Bp(:,j) < Bpmax(j)/50;
    Bpx(weak,j) = Bpmax(j)/50;
  end
  dl = (dC(1:end,:)+dC([2:end 1],:))/2;
  g = dl./Bpx./[ones(n,1)*sum(dl./Bpx)];
  
  % q = sum(dl./Bpx./(2*pi*R0.^2)).*abs(fpol);
  
  
  % Surface integrals done analytically

  % dA0 = integral(r(x)*dzdx) from x = 0 to 1
  dA0 =     R0.*DZ + 1/2*R1.*Z1 + 2/3*R1.*Z2 + 3/4*R1.*Z3 + 1/3*R2.*Z1 + ...
	2/4*R2.*Z2 + 3/5*R2.*Z3 + 1/4*R3.*Z1 + 2/5*R3.*Z2 + 3/6*R3.*Z3;

  % dA1 = integral(x*r(x)*dzdx) from x = 0 to 1
  dA1 = 1/2*R0.*Z1 + 2/3*R0.*Z2 + 3/4*R0.*Z3 + 1/3*R1.*Z1 + 2/4*R1.*Z2 + 3/5*R1.*Z3 + ...
	1/4*R2.*Z1 + 2/5*R2.*Z2 + 3/6*R2.*Z3 + 1/5*R3.*Z1 + 2/6*R3.*Z2 + 3/7*R3.*Z3;

  % dA2 = integral(x^2*r(x)*dzdx) from x = 0 to 1
  dA2 = 1/3*R0.*Z1 + 2/4*R0.*Z2 + 3/5*R0.*Z3 + 1/4*R1.*Z1 + 2/5*R1.*Z2 + 3/6*R1.*Z3 + ...
	1/5*R2.*Z1 + 2/6*R2.*Z2 + 3/7*R2.*Z3 + 1/6*R3.*Z1 + 2/7*R3.*Z2 + 3/8*R3.*Z3;

  % dA3 = integral(x^3*r(x)*dzdx) from x = 0 to 1
  dA3 = 1/4*R0.*Z1 + 2/5*R0.*Z2 + 3/6*R0.*Z3 + 1/5*R1.*Z1 + 2/6*R1.*Z2 + 3/7*R1.*Z3 + ...
	1/6*R2.*Z1 + 2/7*R2.*Z2 + 3/8*R2.*Z3 + 1/7*R3.*Z1 + 2/8*R3.*Z2 + 3/9*R3.*Z3;

  % dRA = integral(r(x)^2/2*dzdx) from x = 0 to 1
  dRA = (R0.*dA0 + R1.*dA1 + R2.*dA2 + R3.*dA3)/2;

  % dZA = integral(r(x)*z(x)*dzdx) from x = 0 to 1
  dZA = (Z0.*dA0 + Z1.*dA1 + Z2.*dA2 + Z3.*dA3);

  % dAR = integral(log(r(x))*dzdx) from x = 0 to 1, where log(r) is Taylor expanded
  dAR = (log(R0)-3/2).*DZ+2*dA0./R0-dRA./R0.^2;

  % Integrate (r1+(r2-r1)*x)^2/2 = x*r1^2/2+x^2/2*r1*(r2-r1)+x^3/6*(r2-r1)^2
  % r1^2/2+r1*(r2-r1)/2+(r2-r1)^2/6
  % (r1^2+r1*r2+r2^2)/6
  dV = [[(rmaxis^2+rmaxis*R(1:n,1)+R(1:n,1).^2)/6.*(Z(1:n,1)-zmaxis)+...
    dRA(:,1)+...
    (rmaxis^2+rmaxis*R(2:n+1,1)+R(2:n+1,1).^2)/6.*(zmaxis-Z(2:n+1,1))] ...
    [(R(1:n,1:m-1).^2 + R(1:n,1:m-1).*R(1:n,2:m) + R(1:n,2:m).^2)/6.* ...
     (Z(1:n,2:m)-Z(1:n,1:m-1)) + ...
     dRA(:,2:m) - dRA(:,1:m-1) + ...
     (R(2:n+1,1:m-1).^2+R(2:n+1,1:m-1).*R(2:n+1,2:m)+R(2:n+1,2:m).^2)/6.* ...
     (Z(2:n+1,1:m-1)-Z(2:n+1,2:m))]]*(2*pi);
  % sum(dRA)*2*pi = cumsum(sum(dV))
  if ismaxis(1)
    dV(:,2) = dV(:,2)+dV(:,1);
    dV(:,1) = 0;
  end
  
  % Calculate Bp^2
  Bp2 = (Yr.^2+Yz.^2)./(2*pi*R0).^2;

  % Calculate Bp^2 in middle of 4 contour points * dV
  B2dV = (Bp2+Bp2([2:n 1],:)+Bp2(:,[2:m 1])+Bp2([2:n 1],[2:m 1])).*dV/4;
  B2V = cumsum(sum(B2dV));
  
  % Toroidal field at all points
  Bt = ones(n,1)*fpol./R0;
  
  % Total field at all points
  Btot2 = Bt.^2 + Bp2;
  Btot = sqrt(Btot2);
    
  % Poloidal field^2 at 4 points
  B4p2 = (Y4r.^2+Y4z.^2)./(2*pi*R4).^2;

  % Toroidal field at 4 points
  B4t = ones(4,1)*fpol./R4;
  
  % Total field at 4 points
  B4tot2 = B4t.^2 + B4p2;
  B4tot = sqrt(B4tot2);
  
  Bmax = max(B4tot);
  
  V = sum(dRA)*2*pi; V(ismaxis) = 0;
  A = sum(dA0);      A(ismaxis) = 0;
  L = sum(dAR);      L(ismaxis) = 0;
  C = sum(dC);       C(ismaxis) = 0;
  q = sum(dq);
  S = sum(dS);
  
  % <Btot^2>
  B2 = sum(g.*Btot2);
  
  % <Btot>
  Babs = sum(g.*Btot);
  
  % Y.R.Lin-Liu and R.L.Miller, Phys. Plasmas 2 (5), May 1995 p.1666, Equation 1
  
  % Fraction circulating particles, fc = 3/4*h2av*Integral(lm*dlm/den)
  fc = zeros(1,m);
  lm = [1:2:99]/100;
  dlm = 2/100;
  for j = 1:m
    h2av = B2(j)/Bmax(j)^2;
    h = Btot(:,j)/Bmax(j);
    lmIntegral = 0;
    for k =1:50
      den = sum(sqrt(1-lm(k)*h).*dS(:,j))/S(j);
      lmIntegral = lmIntegral + lm(k)*dlm/den;
    end
    fc(j) = 3/4*h2av*lmIntegral;
  end
  
  % Trapped particle fraction
  ft = 1-fc;
   
  rsurf = (R4(1,:)+R4(3,:))/2;
  zsurf = (Z4(2,:)+Z4(4,:))/2;
  aminor = (R4(1,:)-R4(3,:))/2;
  bminor = (Z4(2,:)-Z4(4,:))/2;
  elong = bminor./aminor;
  tril = (rsurf-R4(4,:))./aminor;
  triu = (rsurf-R4(2,:))./aminor;

  dWth = [1.5*pres(1)*V(1), 0.75*(pres(1:m-1)+pres(2:m)).*diff(V)];
  Wth = cumsum(dWth); Wth(ismaxis) = 0;

  dIt = [pprime(1)*V(1)/(2*pi) + ffprim(1)*L(1)/mu0, ...
	(pprime(1:m-1)+pprime(2:m)).*diff(V)/(4*pi) + ...
	(ffprim(1:m-1)+ffprim(2:m)).*diff(L)/(2*mu0)];
  It = cumsum(dIt);

  dpsit = [fpol(1)*L(1), (fpol(1:m-1)+fpol(2:m)).*diff(L)/2];
  psit = cumsum(dpsit);

  rhot = sqrt(abs(psit)/abs(psit(end)));

  % Helical flux that drives parallel current
  psis = psip+psit./q;

  betap = 4/3/mu0*(Wth./V-1.5*pres)./(It./C).^2;
  if psibar(1) == 0
    % Wth(r) = Integral(1.5*(p0+p'*yarr*r^2/2)*2*pi*r*dr*2*pi*rmaxis)
    % V(r) = Integral(2*pi*r*dr*2*pi*rmaxis) = pi*r^2
    % Wth(r)/V(r) = 1.5*p0+1.5*p'*yarr*r^2/4
    % Wth(r)/V(r)-1.5*(p0+p'*yarr*r^2/2) = -1.5*p'*yarr*r^2/4
    if kmaxis > 1 % kmaxis^2 = yavv/yauu
      [~,E] = ellipke(1-1/kmaxis^2);
      E = kmaxis*E;
    else
      [~,E] = ellipke(1-kmaxis^2);
    end
    % (I/C)^2 = jmaxis^2*(pi^2*r^4*(yavv/yauu))/(4*r*Emaxis)^2
    % (I/C)^2 = jmaxis^2*(pi^2*r^2*(yavv/yauu))/(16*Emaxis^2)
    betap(1) = -4/mu0/pi^3*pprime(1)*yauu/jmaxis^2*E^2;
  end

  li = B2V./V./(mu0*It./C).^2; % Volume-averaged Bp^2 / bp2flx
  li(ismaxis) = inf;

  % parallel current density
  jdS = ones(n,1)*fprime.*jdSfp+ones(n,1)*pprime.*jdSpp+ones(n,1)*ffprim.*jdSff;
  jSfp = sum(jdSfp);
  jSpp = sum(jdSpp);
  jSff = sum(jdSff);
  jS = jSfp.*fprime +jSpp.*pprime + jSff.*ffprim;
  jpar = jS./S;
  if psibar(1) == 0 % tiny corrections
    dS(:,1) = dS(:,1)*jpar(1)/jmaxis;
    S(1) = S(1)*jpar(1)/jmaxis;
    jpar(1) = jmaxis;
    jdSpp(:,1) = rmaxis*dS(:,1);
    jdSff(:,1) = dS(:,1)/rmaxis/mu0;
    jdSfp(:,1) = 0; % no fprime, i.e. poloidal current on axis
    jdSfpYr1(:,1) = 0;
    jdSfpYr2(:,1) = 0;
    jdSfpYz1(:,1) = 0;
    jdSfpYz2(:,1) = 0;
    jSpp(1) = rmaxis*S(1);
    jSff(1) = S(1)/rmaxis/mu0;
    jSfp(1) = 0;
  end

  L1t = S./q;

  if include_response

    H = ones(1,16);
    YR = Yr(:);
    YZ = Yz(:);
    dR = DR(:);
    dZ = DZ(:);
    t0 = T0(:);
    t1 = T1(:);
    i2 = reshape([2:n 1]'*ones(1,m)+n*ones(n,1)*[0:m-1],nm,1);
    Iv = Iw(i2,:);

    % Contours
    R0w = sign(dba)*YR./Yg2(:)*H.*W;
    R0b = -sign(dba)*YT(:).*YR./Yg2(:);
    R0a = -sign(dba)*(1-YT(:)).*YR./Yg2(:);
    Z0w = sign(dba)*YZ./Yg2(:)*H.*W;
    Z0b = -sign(dba)*YT(:).*YZ./Yg2(:);
    Z0a = -sign(dba)*(1-YT(:)).*YZ./Yg2(:);
    if isx(ibdef)
      drzdgp = -inv([ybrr ybrz; ybrz ybzz]);
      k = n*(m-1)+ibdef;
      R0w(k,:) = drzdgp(1,1)*wbr+drzdgp(1,2)*wbz;
      Z0w(k,:) = drzdgp(2,1)*wbr+drzdgp(2,2)*wbz;
      R0b(k) = 0;
      Z0b(k) = 0;
      R0a(k) = 0; % These two should be zero but
      Z0a(k) = 0; % are nan or inf if Yg2(k) == 0
    end

    Yrw = Wr + Yrr(:)*H.*R0w + Yrz(:)*H.*Z0w;
    Yrb = Yrr(:).*R0b + Yrz(:).*Z0b;
    Yra = Yrr(:).*R0a + Yrz(:).*Z0a;

    Yzw = Wz + Yrz(:)*H.*R0w + Yzz(:)*H.*Z0w;
    Yzb = Yrz(:).*R0b + Yzz(:).*Z0b;
    Yza = Yrz(:).*R0a + Yzz(:).*Z0a;

    DRw = -R0w;
    DRv = R0w(i2,:);
    DRb = R0b(i2)-R0b;
    DRa = R0a(i2)-R0a;
    DZw = -Z0w;
    DZv = Z0w(i2,:);
    DZb = Z0b(i2)-Z0b;
    DZa = Z0a(i2)-Z0a;
    
    % When F>0 & F<1 then its response will influence responses for T0, T1. This is ignored here.

    % T0 = (Yr.*DR + Yz.*DZ)./(Yr.*DZ - Yz.*DR);
    t1w = YR*H.*DRw+dR*H.*Yrw;
    t2w = YZ*H.*DZw+dZ*H.*Yzw;
    t3w = YR*H.*DZw+dZ*H.*Yrw;
    t4w = YZ*H.*DRw+dR*H.*Yzw;
    T0w = 1./(YR.*dZ-YZ.*dR)*H.*(t1w+t2w) - t0./(YR.*dZ-YZ.*dR)*H.*(t3w-t4w);
    t1v = YR*H.*DRv;
    t2v = YZ*H.*DZv;
    t3v = YR*H.*DZv;
    t4v = YZ*H.*DRv;
    T0v = 1./(YR.*dZ-YZ.*dR)*H.*(t1v+t2v) - t0./(YR.*dZ-YZ.*dR)*H.*(t3v-t4v);
    t1b = YR.*DRb+dR.*Yrb;
    t2b = YZ.*DZb+dZ.*Yzb;
    t3b = YR.*DZb+dZ.*Yrb;
    t4b = YZ.*DRb+dR.*Yzb;
    T0b = 1./(YR.*dZ-YZ.*dR)  .*(t1b+t2b) - t0./(YR.*dZ-YZ.*dR)  .*(t3b-t4b);
    t1a = YR.*DRa+dR.*Yra;
    t2a = YZ.*DZa+dZ.*Yza;
    t3a = YR.*DZa+dZ.*Yra;
    t4a = YZ.*DRa+dR.*Yza;
    T0a = 1./(YR.*dZ-YZ.*dR)  .*(t1a+t2a) - t0./(YR.*dZ-YZ.*dR)  .*(t3a-t4a);

    % T1 = (YR(i2).*DR + YZ(i2).*DZ)./(YR(i2).*DZ - YZ(i2).*DR)
    t1w = YR(i2)*H.*DRw;
    t2w = YZ(i2)*H.*DZw;
    t3w = YR(i2)*H.*DZw;
    t4w = YZ(i2)*H.*DRw;
    T1w = 1./(YR(i2).*dZ - YZ(i2).*dR)*H.*(t1w+t2w) - ...
	t1./(YR(i2).*dZ - YZ(i2).*dR)*H.*(t3w-t4w);
    t1v = YR(i2)*H.*DRv + dR*H.*Yrw(i2,:);
    t2v = YZ(i2)*H.*DZv + dZ*H.*Yzw(i2,:);
    t3v = YR(i2)*H.*DZv + dZ*H.*Yrw(i2,:);
    t4v = YZ(i2)*H.*DRv + dR*H.*Yzw(i2,:);
    T1v = 1./(YR(i2).*dZ - YZ(i2).*dR)*H.*(t1v+t2v) - ...
	t1./(YR(i2).*dZ - YZ(i2).*dR)*H.*(t3v-t4v);
    t1b = YR(i2).*DRb+dR.*Yrb(i2);
    t2b = YZ(i2).*DZb+dZ.*Yzb(i2);
    t3b = YR(i2).*DZb+dZ.*Yrb(i2);
    t4b = YZ(i2).*DRb+dR.*Yzb(i2);
    T1b = 1./(YR(i2).*dZ - YZ(i2).*dR).*(t1b+t2b) - ...
	t1./(YR(i2).*dZ - YZ(i2).*dR).*(t3b-t4b);
    t1a = YR(i2).*DRa+dR.*Yra(i2);
    t2a = YZ(i2).*DZa+dZ.*Yza(i2);
    t3a = YR(i2).*DZa+dZ.*Yra(i2);
    t4a = YZ(i2).*DRa+dR.*Yza(i2);
    T1a = 1./(YR(i2).*dZ - YZ(i2).*dR).*(t1a+t2a) - ...
	t1./(YR(i2).*dZ - YZ(i2).*dR).*(t3a-t4a);

    % Coefficients for interpolation
    R1w = DRw - t0*H.*DZw - dZ*H.*T0w ;
    R1v = DRv - t0*H.*DZv - dZ*H.*T0v;
    R1b = DRb - t0  .*DZb - dZ  .*T0b;
    R1a = DRa - t0  .*DZa - dZ  .*T0a;
    R2w = (2*t0+t1)*H.*DZw + dZ*H.*(2*T0w+T1w);
    R2v = (2*t0+t1)*H.*DZv + dZ*H.*(2*T0v+T1v);
    R2b = (2*t0+t1)  .*DZb + dZ  .*(2*T0b+T1b);
    R2a = (2*t0+t1)  .*DZa + dZ  .*(2*T0a+T1a);
    R3w = -(t0+t1)*H.*DZw - dZ*H.*(T0w+T1w);
    R3v = -(t0+t1)*H.*DZv - dZ*H.*(T0v+T1v);
    R3b = -(t0+t1)  .*DZb - dZ  .*(T0b+T1b);
    R3a = -(t0+t1)  .*DZa - dZ  .*(T0a+T1a);
    Z1w = DZw + t0*H.*DRw + dR*H.*T0w;
    Z1v = DZv + t0*H.*DRv + dR*H.*T0v;
    Z1b = DZb + t0  .*DRb + dR  .*T0b;
    Z1a = DZa + t0  .*DRa + dR  .*T0a;
    Z2w = -(2*t0+t1)*H.*DRw - dR*H.*(2*T0w+T1w);
    Z2v = -(2*t0+t1)*H.*DRv - dR*H.*(2*T0v+T1v);
    Z2b = -(2*t0+t1)  .*DRb - dR  .*(2*T0b+T1b);
    Z2a = -(2*t0+t1)  .*DRa - dR  .*(2*T0a+T1a);
    Z3w = (t0+t1)*H.*DRw + dR*H.*(T0w+T1w);
    Z3v = (t0+t1)*H.*DRv + dR*H.*(T0v+T1v);
    Z3b = (t0+t1)  .*DRb + dR  .*(T0b+T1b);
    Z3a = (t0+t1)  .*DRa + dR  .*(T0a+T1a);

    % dA0 = integral(r(x)*dzdx) from x = 0 to 1
    % dA0 =     R0.*DZ + 1/2*R1.*Z1 + 2/3*R1.*Z2 + 3/4*R1.*Z3 + 1/3*R2.*Z1 + ...
    %       2/4*R2.*Z2 + 3/5*R2.*Z3 + 1/4*R3.*Z1 + 2/5*R3.*Z2 + 3/6*R3.*Z3;
    dA0w = (R0(:)*H.*DZw+dZ*H.*R0w) + ...
       1/2*(R1(:)*H.*Z1w+Z1(:)*H.*R1w) + ...
       2/3*(R1(:)*H.*Z2w+Z2(:)*H.*R1w) + ...
       3/4*(R1(:)*H.*Z3w+Z3(:)*H.*R1w) + ...
       1/3*(R2(:)*H.*Z1w+Z1(:)*H.*R2w) + ...
       2/4*(R2(:)*H.*Z2w+Z2(:)*H.*R2w) + ...
       3/5*(R2(:)*H.*Z3w+Z3(:)*H.*R2w) + ...
       1/4*(R3(:)*H.*Z1w+Z1(:)*H.*R3w) + ...
       2/5*(R3(:)*H.*Z2w+Z2(:)*H.*R3w) + ...
       3/6*(R3(:)*H.*Z3w+Z3(:)*H.*R3w);
    dA0v = (R0(:)*H.*DZv             ) + ... % +dZ*H.*R0v = +0
       1/2*(R1(:)*H.*Z1v+Z1(:)*H.*R1v) + ...
       2/3*(R1(:)*H.*Z2v+Z2(:)*H.*R1v) + ...
       3/4*(R1(:)*H.*Z3v+Z3(:)*H.*R1v) + ...
       1/3*(R2(:)*H.*Z1v+Z1(:)*H.*R2v) + ...
       2/4*(R2(:)*H.*Z2v+Z2(:)*H.*R2v) + ...
       3/5*(R2(:)*H.*Z3v+Z3(:)*H.*R2v) + ...
       1/4*(R3(:)*H.*Z1v+Z1(:)*H.*R3v) + ...
       2/5*(R3(:)*H.*Z2v+Z2(:)*H.*R3v) + ...
       3/6*(R3(:)*H.*Z3v+Z3(:)*H.*R3v);
    dA0b = (R0(:) .* DZb+dZ .* R0b) + ...
       1/2*(R1(:) .* Z1b+Z1(:) .* R1b) + ...
       2/3*(R1(:) .* Z2b+Z2(:) .* R1b) + ...
       3/4*(R1(:) .* Z3b+Z3(:) .* R1b) + ...
       1/3*(R2(:) .* Z1b+Z1(:) .* R2b) + ...
       2/4*(R2(:) .* Z2b+Z2(:) .* R2b) + ...
       3/5*(R2(:) .* Z3b+Z3(:) .* R2b) + ...
       1/4*(R3(:) .* Z1b+Z1(:) .* R3b) + ...
       2/5*(R3(:) .* Z2b+Z2(:) .* R3b) + ...
       3/6*(R3(:) .* Z3b+Z3(:) .* R3b);
    dA0a = (R0(:) .* DZa+dZ .* R0a) + ...
       1/2*(R1(:) .* Z1a+Z1(:) .* R1a) + ...
       2/3*(R1(:) .* Z2a+Z2(:) .* R1a) + ...
       3/4*(R1(:) .* Z3a+Z3(:) .* R1a) + ...
       1/3*(R2(:) .* Z1a+Z1(:) .* R2a) + ...
       2/4*(R2(:) .* Z2a+Z2(:) .* R2a) + ...
       3/5*(R2(:) .* Z3a+Z3(:) .* R2a) + ...
       1/4*(R3(:) .* Z1a+Z1(:) .* R3a) + ...
       2/5*(R3(:) .* Z2a+Z2(:) .* R3a) + ...
       3/6*(R3(:) .* Z3a+Z3(:) .* R3a);

    % dA1 = integral(x*r(x)*dzdx) from x = 0 to 1
    % dA1 = 1/2*R0.*Z1 + 2/3*R0.*Z2 + 3/4*R0.*Z3 + 1/3*R1.*Z1 + 2/4*R1.*Z2 + 3/5*R1.*Z3 + ...
    %       1/4*R2.*Z1 + 2/5*R2.*Z2 + 3/6*R2.*Z3 + 1/5*R3.*Z1 + 2/6*R3.*Z2 + 3/7*R3.*Z3;
    dA1w = 1/2*(R0(:)*H.*Z1w+Z1(:)*H.*R0w) + ...
           2/3*(R0(:)*H.*Z2w+Z2(:)*H.*R0w) + ...
           3/4*(R0(:)*H.*Z3w+Z3(:)*H.*R0w) + ...
           1/3*(R1(:)*H.*Z1w+Z1(:)*H.*R1w) + ...
           2/4*(R1(:)*H.*Z2w+Z2(:)*H.*R1w) + ...
           3/5*(R1(:)*H.*Z3w+Z3(:)*H.*R1w) + ...
           1/4*(R2(:)*H.*Z1w+Z1(:)*H.*R2w) + ...
           2/5*(R2(:)*H.*Z2w+Z2(:)*H.*R2w) + ...
           3/6*(R2(:)*H.*Z3w+Z3(:)*H.*R2w) + ...
           1/5*(R3(:)*H.*Z1w+Z1(:)*H.*R3w) + ...
           2/6*(R3(:)*H.*Z2w+Z2(:)*H.*R3w) + ...
           3/7*(R3(:)*H.*Z3w+Z3(:)*H.*R3w);
    dA1v = 1/2*(R0(:)*H.*Z1v             ) + ... % +Z1(:)*H.*R0v = +0
           2/3*(R0(:)*H.*Z2v             ) + ... % +Z2(:)*H.*R0v = +0
           3/4*(R0(:)*H.*Z3v             ) + ... % +Z3(:)*H.*R0v = +0
           1/3*(R1(:)*H.*Z1v+Z1(:)*H.*R1v) + ...
           2/4*(R1(:)*H.*Z2v+Z2(:)*H.*R1v) + ...
           3/5*(R1(:)*H.*Z3v+Z3(:)*H.*R1v) + ...
           1/4*(R2(:)*H.*Z1v+Z1(:)*H.*R2v) + ...
           2/5*(R2(:)*H.*Z2v+Z2(:)*H.*R2v) + ...
           3/6*(R2(:)*H.*Z3v+Z3(:)*H.*R2v) + ...
           1/5*(R3(:)*H.*Z1v+Z1(:)*H.*R3v) + ...
           2/6*(R3(:)*H.*Z2v+Z2(:)*H.*R3v) + ...
           3/7*(R3(:)*H.*Z3v+Z3(:)*H.*R3v);
    dA1b = 1/2*(R0(:) .* Z1b+Z1(:) .* R0b) + ...
           2/3*(R0(:) .* Z2b+Z2(:) .* R0b) + ...
           3/4*(R0(:) .* Z3b+Z3(:) .* R0b) + ...
           1/3*(R1(:) .* Z1b+Z1(:) .* R1b) + ...
           2/4*(R1(:) .* Z2b+Z2(:) .* R1b) + ...
           3/5*(R1(:) .* Z3b+Z3(:) .* R1b) + ...
           1/4*(R2(:) .* Z1b+Z1(:) .* R2b) + ...
           2/5*(R2(:) .* Z2b+Z2(:) .* R2b) + ...
           3/6*(R2(:) .* Z3b+Z3(:) .* R2b) + ...
           1/5*(R3(:) .* Z1b+Z1(:) .* R3b) + ...
           2/6*(R3(:) .* Z2b+Z2(:) .* R3b) + ...
           3/7*(R3(:) .* Z3b+Z3(:) .* R3b);
    dA1a = 1/2*(R0(:) .* Z1a+Z1(:) .* R0a) + ...
           2/3*(R0(:) .* Z2a+Z2(:) .* R0a) + ...
           3/4*(R0(:) .* Z3a+Z3(:) .* R0a) + ...
           1/3*(R1(:) .* Z1a+Z1(:) .* R1a) + ...
           2/4*(R1(:) .* Z2a+Z2(:) .* R1a) + ...
           3/5*(R1(:) .* Z3a+Z3(:) .* R1a) + ...
           1/4*(R2(:) .* Z1a+Z1(:) .* R2a) + ...
           2/5*(R2(:) .* Z2a+Z2(:) .* R2a) + ...
           3/6*(R2(:) .* Z3a+Z3(:) .* R2a) + ...
           1/5*(R3(:) .* Z1a+Z1(:) .* R3a) + ...
           2/6*(R3(:) .* Z2a+Z2(:) .* R3a) + ...
           3/7*(R3(:) .* Z3a+Z3(:) .* R3a);

    % dA2 = integral(x^2*r(x)*dzdx) from x = 0 to 1
    % dA2 = 1/3*R0.*Z1 + 2/4*R0.*Z2 + 3/5*R0.*Z3 + 1/4*R1.*Z1 + 2/5*R1.*Z2 + 3/6*R1.*Z3 + ...
    %       1/5*R2.*Z1 + 2/6*R2.*Z2 + 3/7*R2.*Z3 + 1/6*R3.*Z1 + 2/7*R3.*Z2 + 3/8*R3.*Z3;
    dA2w = 1/3*(R0(:)*H.*Z1w+Z1(:)*H.*R0w) + ...
           2/4*(R0(:)*H.*Z2w+Z2(:)*H.*R0w) + ...
           3/5*(R0(:)*H.*Z3w+Z3(:)*H.*R0w) + ...
           1/4*(R1(:)*H.*Z1w+Z1(:)*H.*R1w) + ...
           2/5*(R1(:)*H.*Z2w+Z2(:)*H.*R1w) + ...
           3/6*(R1(:)*H.*Z3w+Z3(:)*H.*R1w) + ...
           1/5*(R2(:)*H.*Z1w+Z1(:)*H.*R2w) + ...
           2/6*(R2(:)*H.*Z2w+Z2(:)*H.*R2w) + ...
           3/7*(R2(:)*H.*Z3w+Z3(:)*H.*R2w) + ...
           1/6*(R3(:)*H.*Z1w+Z1(:)*H.*R3w) + ...
           2/7*(R3(:)*H.*Z2w+Z2(:)*H.*R3w) + ...
           3/8*(R3(:)*H.*Z3w+Z3(:)*H.*R3w);
    dA2v = 1/3*(R0(:)*H.*Z1v             ) + ... % +Z1(:)*H.*R0v = +0
           2/4*(R0(:)*H.*Z2v             ) + ... % +Z2(:)*H.*R0v = +0
           3/5*(R0(:)*H.*Z3v             ) + ... % +Z3(:)*H.*R0v = +0
           1/4*(R1(:)*H.*Z1v+Z1(:)*H.*R1v) + ...
           2/5*(R1(:)*H.*Z2v+Z2(:)*H.*R1v) + ...
           3/6*(R1(:)*H.*Z3v+Z3(:)*H.*R1v) + ...
           1/5*(R2(:)*H.*Z1v+Z1(:)*H.*R2v) + ...
           2/6*(R2(:)*H.*Z2v+Z2(:)*H.*R2v) + ...
           3/7*(R2(:)*H.*Z3v+Z3(:)*H.*R2v) + ...
           1/6*(R3(:)*H.*Z1v+Z1(:)*H.*R3v) + ...
           2/7*(R3(:)*H.*Z2v+Z2(:)*H.*R3v) + ...
           3/8*(R3(:)*H.*Z3v+Z3(:)*H.*R3v);
    dA2b = 1/3*(R0(:) .* Z1b+Z1(:) .* R0b) + ...
           2/4*(R0(:) .* Z2b+Z2(:) .* R0b) + ...
           3/5*(R0(:) .* Z3b+Z3(:) .* R0b) + ...
           1/4*(R1(:) .* Z1b+Z1(:) .* R1b) + ...
           2/5*(R1(:) .* Z2b+Z2(:) .* R1b) + ...
           3/6*(R1(:) .* Z3b+Z3(:) .* R1b) + ...
           1/5*(R2(:) .* Z1b+Z1(:) .* R2b) + ...
           2/6*(R2(:) .* Z2b+Z2(:) .* R2b) + ...
           3/7*(R2(:) .* Z3b+Z3(:) .* R2b) + ...
           1/6*(R3(:) .* Z1b+Z1(:) .* R3b) + ...
           2/7*(R3(:) .* Z2b+Z2(:) .* R3b) + ...
           3/8*(R3(:) .* Z3b+Z3(:) .* R3b);
    dA2a = 1/3*(R0(:) .* Z1a+Z1(:) .* R0a) + ...
           2/4*(R0(:) .* Z2a+Z2(:) .* R0a) + ...
           3/5*(R0(:) .* Z3a+Z3(:) .* R0a) + ...
           1/4*(R1(:) .* Z1a+Z1(:) .* R1a) + ...
           2/5*(R1(:) .* Z2a+Z2(:) .* R1a) + ...
           3/6*(R1(:) .* Z3a+Z3(:) .* R1a) + ...
           1/5*(R2(:) .* Z1a+Z1(:) .* R2a) + ...
           2/6*(R2(:) .* Z2a+Z2(:) .* R2a) + ...
           3/7*(R2(:) .* Z3a+Z3(:) .* R2a) + ...
           1/6*(R3(:) .* Z1a+Z1(:) .* R3a) + ...
           2/7*(R3(:) .* Z2a+Z2(:) .* R3a) + ...
           3/8*(R3(:) .* Z3a+Z3(:) .* R3a);

    % dA3 = integral(x^3*r(x)*dzdx) from x = 0 to 1
    % dA3 = 1/4*R0.*Z1 + 2/5*R0.*Z2 + 3/6*R0.*Z3 + 1/5*R1.*Z1 + 2/6*R1.*Z2 + 3/7*R1.*Z3 + ...
    %       1/6*R2.*Z1 + 2/7*R2.*Z2 + 3/8*R2.*Z3 + 1/7*R3.*Z1 + 2/8*R3.*Z2 + 3/9*R3.*Z3;
    dA3w = 1/4*(R0(:)*H.*Z1w+Z1(:)*H.*R0w) + ...
           2/5*(R0(:)*H.*Z2w+Z2(:)*H.*R0w) + ...
           3/6*(R0(:)*H.*Z3w+Z3(:)*H.*R0w) + ...
           1/5*(R1(:)*H.*Z1w+Z1(:)*H.*R1w) + ...
           2/6*(R1(:)*H.*Z2w+Z2(:)*H.*R1w) + ...
           3/7*(R1(:)*H.*Z3w+Z3(:)*H.*R1w) + ...
           1/6*(R2(:)*H.*Z1w+Z1(:)*H.*R2w) + ...
           2/7*(R2(:)*H.*Z2w+Z2(:)*H.*R2w) + ...
           3/8*(R2(:)*H.*Z3w+Z3(:)*H.*R2w) + ...
           1/7*(R3(:)*H.*Z1w+Z1(:)*H.*R3w) + ...
           2/8*(R3(:)*H.*Z2w+Z2(:)*H.*R3w) + ...
           3/9*(R3(:)*H.*Z3w+Z3(:)*H.*R3w);
    dA3v = 1/4*(R0(:)*H.*Z1v             ) + ... % +Z1(:)*H.*R0v = +0
           2/5*(R0(:)*H.*Z2v             ) + ... % +Z2(:)*H.*R0v = +0
           3/6*(R0(:)*H.*Z3v             ) + ... % +Z3(:)*H.*R0v = +0
           1/5*(R1(:)*H.*Z1v+Z1(:)*H.*R1v) + ...
           2/6*(R1(:)*H.*Z2v+Z2(:)*H.*R1v) + ...
           3/7*(R1(:)*H.*Z3v+Z3(:)*H.*R1v) + ...
           1/6*(R2(:)*H.*Z1v+Z1(:)*H.*R2v) + ...
           2/7*(R2(:)*H.*Z2v+Z2(:)*H.*R2v) + ...
           3/8*(R2(:)*H.*Z3v+Z3(:)*H.*R2v) + ...
           1/7*(R3(:)*H.*Z1v+Z1(:)*H.*R3v) + ...
           2/8*(R3(:)*H.*Z2v+Z2(:)*H.*R3v) + ...
           3/9*(R3(:)*H.*Z3v+Z3(:)*H.*R3v);
    dA3b = 1/4*(R0(:) .* Z1b+Z1(:) .* R0b) + ...
           2/5*(R0(:) .* Z2b+Z2(:) .* R0b) + ...
           3/6*(R0(:) .* Z3b+Z3(:) .* R0b) + ...
           1/5*(R1(:) .* Z1b+Z1(:) .* R1b) + ...
           2/6*(R1(:) .* Z2b+Z2(:) .* R1b) + ...
           3/7*(R1(:) .* Z3b+Z3(:) .* R1b) + ...
           1/6*(R2(:) .* Z1b+Z1(:) .* R2b) + ...
           2/7*(R2(:) .* Z2b+Z2(:) .* R2b) + ...
           3/8*(R2(:) .* Z3b+Z3(:) .* R2b) + ...
           1/7*(R3(:) .* Z1b+Z1(:) .* R3b) + ...
           2/8*(R3(:) .* Z2b+Z2(:) .* R3b) + ...
           3/9*(R3(:) .* Z3b+Z3(:) .* R3b);
    dA3a = 1/4*(R0(:) .* Z1a+Z1(:) .* R0a) + ...
           2/5*(R0(:) .* Z2a+Z2(:) .* R0a) + ...
           3/6*(R0(:) .* Z3a+Z3(:) .* R0a) + ...
           1/5*(R1(:) .* Z1a+Z1(:) .* R1a) + ...
           2/6*(R1(:) .* Z2a+Z2(:) .* R1a) + ...
           3/7*(R1(:) .* Z3a+Z3(:) .* R1a) + ...
           1/6*(R2(:) .* Z1a+Z1(:) .* R2a) + ...
           2/7*(R2(:) .* Z2a+Z2(:) .* R2a) + ...
           3/8*(R2(:) .* Z3a+Z3(:) .* R2a) + ...
           1/7*(R3(:) .* Z1a+Z1(:) .* R3a) + ...
           2/8*(R3(:) .* Z2a+Z2(:) .* R3a) + ...
           3/9*(R3(:) .* Z3a+Z3(:) .* R3a);

    % dRA = integral(r(x)^2/2*dzdx) from x = 0 to 1
    % dRA = (R0.*dA0 + R1.*dA1 + R2.*dA2 + R3.*dA3)/2;
    dRAw = 1/2*((R0(:)*H.*dA0w+dA0(:)*H.*R0w) + ...
        	(R1(:)*H.*dA1w+dA1(:)*H.*R1w) + ...
        	(R2(:)*H.*dA2w+dA2(:)*H.*R2w) + ...
        	(R3(:)*H.*dA3w+dA3(:)*H.*R3w));
    dRAv = 1/2*((R0(:)*H.*dA0v              ) + ... % +dA0(:)*H.*R0v = +0
        	(R1(:)*H.*dA1v+dA1(:)*H.*R1v) + ...
        	(R2(:)*H.*dA2v+dA2(:)*H.*R2v) + ...
        	(R3(:)*H.*dA3v+dA3(:)*H.*R3v));
    dRAb = 1/2*((R0(:) .* dA0b+dA0(:) .* R0b) + ...
        	(R1(:) .* dA1b+dA1(:) .* R1b) + ...
        	(R2(:) .* dA2b+dA2(:) .* R2b) + ...
        	(R3(:) .* dA3b+dA3(:) .* R3b));
    dRAa = 1/2*((R0(:) .* dA0a+dA0(:) .* R0a) + ...
        	(R1(:) .* dA1a+dA1(:) .* R1a) + ...
        	(R2(:) .* dA2a+dA2(:) .* R2a) + ...
        	(R3(:) .* dA3a+dA3(:) .* R3a));

    % dZA = integral(r(x)*z(x)*dzdx) from x = 0 to 1
    % dZA = (Z0.*dA0 + Z1.*dA1 + Z2.*dA2 + Z3.*dA3);
    dZAw = (Z0(:)*H.*dA0w+dA0(:)*H.*Z0w) + ...
           (Z1(:)*H.*dA1w+dA1(:)*H.*Z1w) + ...
           (Z2(:)*H.*dA2w+dA2(:)*H.*Z2w) + ...
           (Z3(:)*H.*dA3w+dA3(:)*H.*Z3w);
    dZAv = (Z0(:)*H.*dA0v              ) + ... % +dA0(:)*H.*Z0v = +0
           (Z1(:)*H.*dA1v+dA1(:)*H.*Z1v) + ...
           (Z2(:)*H.*dA2v+dA2(:)*H.*Z2v) + ...
           (Z3(:)*H.*dA3v+dA3(:)*H.*Z3v);
    dZAb = (Z0(:) .* dA0b+dA0(:) .* Z0b) + ...
           (Z1(:) .* dA1b+dA1(:) .* Z1b) + ...
           (Z2(:) .* dA2b+dA2(:) .* Z2b) + ...
           (Z3(:) .* dA3b+dA3(:) .* Z3b);
    dZAa = (Z0(:) .* dA0a+dA0(:) .* Z0a) + ...
           (Z1(:) .* dA1a+dA1(:) .* Z1a) + ...
           (Z2(:) .* dA2a+dA2(:) .* Z2a) + ...
           (Z3(:) .* dA3a+dA3(:) .* Z3a);

    % dAR = integral(log(r(x))*dzdx) from x = 0 to 1, where log(r) is Taylor expanded
    % dAR = (log(R0)-3/2).*DZ+2*dA0./R0-dRA./R0.^2;
    dARw = ((log(R0(:))-3/2)*H.*DZw+dZ./R0(:)*H.*R0w) + ...
           (2./R0(:)*H.*dA0w-2*dA0(:)./R0(:).^2*H.*R0w) - ...
	   (1./R0(:).^2*H.*dRAw-2*dRA(:)./R0(:).^3*H.*R0w);
    dARv = ((log(R0(:))-3/2)*H.*DZv) + ... % +dZ./R0(:)*H.*R0v = +0
           (2./R0(:)*H.*dA0v) - ... % -2*dA0(:)./R0(:).^2*H.*R0v = -0
	   (1./R0(:).^2*H.*dRAv); % -2*dRA(:)./R0(:).^3*H.*R0v = -0
    dARb = ((log(R0(:))-3/2) .* DZb+dZ./R0(:) .* R0b) + ...
           (2./R0(:) .* dA0b-2*dA0(:)./R0(:).^2 .* R0b) - ...
	   (1./R0(:).^2 .* dRAb-2*dRA(:)./R0(:).^3 .* R0b);
    dARa = ((log(R0(:))-3/2) .* DZa+dZ./R0(:) .* R0a) + ...
           (2./R0(:) .* dA0a-2*dA0(:)./R0(:).^2 .* R0a) - ...
	   (1./R0(:).^2 .* dRAa-2*dRA(:)./R0(:).^3 .* R0a);


    % V,A,L
    Vp = zeros(m,ngg);
    Ap = zeros(m,ngg);
    Lp = zeros(m,ngg);
    k = 0;
    for j = 1:m
      for i = 1:n
	k = k+1;
	if j > 1 | ~ismaxis(1)
          Vp(j,Iw(k,:)) = Vp(j,Iw(k,:)) + 2*pi*dRAw(k,:);
          Vp(j,Iv(k,:)) = Vp(j,Iv(k,:)) + 2*pi*dRAv(k,:);
          Vp(j,iib    ) = Vp(j,iib    ) + 2*pi*dRAb(k)*wb;
          Vp(j,iia    ) = Vp(j,iia    ) + 2*pi*dRAa(k)*wa;
	  Ap(j,Iw(k,:)) = Ap(j,Iw(k,:)) + dA0w(k,:);
	  Ap(j,Iv(k,:)) = Ap(j,Iv(k,:)) + dA0v(k,:);
          Ap(j,iib    ) = Ap(j,iib    ) + dA0b(k)*wb;
          Ap(j,iia    ) = Ap(j,iia    ) + dA0a(k)*wa;
	  Lp(j,Iw(k,:)) = Lp(j,Iw(k,:)) + dARw(k,:);
	  Lp(j,Iv(k,:)) = Lp(j,Iv(k,:)) + dARv(k,:);
          Lp(j,iib    ) = Lp(j,iib    ) + dARb(k)*wb;
          Lp(j,iia    ) = Lp(j,iia    ) + dARa(k)*wa;
	end
      end
    end

    % psit
    % dpsit = [fpol(1)*L(1), (fpol(1:m-1)+fpol(2:m)).*diff(L)/2];
    % psit = cumsum(dpsit);
    psitfp = zeros(m,m);
    psitfp(1,1) = L(1);
    for j = 2:m
      psitfp(j,j-1:j) = (L(j)-L(j-1))/2;
    end
    psitfp = cumsum(psitfp);

    % psis
    psisps = zeros(m,ngg);
    % Shift dq intervals back a half step to center on points
    dq2 = (dq+dq([n 1:n-1],:))/2;
    df = dq2./(ones(n,1)*q);
    k = 0;
    for j = 1:m
      for i = 1:n
	k = k+1;
	psisps(j,Iw(k,:)) = psisps(j,Iw(k,:)) + df(k)*W(k,:);
      end
    end
    if ismaxis(1)
      psisps(1,:) = 0;
      psisps(1,iia) = wa;
    end
    psisfp = psitfp./(q'*ones(1,m));

    % Response of jpar to psizr for stationary helical paths
    % Yr, Yz affects jdSdp along segments on both sides of the point
    jYr = ones(n,1)*[fprime./S].*[jdSfpYr1 + jdSfpYr2([n 1:n-1],:)];
    jYz = ones(n,1)*[fprime./S].*[jdSfpYz1 + jdSfpYz2([n 1:n-1],:)];  
    jW = jYr(:)*H.*Wr + jYz(:)*H.*Wz;
    jFps = zeros(m,ngg); % (j)par, (F)ixed-paths, (ps)izr 
    k = 0;
    for j = 1:m
      for i = 1:n
	k = k+1;
	jFps(j,Iw(k,:)) = jFps(j,Iw(k,:)) + jW(k,:);
      end
    end
    % jSfp'./S'.*fprime' = jFps*psizr(:) = fprime contribution to jpar

  end

else
  % No plasma
  R = nan(n+1,m);
  Z = nan(n+1,m);
  R0 = nan(n,m);
  Z0 = nan(n,m);
  R1 = nan(n,m);
  Z1 = nan(n,m);
  R2 = nan(n,m);
  Z2 = nan(n,m);
  R3 = nan(n,m);
  Z3 = nan(n,m);
  Yr = nan(n,m);
  Yz = nan(n,m);
  dC = nan(n,m);
  dS = nan(n,m);
  jdSfpYr1 = nan(n,m);
  jdSfpYr2 = nan(n,m);
  jdSfpYz1 = nan(n,m);
  jdSfpYz2 = nan(n,m);
  jdSfp = nan(n,m);
  jdSpp = nan(n,m);
  jdSff = nan(n,m);
  jdS = nan(n,m);
  dV = nan(n,m);
  R4 = nan(4,m);
  Z4 = nan(4,m);
  rhot = psibar;
  V = zeros(1,m);
  A = zeros(1,m);
  L = zeros(1,m);
  C = zeros(1,m);
  S = zeros(1,m);
  L1t = zeros(1,m);
  q = nan(1,m);
  ft = zeros(1,m);
  rsurf = nan(1,m);
  zsurf = nan(1,m);
  aminor = zeros(1,m);
  bminor = zeros(1,m);
  elong = bminor./aminor;
  tril = (rsurf-R4(4,:))./aminor;
  triu = (rsurf-R4(2,:))./aminor;
  Wth = zeros(1,m);
  It = zeros(1,m);
  betap = zeros(1,m);
  li = ones(1,m);
  psip = nan(1,m);
  psit = nan(1,m);
  psis = nan(1,m);
  jSfp = nan(1,m);
  jSpp = nan(1,m);
  jSff = nan(1,m);
  jS = nan(1,m);
  jpar = nan(1,m);
  Babs = nan(1,m);
  g = nan(1,m);
  if include_response
    W = zeros(nm,16);
    Iw = ones(nm,16);
    Iv = ones(nm,16);
    R0w = nan(nm,16);
    R0b = nan(nm,1);
    R0a = nan(nm,1);
    R1w = nan(nm,16);
    R1v = nan(nm,16);
    R1b = nan(nm,1);
    R1a = nan(nm,1);
    R2w = nan(nm,16);
    R2v = nan(nm,16);
    R2b = nan(nm,1);
    R2a = nan(nm,1);
    R3w = nan(nm,16);
    R3v = nan(nm,16);
    R3b = nan(nm,1);
    R3a = nan(nm,1);
    Z0w = nan(nm,16);
    Z0b = nan(nm,1);
    Z0a = nan(nm,1);
    Z1w = nan(nm,16);
    Z1v = nan(nm,16);
    Z1b = nan(nm,1);
    Z1a = nan(nm,1);
    Z2w = nan(nm,16);
    Z2v = nan(nm,16);
    Z2b = nan(nm,1);
    Z2a = nan(nm,1);
    Z3w = nan(nm,16);
    Z3v = nan(nm,16);
    Z3b = nan(nm,1);
    Z3a = nan(nm,1);
    Yrw = nan(nm,16);
    Yrb = nan(nm,1);
    Yra = nan(nm,1);
    Yzw = nan(nm,16);
    Yzb = nan(nm,1);
    Yza = nan(nm,1);
    Vp = zeros(m,ngg);
    Ap = zeros(m,ngg);
    Lp = zeros(m,ngg);
    psitfp = zeros(m,m);
    psisfp = zeros(m,m);
    psisps = zeros(m,ngg);
    jFps = zeros(m,ngg);
  end
end

if include_response
  % Archive response
  r.W = W;
  d.W = 'psi for contour points = sum(W.*psizr(Iw),2)';
  r.Iw = Iw;
  d.Iw = 'Used with response objects *w';
  r.Iv = Iv;
  d.Iv = 'Used with response objects *v';

  r.R0w = R0w;
  d.R0w = 'dR0 = sum(R0w.*dpsizr(Iw),2)';
  r.R0b = R0b;
  d.R0b = 'dR0 = R0b*dpsibry';
  r.R0a = R0a;
  d.R0a = 'dR0 = R0a*dpsimag';
  r.R1w = R1w;
  d.R1w = 'dR1 = sum(R1w.*dpsizr(Iw),2)';
  r.R1v = R1v;
  d.R1v = 'dR1 = sum(R1v.*dpsizr(Iv),2)';
  r.R1b = R1b;
  d.R1b = 'dR1 = R1b*dpsibry';
  r.R1a = R1a;
  d.R1a = 'dR1 = R1a*dpsimag';
  r.R2w = R2w;
  d.R2w = 'dR2 = sum(R2w.*dpsizr(Iw),2)';
  r.R2v = R2v;
  d.R2v = 'dR2 = sum(R2v.*dpsizr(Iv),2)';
  r.R2b = R2b;
  d.R2b = 'dR2 = R2b*dpsibry';
  r.R2a = R2a;
  d.R2a = 'dR2 = R2a*dpsimag';
  r.R3w = R3w;
  d.R3w = 'dR3 = sum(R3w.*dpsizr(Iw),2)';
  r.R3v = R3v;
  d.R3v = 'dR3 = sum(R3v.*dpsizr(Iv),2)';
  r.R3b = R3b;
  d.R3b = 'dR3 = R3b*dpsibry';
  r.R3a = R3a;
  d.R3a = 'dR3 = R3a*dpsimag';
  
  r.Z0w = Z0w;
  d.Z0w = 'dZ0 = sum(Z0w.*dpsizr(Iw),2)';
  r.Z0b = Z0b;
  d.Z0b = 'dZ0 = Z0b*dpsibry';
  r.Z0a = Z0a;
  d.Z0a = 'dZ0 = Z0a*dpsimag';
  r.Z1w = Z1w;
  d.Z1w = 'dZ1 = sum(Z1w.*dpsizr(Iw),2)';
  r.Z1v = Z1v;
  d.Z1v = 'dZ1 = sum(Z1v.*dpsizr(Iv),2)';
  r.Z1b = Z1b;
  d.Z1b = 'dZ1 = Z1b*dpsibry';
  r.Z1a = Z1a;
  d.Z1a = 'dZ1 = Z1a*dpsimag';
  r.Z2w = Z2w;
  d.Z2w = 'dZ2 = sum(Z2w.*dpsizr(Iw),2)';
  r.Z2v = Z2v;
  d.Z2v = 'dZ2 = sum(Z2v.*dpsizr(Iv),2)';
  r.Z2b = Z2b;
  d.Z2b = 'dZ2 = Z2b*dpsibry';
  r.Z2a = Z2a;
  d.Z2a = 'dZ2 = Z2a*dpsimag';
  r.Z3w = Z3w;
  d.Z3w = 'dZ3 = sum(Z3w.*dpsizr(Iw),2)';
  r.Z3v = Z3v;
  d.Z3v = 'dZ3 = sum(Z3v.*dpsizr(Iv),2)';
  r.Z3b = Z3b;
  d.Z3b = 'dZ3 = Z3b*dpsibry';
  r.Z3a = Z3a;
  d.Z3a = 'dZ3 = Z3a*dpsimag';
    
  r.Yrw = Yrw;
  d.Yrw = 'dYr = sum(Yrw.*dpsizr(Iw),2)';
  r.Yrb = Yrb;
  d.Yrb = 'dYr = Yrb*dpsibry';
  r.Yra = Yra;
  d.Yra = 'dYr = Yrb*dpsimag';
    
  r.Yzw = Yzw;
  d.Yzw = 'dYz = sum(Yzw.*dpsizr(Iw),2)';
  r.Yzb = Yzb;
  d.Yzb = 'dYz = Yzb*dpsibry';
  r.Yza = Yza;
  d.Yza = 'dYz = Yzb*dpsimag';

  r.Vp = Vp;
  d.Vp = 'dV = Vp*dpsizr(:)';
  r.Ap = Ap;
  d.Ap = 'dA = Ap*dpsizr(:)';
  r.Lp = Lp;
  d.Lp = 'dL = Lp*dpsizr(:)';

  r.psitfp = psitfp;
  d.psitfp = 'dpsit = psitfp*dfpol(:)';

  r.psisfp = psisfp;
  d.psisfp = 'dpsis = psisfp*dfpol(:)';

  r.psisps = psisps;
  d.psisps = 'How psis for stationary helical paths changes when psizr changes';

  r.jFps = jFps;
  d.jFps = 'd(jpar)/d(psizr) due to Yr, Yz changes along stationary paths';

  r.info = d;
end

ds.R = 'Radius of contour points + first points again';
p.R = R;

ds.Z = 'Height of contour points + first points again';
p.Z = Z;

ds.R0 = 'Radius of contour points';
p.R0 = R0;

ds.R1 = 'To find contour between points in R0:';
p.R1 = R1;

ds.R2 = 'R = R0 + R1*x^1 + R2*x^2 + R3*x^3;';
p.R2 = R2;

ds.R3 = 'for x in range 0 to 1';
p.R3 = R3;

ds.Z0 = 'Height of contour points';
p.Z0 = Z0;

ds.Z1 = 'To find contour between points in Z0:';
p.Z1 = Z1;

ds.Z2 = 'Z = Z0 + Z1*x^1 + Z2*x^2 + Z3*x^3;';
p.Z2 = Z2;

ds.Z3 = 'for x in range 0 to 1';
p.Z3 = Z3;

ds.Yr = 'd(psi)/dR at contour points';
p.Yr = Yr;

ds.Yz = 'd(psi)/dZ at contour points';
p.Yz = Yz;

ds.dC = 'Distances between points along contour';
p.dC = dC;

ds.dS = 'Field line distances between contour points';
p.dS = dS;

ds.jdSfpYr1 = 'jdSfp = jdSfpYr1.*Yr+jdSfpYr2.*Yr([2:n 1],:)+jdSfpYz1.*Yz+jdSfpYz2.*Yz([2:n 1],:)';
p.jdSfpYr1 = jdSfpYr1;

ds.jdSfpYr2 = 'jdSfp = jdSfpYr1.*Yr+jdSfpYr2.*Yr([2:n 1],:)+jdSfpYz1.*Yz+jdSfpYz2.*Yz([2:n 1],:)';
p.jdSfpYr2 = jdSfpYr2;

ds.jdSfpYz1 = 'jdSfp = jdSfpYr1.*Yr+jdSfpYr2.*Yr([2:n 1],:)+jdSfpYz1.*Yz+jdSfpYz2.*Yz([2:n 1],:)';
p.jdSfpYz1 = jdSfpYz1;

ds.jdSfpYz2 = 'jdSfp = jdSfpYr1.*Yr+jdSfpYr2.*Yr([2:n 1],:)+jdSfpYz1.*Yz+jdSfpYz2.*Yz([2:n 1],:)';
p.jdSfpYz2 = jdSfpYz2;

ds.jdSfp = 'sum(jdSfp).*fprime = term in the integral(j*dS)';
p.jdSfp = jdSfp;

ds.jdSpp = 'sum(jdSpp).*pprime = term in the integral(j*dS)';
p.jdSpp = jdSpp;

ds.jdSff = 'sum(jdSff).*ffprim = term in the integral(j*dS)';
p.jdSff = jdSff;

ds.jdS = 'Integral(j*dS) along field between contour points';
p.jdS = jdS;

ds.dV = 'Volume between contour points, V = cumsum(sum(dV))';
p.dV = dV;

ds.R4 = 'Radius at outer, upper, inner, lower extremes';
p.R4 = R4;

ds.Z4 = 'Height at outer, upper, inner, lower extremes';
p.Z4 = Z4;

ds.pres = 'pressure profile';
p.pres = pres;

ds.fpol = 'fpol profile';
p.fpol = fpol;

ds.pprime = 'pprime = 2*pi*d(pres)/d(psi)';
p.pprime = pprime;

ds.ffprim = 'ffprim = pi*d(fpol^2)/d(psi)';
p.ffprim = ffprim;

ds.fprime = 'fprime = 2*pi*d(fpol)/d(psi)';
p.fprime = fprime;

ds.psibar = 'normalized poloidal flux at contours';
p.psibar = psibar;

ds.rhot = 'normalized radius = sqrt(psit/psit(end))';
p.rhot = rhot;

ds.V = 'Volume within flux surfaces';
p.V = V;

ds.A = 'Area within contours';
p.A = A;

ds.L = 'Integral of 1/R within contours';
p.L = L;

ds.C = 'Contour lengths';
p.C = C;

ds.S = 'Field line lengths through 1 poloidal turn';
p.S = S;

ds.Babs = 'Flux-surface-averaged B field magnitude, <B>';
p.Babs = Babs;

ds.L1t = 'Average field line lengths through 1 toroidal turn = S./q';
p.L1t = L1t;

ds.rsurf = 'R for geometric center of flux surfaces';
p.rsurf = rsurf;

ds.zsurf = 'Z for geometric center of flux surfaces';
p.zsurf = zsurf;

ds.aminor = 'Half the radial width of flux surfaces';
p.aminor = aminor;

ds.bminor = 'Half the vertical height of flux surfaces';
p.bminor = bminor;

ds.elong = 'Elongation profile';
p.elong = elong;

ds.tril = 'Lower triangularity';
p.tril = tril;

ds.triu = 'Upper triangularity';
p.triu = triu;

ds.q = 'Safety factor';
p.q = q;

ds.Wth = 'Thermal energy within contours';
p.Wth = Wth;

ds.It = 'Toroidal current within contours';
p.It = It;

ds.betap = 'Poloidal beta profile';
p.betap = betap;

ds.li = 'Normalized inductance profile';
p.li = li;

ds.psip = 'Poloidal flux';
p.psip = psip;

ds.psit = 'Toroidal flux';
p.psit = psit;

ds.psis = 'Integral(psi*dS)./S = psip+psit./q';
p.psis = psis;

ds.jSfp = '=sum(jdSfp)';
p.jSfp = jSfp;

ds.jSpp = '=sum(jdSpp)';
p.jSpp = jSpp;

ds.jSff = '=sum(jdSff)';
p.jSff = jSff;

ds.jS = 'Integral(j*dS) = jSfp.*fprime+jSpp.*pprime+jSff.*ffprim';
p.jS = jS;

ds.jpar = ['<j' 183 'B>/<B>'];
p.jpar = jpar;

ds.ft = 'Effective trapped particle fraction';
p.ft = ft;

ds.Vres = 'Resistive loop voltage = eta1D.*jS./q or sum(eta2D.*jdS)./q';

ds.g = 'Flux-surface average of X, <X> = sum(g.*X(R0,Z0))';
p.g = g;

if include_response
ds.resp = 'Responses to psizr, psibry, psimag';
p.resp = r;
end

p.info = ds;

