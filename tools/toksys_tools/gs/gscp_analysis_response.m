%  USAGE:   gs_cp_analysis_response
%
%  PURPOSE: Analyze a plasma described by circle model in same manner as
%           gs_eq_analysis does for grid model
%           Also calculate response immediately for each quantity
%
%  INPUTS: r0, z0, a0 = geometric center and minor radius
%          All r, z that are derived from r0,z0,a0 including rbdef, zbdef
%          psih = total flux at rh, z0
%          sp, sf = parameters for pres and fpol
%          ic, iv = coil and vessel currents
%
%  OUTPUTS: Many equilibrium quantities and how they respond to p,
%           where p = [psih r0 z0 a0]

%	
%  METHOD: 
	
%  NOTES:  
	
%
%  WRITTEN BY:  Anders Welander ON 2015-02-28
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

response_count = response_count+1;
if isfield(index_in_y,'response_count')
  lae.y(index_in_y.response_count) = response_count;
end

shape = 'LIM';

% Some extra numbers for the boundary
rmin = r0 - a0;
rmax = r0 + a0;
zmin = z0 - a0;
zmax = z0 + a0;

% Find boundary-defining point
xl =  -((rl-r0).*drl+(zl-z0).*dzl)./dl.^2';
xl(xl < 0) = 0;
xl(xl > 1) = 1;
[dum, i] = min((rl+xl.*drl-r0).^2+(zl+xl.*dzl-z0).^2);
a0 = sqrt(dum);
rbdef = rl(i)+xl(i)*drl(i);
zbdef = zl(i)+xl(i)*dzl(i);
vlim = [r0-rbdef z0-zbdef]/sqrt((r0-rbdef)^2+(z0-zbdef)^2);

% The moving coordinates where flux is measured
rh = linspace(rmin,rmax,nrh);
drh = (rmax-rmin)/(nrh-1);

% Applied flux on grid
psizr_app = reshape(mpc*ic+mpv*iv,nz,nr);

% Applied flux at rh
kh = floor((rh-rg(1))/dr)*nz + floor((z0-zg(1))/dz) + 1;
xh = (rh - rgg(kh))/dr;
wr = [ones(nrh,1) xh(:) xh(:).^2 xh(:).^3]*mx;
wrr = [zeros(nrh,1) ones(nrh,1) 2*xh(:) 3*xh(:).^2]*mx/dr;
xh = (z0-zgg(kh))/dz;
wz = [ones(nrh,1) xh(:) xh(:).^2 xh(:).^3]*mx;
wzz = [zeros(nrh,1) ones(nrh,1) 2*xh(:) 3*xh(:).^2]*mx/dz;
ih = [kh(:)-(nz+1)   kh(:)-nz     kh(:)-(nz-1)   kh(:)-(nz-2) ...
      kh(:)-1        kh(:)        kh(:)+1        kh(:)+2      ...
      kh(:)+(nz-1)   kh(:)+nz     kh(:)+(nz+1)   kh(:)+(nz+2) ...
      kh(:)+(2*nz-1) kh(:)+(2*nz) kh(:)+(2*nz+1) kh(:)+(2*nz+2)];
wh = ...
  [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
   wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
   wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
   wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
wh_r = ...
  [wz(:,1).*wrr(:,1) wz(:,2).*wrr(:,1) wz(:,3).*wrr(:,1) wz(:,4).*wrr(:,1) ...
   wz(:,1).*wrr(:,2) wz(:,2).*wrr(:,2) wz(:,3).*wrr(:,2) wz(:,4).*wrr(:,2) ...
   wz(:,1).*wrr(:,3) wz(:,2).*wrr(:,3) wz(:,3).*wrr(:,3) wz(:,4).*wrr(:,3) ...
   wz(:,1).*wrr(:,4) wz(:,2).*wrr(:,4) wz(:,3).*wrr(:,4) wz(:,4).*wrr(:,4)];
wh_z = ...
  [wzz(:,1).*wr(:,1) wzz(:,2).*wr(:,1) wzz(:,3).*wr(:,1) wzz(:,4).*wr(:,1) ...
   wzz(:,1).*wr(:,2) wzz(:,2).*wr(:,2) wzz(:,3).*wr(:,2) wzz(:,4).*wr(:,2) ...
   wzz(:,1).*wr(:,3) wzz(:,2).*wr(:,3) wzz(:,3).*wr(:,3) wzz(:,4).*wr(:,3) ...
   wzz(:,1).*wr(:,4) wzz(:,2).*wr(:,4) wzz(:,3).*wr(:,4) wzz(:,4).*wr(:,4)];
psihapp = sum(wh'.*psizr_app(ih)');
psihapp_r = sum(wh_r'.*psizr_app(ih)');
psihapp_z = sum(wh_z'.*psizr_app(ih)');
kr0 = min(nr-3,max(1,floor((r0-rg(1))/dr)));
kz1 = min(nz-3,max(1,floor((zmin-zg(1))/dz)))+1;
k = kr0*nz+kz1;
tr = (r0-rgg(k))/dr;
tz = (zmin-zgg(k))/dz;
wbot = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
wbotr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
wbotz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
ibot = k+neighbors';
psibotapp = wbot*psizr_app(ibot);
psibotapp_r = wbotr*psizr_app(ibot);
psibotapp_z = wbotz*psizr_app(ibot);
kz1 = min(nz-3,max(1,floor((zmax-zg(1))/dz)))+1;
k = kr0*nz+kz1;
tz = (zmax-zgg(k))/dz;
wtop = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
wtopr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
wtopz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
itop = k+neighbors';
psitopapp = wtop*psizr_app(itop);
psitopapp_r = wtopr*psizr_app(itop);
psitopapp_z = wtopz*psizr_app(itop);
dpsihappdp(:,nrh+1) = psihapp_r;
dpsihappdp(:,nrh+2) = psihapp_z;
dpsihappdp(:,nrh+3) = (linspace(0,2,nrh)-1).*psihapp_r;
dpsibotappdp(nrh+1) = psibotapp_r;
dpsibotappdp(nrh+2) = psibotapp_z;
dpsibotappdp(nrh+3) = -psibotapp_z;
dpsitopappdp(nrh+1) = psitopapp_r;
dpsitopappdp(nrh+2) = psitopapp_z;
dpsitopappdp(nrh+3) = +psitopapp_z;

% Find magnetic axis
[~, ia] = max(abs(psih-psih(1)));
ia = min(nrh-1,max(2,ia));
psib = (psih(ia-1)-2*psih(ia)+psih(ia+1));
psid = (psih(ia+1)-psih(ia-1))/2;
xa = -psid/psib;
rmaxis = rh(ia)+xa*drh;
zmaxis = z0;
psimag = psih(ia)-psid^2/psib/2;
psibry = psih(1);

% How axis responds to flux changes and moving coordinate system
drmaxisdp = zeros(1,np);
drmaxisdp(ia-1) = -(1/psib+1/2/psid)*xa*drh;
drmaxisdp(ia  ) = +(2/psib         )*xa*drh;
drmaxisdp(ia+1) = -(1/psib-1/2/psid)*xa*drh;
drmaxisdp(nrh+1) = 1;
drmaxisdp(nrh+3) = (rmaxis-r0)/a0;
dpsimagdp = zeros(1,np);
dpsimagdp(ia-1) = +psid/psib/2 + psid^2/psib^2/2;
dpsimagdp(ia  ) = 1 - psid^2/psib^2;
dpsimagdp(ia+1) = -psid/psib/2 + psid^2/psib^2/2;

% Response of "flux compression", i.e. psimag-psibry
dfcdp = (dpsimagdp-dpsibrydp)/(psibry-psimag);

% Calculate normalized flux
psihbari = (psih-psimag)/(psih( 1 )-psimag);
psihbaro = (psih-psimag)/(psih(end)-psimag);

% Coefficients for the nkn third degree polynomials for pres
p0 = c0*sp;
p1 = c1*sp;
p2 = c2*sp;
p3 = c3*sp;

% Coefficients for the nkn 3rd degree polys for fpol^2/2-fpol(nr)^2/2
f0 = c0*sf;
f1 = c1*sf;
f2 = c2*sf;
f3 = c3*sf;

pres = p0(iknotg) + ...
          p1(iknotg).*psibar + ...
          p2(iknotg).*psibar.^2 + ...
          p3(iknotg).*psibar.^3;
pres(nr) = 0;

dpresds = [c0(iknotg,:) + ...
              c1(iknotg,:).*(psibar*ones(1,nkn+2)) + ...
              c2(iknotg,:).*(psibar.^2*ones(1,nkn+2)) + ...
              c3(iknotg,:).*(psibar.^3*ones(1,nkn+2)), ...
	      zeros(nr,nkn+2)];

pprime = (p1(iknotg) + ...
             p2(iknotg)*2.*psibar + ...
	     p3(iknotg)*3.*psibar.^2)*twopi/(psibry-psimag);

dpprimeds = [...
  c1(iknotg,:) + ...
  c2(iknotg,:).*(2*psibar*ones(1,nkn+2)) + ...
  c3(iknotg,:).*(3*psibar.^2*ones(1,nkn+2)), ...
  zeros(nr,nkn+2)]*twopi/(psibry-psimag);

dpprimedp = pprime*dfcdp;

halffpolsquared = f0(iknotg) + ...
                     f1(iknotg).*psibar + ...
		     f2(iknotg).*psibar.^2+...
                     f3(iknotg).*psibar.^3 + ...
		     rzero^2*bzero^2/2;

halffpolsquared(halffpolsquared < 0) = 0;

dhalffpolsquaredds = [zeros(nr,nkn+2), ...
  c0(iknotg,:) + ...
  c1(iknotg,:).*(psibar*ones(1,nkn+2)) + ...
  c2(iknotg,:).*(psibar.^2*ones(1,nkn+2)) + ...
  c3(iknotg,:).*(psibar.^3*ones(1,nkn+2))];

fpol = sqrt(2*halffpolsquared);

dfpolds = 1./fpol*ones(1,ns).*dhalffpolsquaredds;

ffprim = (f1(iknotg) + ...
             f2(iknotg)*2.*psibar + ...
	     f3(iknotg)*3.*psibar.^2)*...
	     twopi/(psibry-psimag);

dffprimds = [zeros(nr,nkn+2), ...
  c1(iknotg,:) + ...
  c2(iknotg,:).*(2*psibar*ones(1,nkn+2)) + ...
  c3(iknotg,:).*(3*psibar.^2*ones(1,nkn+2))] * ...
  twopi/(psibry-psimag);

dffprimdp = ffprim*dfcdp;

ffp = [(ffprim(1:nr-1)+ffprim(2:nr))/2; 0];
dffpds = [(dffprimds(1:nr-1,:)+dffprimds(2:nr,:))/2; ...
          zeros(1,ns)];
ppr = [(pprime(1:nr-1)+pprime(2:nr))/2; 0];
dpprds = [(dpprimeds(1:nr-1,:)+dpprimeds(2:nr,:))/2; ...
          zeros(1,ns)];

% Contour 
i1 = ia; % Index nearest lower rh contour point inside axis
i2 = ia; % Index nearest lower rh contour point outside axis
rconti(1) = rmaxis;
rconti(nr) = rh(1);
rconto(1) = rmaxis;
rconto(nr) = rh(end);
drcontidp = zeros(nr,np);
drcontodp = zeros(nr,np);
drcontidp(1,:) = drmaxisdp;
drcontodp(1,:) = drmaxisdp;
i3 = ia-1:ia+1;
% Need 3:rd degree to solve properly for contour points
for i = 2:nr-1
  dum = (1-psibar(i))*psimag + psibar(i)*psih(1);
  while i1 > 1 & psihbari(i1-1) < psibar(i)
    i1 = i1-1;
  end
  psib = (psih(i1-1)-2*psih(i1)+psih(i1+1));
  psid = (psih(i1+1)-psih(i1-1))/2;
  dps = dum-psih(i1);
  xa = -psid/psib;
  xb = real(sqrt(xa^2+2*dps/psib));
  dxadp = ([-1 0 1]/2/psid-[1 -2 1]/psib)*xa;
  dxbdp = (2*xa*dxadp-2/psib^2*[1 -2 1]*dps-2/psib*[0 1 0])/xb/2;
  dxbda = (1-psibar(i))/psib/xb;
  dxbdb = psibar(i)/psib/xb;
  if abs(xa+xb) < 1 & rh(i1)+(xa+xb)*drh < rconti(i-1)
    rconti(i) = rh(i1)+(xa+xb)*drh;
    drcontidp(i,i1-1:i1+1) = (dxadp+dxbdp)*drh;
    drcontidp(i,1 ) = drcontidp(i,1 ) + dxbdb*drh;
    drcontidp(i,i3) = drcontidp(i,i3) + dxbda*drh*dpsimagdp(i3);
  elseif abs(xa-xb) < 1 & rh(i1)+(xa-xb)*drh < rconti(i-1)
    rconti(i) = rh(i1)+(xa-xb)*drh;
    drcontidp(i,i1-1:i1+1) = (dxadp-dxbdp)*drh;
    drcontidp(i,1 ) = drcontidp(i,1 ) - dxbdb*drh;
    drcontidp(i,i3) = drcontidp(i,i3) - dxbda*drh*dpsimagdp(i3);
  else
    rconti(i) = rconti(i-1)-drh/2;
  end
  dum = (1-psibar(i))*psimag + psibar(i)*psih(end);
  while i2 < nrh-1 & psihbaro(i2+1) < psibar(i)
    i2 = i2+1;
  end
  psib = (psih(i2-1)-2*psih(i2)+psih(i2+1));
  psid = (psih(i2+1)-psih(i2-1))/2;
  dps = dum-psih(i2);
  xa = -psid/psib;
  xb = real(sqrt(xa^2+2*dps/psib));
  dxadp = ([-1 0 1]/2/psid-[1 -2 1]/psib)*xa;
  dxbdp = (2*xa*dxadp-2/psib^2*[1 -2 1]*dps-2/psib*[0 1 0])/xb/2;
  dxbda = (1-psibar(i))/psib/xb;
  dxbdb = psibar(i)/psib/xb;
  if abs(xa+xb) < 1 & rh(i2)+(xa+xb)*drh > rconto(i-1)
    rconto(i) = rh(i2)+(xa+xb)*drh;
    drcontodp(i,i2-1:i2+1) = (dxadp+dxbdp)*drh;
    drcontodp(i,nrh) = drcontodp(i,nrh) + dxbdb*drh;
    drcontodp(i,i3) = drcontodp(i,i3) + dxbda*drh*dpsimagdp(i3);
  elseif abs(xa-xb) < 1 & rh(i2)+(xa-xb)*drh > rconto(i-1)
    rconto(i) = rh(i2)+(xa-xb)*drh;
    drcontodp(i,i2-1:i2+1) = (dxadp-dxbdp)*drh;
    drcontodp(i,nrh) = drcontodp(i,nrh) - dxbdb*drh;
    drcontodp(i,i3) = drcontodp(i,i3) - dxbda*drh*dpsimagdp(i3);
  else
    rconto(i) = rh(i2)+drh/2;
  end
end
drcontidp(:,nrh+1) = 1;
drcontidp(:,nrh+3) = (rconti'-r0)/a0;
drcontodp(:,nrh+1) = 1;
drcontodp(:,nrh+3) = (rconto'-r0)/a0;

% Major radii of current sources (with different pprime, ffprim)
Rmajor = (rconto+rconti)/2;
dRmajordp = (drcontodp+drcontidp)/2;

% Minor radii of current sources
rminor = (rconto-rconti)/2;
drminordp = (drcontodp-drcontidp)/2;

onep = ones(1,np);

% Integral of dA/R within circles
L = 2*pi*(Rmajor - sqrt(rconti.*rconto));
dLdp = pi*(1 - sqrt(rconti./rconto)')*onep.*drcontodp + ...
          pi*(1 - sqrt(rconto./rconti)')*onep.*drcontidp;
	  
% Area of current sources
A = pi*rminor.^2;
dAdp = 2*pi*rminor'*onep.*drminordp;

% Volumes within circles
V = 2*pi*Rmajor.*A;
dVdp = 2*pi*Rmajor'*onep.*dAdp + 2*pi*A'*onep.*dRmajordp;

% Calculate plasma flux & current profile
psihpla = zeros(1,nrh);
jp = zeros(1,nr);
jf = zeros(1,nr);
djpdp = zeros(nr,np);
djfdp = zeros(nr,np);
djpds = zeros(nr,ns);
djfds = zeros(nr,ns);
psihpla = zeros(1,nrh);
dpsihpladp = zeros(nrh,np);
dpsihplads = zeros(nrh,ns);
ffp = [(ffprim(1:nr-1)+ffprim(2:nr))/2; 0];
dffpds = [(dffprimds(1:nr-1,:)+dffprimds(2:nr,:))/2; zeros(1,ns)];
ppr = [(pprime(1:nr-1)+pprime(2:nr))/2; 0];
dpprds = [(dpprimeds(1:nr-1,:)+dpprimeds(2:nr,:))/2; zeros(1,ns)];
Ip = cumsum([0 diff(V).*ppr(1:nr-1)'/twopi]);
If = cumsum([0 diff(L).*ffp(1:nr-1)'/mu0]);
dIpds = cumsum([zeros(1,ns); diff(V)'*ones(1,ns).*dpprds(1:nr-1,:)/twopi]);
dIfds = cumsum([zeros(1,ns);  diff(L)'*ones(1,ns).*dffpds(1:nr-1,:)/mu0]);
dIpdp = Ip'*dfcdp+cumsum([zeros(1,np); ppr(1:nr-1)*onep.*diff(dVdp)/twopi]);
dIfdp = If'*dfcdp+cumsum([zeros(1,np); ffp(1:nr-1)*onep.*diff(dLdp)/mu0]);
jp = diff(Ip);
djpdp = diff(dIpdp);
djpds = diff(dIpds);
jp = -V(2:nr)'.*diff(ppr)/twopi;
djpdp = jp*dfcdp-diff(ppr)/twopi*onep.*dVdp(2:nr,:);
djpds = -V(2:nr)'*ones(1,ns).*diff(dpprds)/twopi;
jf = diff(If);
djfdp = diff(dIfdp);
djfds = diff(dIfds);
jf = -L(2:nr)'.*diff(ffp)/mu0;
djfdp = jf*dfcdp-diff(ffp)/mu0*onep.*dLdp(2:nr,:);
djfds = -L(2:nr)'*ones(1,ns).*diff(dffpds)/mu0;
g = drh/1000;
for i = 1:nr-1
  fph = jp(i)*mpcircle2point(rh,0,Rmajor(i+1),0,rminor(i+1));
  fphR = (jp(i)*mpcircle2point(rh,0,Rmajor(i+1)+g,0,rminor(i+1)) - fph)'/g;
  fphr = (jp(i)*mpcircle2point(rh,0,Rmajor(i+1),0,rminor(i+1)+g) - fph)'/g;  
  fphrh = (jp(i)*mpcircle2point(rh+g,0,Rmajor(i+1),0,rminor(i+1)) - fph)'/g;  
  ffh = jf(i)*mfcircle2point(rh,0,Rmajor(i+1),0,rminor(i+1));
  ffhR = (jf(i)*mfcircle2point(rh,0,Rmajor(i+1)+g,0,rminor(i+1)) - ffh)'/g;
  ffhr = (jf(i)*mfcircle2point(rh,0,Rmajor(i+1),0,rminor(i+1)+g) - ffh)'/g;
  ffhrh = (jf(i)*mfcircle2point(rh+g,0,Rmajor(i+1),0,rminor(i+1)) - ffh)'/g;
  psihpla = psihpla + fph + ffh;
  dpsihpladp = dpsihpladp + ...
    fphR*dRmajordp(i+1,:) + fphr*drminordp(i+1,:) + ...
    ffhR*dRmajordp(i+1,:) + ffhr*drminordp(i+1,:);
  dpsihpladp(:,nrh+1) = dpsihpladp(:,nrh+1) + fphrh + ffhrh;
  dpsihpladp(:,nrh+3) = dpsihpladp(:,nrh+3) + ...
    (linspace(0,2,nrh)'-1).*(fphrh+ffhrh);
  if jp(i) ~= 0
    dpsihpladp = dpsihpladp + fph(:)/jp(i)*djpdp(i,:);
    dpsihplads = dpsihplads + fph(:)/jp(i)*djpds(i,:);
  end
  if jf(i) ~= 0
    dpsihpladp = dpsihpladp + ffh(:)/jf(i)*djfdp(i,:);
    dpsihplads = dpsihplads + ffh(:)/jf(i)*djfds(i,:);
  end
end
I = Ip+If;
dIdp = dIpdp +dIfdp ;
dIds = dIpds +dIfds ;

psiherr   = psih - psihpla - psihapp;

cpasma = I(end);
dcpasmadp = dIdp(end,:);
dcpasmads = dIds(end,:);
	    
integralRjphidA = ...
  -A(2:nr).*(Rmajor(2:nr).^2+rminor(2:nr).^2/4)*diff(ppr) + ...
  -A(2:nr)*diff(ffp)/mu0;

dintegralRjphidAds = ...
  -A(2:nr).*(Rmajor(2:nr).^2+rminor(2:nr).^2/4)*diff(dpprds) + ...
  -A(2:nr)*diff(dffpds)/mu0;

dintegralRjphidAdp = ...
  -(diff(ffp')/mu0+(Rmajor(2:nr).^2+rminor(2:nr).^2/4).*diff(ppr'))*dAdp(2:nr,:) + ...
  -A(2:nr).*Rmajor(2:nr)*2.*diff(ppr')*dRmajordp(2:nr,:) + ...
  -A(2:nr).*rminor(2:nr)/2.*diff(ppr')*drminordp(2:nr,:) + ...
  integralRjphidA*dfcdp;
  
rcur = integralRjphidA/cpasma;

drcurds = dintegralRjphidAds/cpasma - ...
             rcur/cpasma*dcpasmads;

drcurdp = dintegralRjphidAdp/cpasma - ...
             rcur/cpasma*dcpasmadp;

zcur = z0;
dzcurdp = [zeros(1,nrh) 0 1 0];

psigrid = linspace(psimag,psibry,nr)';

dpsigriddp = linspace(1,0,nr)'*dpsimagdp + linspace(0,1,nr)'*dpsibrydp;

psig = (psigrid(1:nr-1)+psigrid(2:nr))/2;

dpsigdp = (dpsigriddp(1:nr-1,:)+dpsigriddp(2:nr,:))/2;

ipjdA = diff(V)*(ppr(1:nr-1).*psig)/twopi + ...
           diff(L)*(ffp(1:nr-1).*psig)/mu0;

dipjdAds = ...
  diff(V)*(psig*ones(1,ns).*dpprds(1:nr-1,:))/twopi + ...
  diff(L)*(psig*ones(1,ns).*dffpds(1:nr-1,:))/mu0;
	
dipjdAdp = (ppr(1:nr-1).*psig)'/twopi*diff(dVdp) + ...
              diff(V).*ppr(1:nr-1)'/twopi*dpsigdp + ...
	      (ffp(1:nr-1).*psig)'/mu0*diff(dLdp) + ...
              diff(L).*ffp(1:nr-1)'/mu0*dpsigdp + ...
              ipjdA*dfcdp;	
	
psipla = ipjdA/cpasma;

dpsiplads = dipjdAds/cpasma - ...
               psipla/cpasma*dcpasmads;

dpsipladp = dipjdAdp/cpasma - ...
               psipla/cpasma*dcpasmadp;

T = [0 cumsum(diff(L).*(fpol(1:nr-1)+fpol(2:nr))')]/2;

dTdp = [zeros(1,np); ...
  cumsum((fpol(1:nr-1)+fpol(2:nr))/2*onep.*diff(dLdp))];

dTds = [zeros(1,ns); ...
  cumsum(diff(L')*ones(1,ns).*...
  (dfpolds(1:nr-1,:)+dfpolds(2:nr,:)))]/2;

rhot = sqrt(T/T(end));

drhotdp = rhot'/2./T'*onep.*dTdp - ...
             rhot'/2/T(end)*dTdp(end,:);
drhotdp(1,:) = 0;

drhotds = rhot'/2./T'*ones(1,ns).*dTds - ...
             rhot'/2/T(end)*dTds(end,:);
drhotds(1,:) = 0;

% The derivative of T is q, evaluate for psibar values
qpsi = [4*T(2)-3*T(1)-T(3), T(3:nr)-T(1:nr-2), ...
  3*T(nr)+T(nr-2)-4*T(nr-1)]/2/(psimag-psibry)*(nr-1);

dqpsidp = [4*dTdp(2,:)-3*dTdp(1,:)-dTdp(3,:); ...
              dTdp(3:nr,:)-dTdp(1:nr-2,:); ...
	      3*dTdp(nr,:)+dTdp(nr-2,:)-4*dTdp(nr-1,:)]/...
              2/(psimag-psibry)*(nr-1) + ...
	      qpsi'*dfcdp;

dqpsids = [4*dTds(2,:)-3*dTds(1,:)-dTds(3,:); ...
              dTds(3:nr,:)-dTds(1:nr-2,:); ...
	      3*dTds(nr,:)+dTds(nr-2,:)-4*dTds(nr-1,:)]/...
              2/(psimag-psibry)*(nr-1);

torflux = T(nr);

W = [0 cumsum(diff(V).*(pres(1:nr-1)+pres(2:nr))')]/4*3;

dWdp = [zeros(1,np); ...
  cumsum((pres(1:nr-1)+pres(2:nr))/4*3*onep.*diff(dVdp))];

dWds = [zeros(1,ns); ...
  cumsum(diff(V')*ones(1,ns).*...
  (dpresds(1:nr-1,:)+dpresds(2:nr,:)))]/4*3;

% The circumference, Cl
Cl = 2*pi*a0;
dCldp = [zeros(1,nrh+2) 2*pi];

bp2flx = (mu0*cpasma/Cl)^2;
dbp2flxds = 2*bp2flx*(dcpasmads/cpasma);
dbp2flxdp = 2*bp2flx*(dcpasmadp/cpasma - dCldp/Cl);

betap = 4/3*mu0*W(nr)/V(nr)/bp2flx;

dbetapds = betap*(dWds(nr,:)/W(nr) - dbp2flxds/bp2flx);
	      
dbetapdp = betap*(dWdp(nr,:)/W(nr) - ...
                        dVdp(nr,:)/V(nr) - ...
                        dbp2flxdp/bp2flx);
	      
% Solution for circular concentric rings used to calculate Bp2V3
cbp2v =  (rminor(1:nr-1)+rminor(2:nr))./(Rmajor(1:nr-1)+Rmajor(2:nr));

dcbp2vdp = 1./(Rmajor(1:nr-1)+Rmajor(2:nr))'*onep.*...
              (drminordp(1:nr-1,:)+drminordp(2:nr,:)) - ...
	      cbp2v'./(Rmajor(1:nr-1)+Rmajor(2:nr))'*onep.*...
	      (dRmajordp(1:nr-1,:)+dRmajordp(2:nr,:));

Bp2V3 = cbp2v./diff(rminor)*diff(psigrid).^2;

dBp2V3dp = diff(psigrid').^2./diff(rminor)*dcbp2vdp + ...
  2*diff(psigrid').*cbp2v./diff(rminor)*diff(dpsigriddp) - ...
  diff(psigrid').^2.*cbp2v./diff(rminor).^2*diff(drminordp);

li3 = Bp2V3/V(nr)/bp2flx;

dli3ds = -li3/bp2flx*dbp2flxds;

dli3dp = li3*(dBp2V3dp/Bp2V3 - ...
                    dVdp(nr,:)/V(nr) - ...
	            dbp2flxdp/bp2flx);

% Solution for Shafranov-shifted flux surfaces
Bp2V = 0;
dBp2Vdp = zeros(1,np);
for i = 1:nr-1
  t = rminor(i+1)-rminor(i);
  r = rminor(i+1)+rminor(i);
  R = Rmajor(i+1)+Rmajor(i);
  Ro = R+r;            % a.k.a. rconto(i+1)+rconto(i)
  s = sqrt(2*R/Ro-1);  % Between 1 & ~0.8
  g = R*(R/Ro-s) + R*(s-1) + r;
  h = g/r/s/t;
  Bp2V = Bp2V + h;
  
  dsdr = -R/Ro^2/s;
  dsdR = r/Ro^2/sqrt((R-r)/Ro);
  dgdr = -R^2/Ro^2 + 1;
  dgdR = 2*R/Ro-R^2/Ro^2-1;
  dhdr = h/g*dgdr - h/r - h/s*dsdr;
  dhdR = h/g*dgdR - h/s*dsdR;
  dBp2Vdp = dBp2Vdp + ...
    (dhdr-h/t)*drminordp(i+1,:) + ...
    (dhdr+h/t)*drminordp(i  ,:) + ...
    dhdR*(dRmajordp(i+1,:)+dRmajordp(i,:));
end
Bp2V = Bp2V*((psibry-psimag)/(nr-1))^2;
dBp2Vdp = dBp2Vdp*((psibry-psimag)/(nr-1))^2 + ...
  2*Bp2V/(psibry-psimag)*dpsibrydp - ...
  2*Bp2V/(psibry-psimag)*dpsimagdp;

li = Bp2V/V(nr)/bp2flx;

dlids = -li/bp2flx*dbp2flxds;

dlidp = li*(dBp2Vdp/Bp2V - ...
                  dVdp(nr,:)/V(nr) - ...
	          dbp2flxdp/bp2flx);
		  
% The touch point response
drbdefdp = [zeros(1,nrh) 1 0 -vlim(1)];
dzbdefdp = [zeros(1,nrh) 0 1 -vlim(2)];
% Here touch point is where circle is parallel to limiter

% Displacement of touch point normal into the limiter
nbdef = 0;
dnbdefdp = vlim*[drbdefdp; dzbdefdp];

% The variables in equation M are: total flux change at rh points 
% and shifts of r0, z0, a0 for a total of nrh+3
% The equations are: changes in (total - plasma = applied flux), and:
% 1) flux difference psih(end)-psih(1), 
% 2) psitopapp-psibotapp
% 3) displacement into the wall
% The right-hand-side to 1 & 2 can be non-zero when fixing errors

cpM = ...
  [eye(nrh,np)-dpsihappdp-dpsihpladp; ... % tot- pla = app flux
  -1 zeros(1,nrh-2) 1 0 0 0 ; ...   % psih(end) - psih(1)
  dpsitopappdp-dpsibotappdp; ... % psitopapp-psibotapp
  dnbdefdp];                        % = 0 to not go through wall

perr = [psiherr(:);psih(end)-psih(1);psitopapp-psibotapp;0];
er = 1;

cpMi = inv(cpM);

for i = 1:nc
  dpsihappdis(:,i) = sum(wh'.*reshape(mpc(ih,i),nrh,16)');
end
for i = 1:nv
  dpsihappdis(:,nc+i) = sum(wh'.*reshape(mpv(ih,i),nrh,16)');
end
dpsibotappdis = wbot*[mpc(ibot,:) mpv(ibot,:)];
dpsitopappdis = wtop*[mpc(itop,:) mpv(itop,:)];

dpdx = cpMi*[[dpsihappdis dpsihplads; ...
  zeros(1,nx-1); ...
  dpsibotappdis-dpsitopappdis zeros(1,2*nkn+4); ...
  zeros(1,nx-1)] perr(:)];

dcpasmadx = dcpasmadp*dpdx;
dcpasmadx([indsp indsf]) = dcpasmadx([indsp indsf])+dcpasmads;

dlidx = dlidp*dpdx;
dlidx([indsp indsf]) = dlidx([indsp indsf])+dlids;

dbetapdx = dbetapdp*dpdx;
dbetapdx([indsp indsf]) = dbetapdx([indsp indsf])+dbetapds;

dWdx = dWdp*dpdx;
dWdx(:,[indsp indsf]) = dWdx(:,[indsp indsf])+dWds;

Wth = W(nr);
dWthdx = dWdx(nr,:);

dLdx = dLdp*dpdx;
dAdx = dAdp*dpdx;
dVdx = dVdp*dpdx;

Ltot = L(nr);
Atot = A(nr);
Vtot = V(nr);
dLtotdx = dLdx(nr,:);
dAtotdx = dAdx(nr,:);
dVtotdx = dVdx(nr,:);

if 0
plot(rh,psihpla)
title(num2str(time))
drawnow
end

if constraints == 1 % Constrain with user-supplied cpasma, li, betap

% Target specific ip, li, betap using spline vectors sp0, sf0, sg0

% Reduce dimension of full state vector to x3 = [ic;iv;cp;cf;cg;er]
% which maps onto the full vector as [ic;iv;cp*sp0;cf*sf0+cg*sg0;er]

% The constrained state vector is: xc = [ic;iv;cpasma;li;betap;er]

  dcpasmadx3 = [dcpasmadx(indis)     ... 
                dcpasmadx(indsp)*sp0 ...
                dcpasmadx(indsf)*sf0 ...
	        dcpasmadx(indsf)*sg0 ...
	        dcpasmadx(inder)];
	       
  dlidx3 = [dlidx(indis)     ...
            dlidx(indsp)*sp0 ...
            dlidx(indsf)*sf0 ...
	    dlidx(indsf)*sg0 ...
	    dlidx(inder)];
	   
  dbetapdx3 = [dbetapdx(1:nc+nv)   ...
               dbetapdx(indsp)*sp0 ...
               dbetapdx(indsf)*sf0 ...
	       dbetapdx(indsf)*sg0 ...
	       dbetapdx(inder)];

  dx1dx3 = [[eye(nc+nv) zeros(nc+nv,4)]; ...
            dcpasmadx3; ...
	    dlidx3; ...
	    dbetapdx3; ...
	    [zeros(1,nc+nv+3) 1]];

  dx3dx1 = inv(dx1dx3);

  dxdx3 = [[eye(nc+nv) zeros(nc+nv,4)]; ...
	  sp0*[zeros(1,nc+nv) 1 0 0 0]; ...
	  sf0*[zeros(1,nc+nv) 0 1 0 0]+ ...
	  sg0*[zeros(1,nc+nv) 0 0 1 0]; ...
	      [zeros(1,nc+nv+3) 1]   ];

  dxdxc = dxdx3*dx3dx1;
  
elseif constraints == 2 % Constrain with user-supplied cpasma, li, Wth

% Target specific ip, li, Wth using spline vectors sp0, sf0, sg0

% Reduce dimension of full state vector to x3 = [ic;iv;cp;cf;cg;er]
% which maps onto the full vector as [ic;iv;cp*sp0;cf*sf0+cg*sg0;er]

% The constrained state vector is: xc = [ic;iv;cpasma;li;Wth;er]

  dcpasmadx3 = [dcpasmadx(indis)     ... 
                dcpasmadx(indsp)*sp0 ...
                dcpasmadx(indsf)*sf0 ...
	        dcpasmadx(indsf)*sg0 ...
	        dcpasmadx(inder)];
	       
  dlidx3 = [dlidx(indis)     ...
            dlidx(indsp)*sp0 ...
            dlidx(indsf)*sf0 ...
	    dlidx(indsf)*sg0 ...
	    dlidx(inder)];
	   
  dWthdx3 = [dWthdx(1:nc+nv)   ...
               dWthdx(indsp)*sp0 ...
               dWthdx(indsf)*sf0 ...
	       dWthdx(indsf)*sg0 ...
	       dWthdx(inder)];

  dx1dx3 = [[eye(nc+nv) zeros(nc+nv,4)]; ...
            dcpasmadx3; ...
	    dlidx3; ...
	    dWthdx3; ...
	    [zeros(1,nc+nv+3) 1]];

  dx3dx1 = inv(dx1dx3);

  dxdx3 = [[eye(nc+nv) zeros(nc+nv,4)]; ...
	  sp0*[zeros(1,nc+nv) 1 0 0 0]; ...
	  sf0*[zeros(1,nc+nv) 0 1 0 0]+ ...
	  sg0*[zeros(1,nc+nv) 0 0 1 0]; ...
	      [zeros(1,nc+nv+3) 1]   ];

  dxdxc = dxdx3*dx3dx1;
  
end

dpsihpladx = dpsihpladp*dpdx;
dpsihpladx(:,[indsp indsf]) = dpsihpladx(:,[indsp indsf])+dpsihplads;

dpsihappdx = dpsihappdp*dpdx;
dpsihappdx(:,indis) = dpsihappdx(:,indis) + dpsihappdis;

dIdx = dIdp*dpdx;
dIdx(:,[indsp indsf]) = dIdx(:,[indsp indsf])+dIds;


drcurdx = drcurdp*dpdx;
drcurdx([indsp indsf]) = drcurdx([indsp indsf])+drcurds;

dzcurdx = dpdx(nrh+2,:);

drmaxisdx = drmaxisdp*dpdx;
dzmaxisdx = dpdx(nrh+2,:);
dpsimagdx = dpsimagdp*dpdx;

% The response of bdef position
drbdefdx = dpdx(nrh+1,:)-vlim(1)*dpdx(nrh+3,:);
dzbdefdx = dpdx(nrh+2,:)-vlim(2)*dpdx(nrh+3,:);
dpsibrydx = dpdx(1,:);

psihbar = (psih-psibry)/(psimag-psibry);
dpsihbardx = (dpdx(1:nrh,:)-ones(nrh,1)*dpsibrydx)/(psimag-psibry) - ...
  psihbar(:)/(psimag-psibry)*(dpsimagdx-dpsibrydx);

% The response of plasma flux
dipjdAdx = dipjdAdp*dpdx;
dipjdAdx([indsp indsf]) = dipjdAdx([indsp indsf])+dipjdAds;
psipla = ipjdA/cpasma;
dpsipladx = dipjdAdx/cpasma - psipla/cpasma*dcpasmadx;

dCldx = dCldp*dpdx;
dbp2flxdx = 2*bp2flx*(dcpasmadx/cpasma - dCldx/Cl);

nbbbs = nbbbs_max;
thbbbs = angle(rbdef-r0 + 1i*zbdef-1i*z0) + linspace(0,2*pi,nbbbs_max)';
rbbbs = r0 + a0*cos(thbbbs);
zbbbs = z0 + a0*sin(thbbbs);

dr0dx = dpdx(nrh+1,:);
dz0dx = dpdx(nrh+2,:);
da0dx = dpdx(nrh+3,:);
drsurfdx = dpdx(nrh+1,:);
daminordx = dpdx(nrh+3,:);

% For gs_rtplot
drbdx = ones(nbbbs_max,1)*dr0dx + cos(thbbbs(:))*da0dx;
dzbdx = ones(nbbbs_max,1)*dz0dx + sin(thbbbs(:))*da0dx;

% Map the circle onto the grid
% pcurrt will be filled with currents according to
% local current density at the grid point and on top of that
% a uniform current on all grid elements which is designed
% to make a match for cpasma, rcur, zcur between ring and grid models
psibarzr = sqrt((rgg-r0).^2+(zgg-z0).^2)/a0;
iplasma = psibarzr < 1;
dpsibarzrdx = zeros(ngg,nx);
for i = 1:nz
  for j = 1:nr
    if iplasma(i,j)
      for k = 2:nr
	r = sqrt((rg(j)-Rmajor(k))^2+(zg(i)-z0)^2);
	if rminor(k-1) < r & r < rminor(k)
	  g = (rminor(k) - rminor(k-1));
	  h = (   r -    rminor(k-1));
	  psibarzr(i,j) = psibar(k-1) + h/g*(psibar(k)-psibar(k-1));
	  dRdx = dRmajordp(k,:)*dpdx;
	  drdx = -((rg(j)-Rmajor(k))*dRdx+(zg(i)-z0)*dzcurdx)/r;
	  dgdx = drminordp(k,:)*dpdx - drminordp(k-1,:)*dpdx;
	  dhdx = drdx - drminordp(k-1,:)*dpdx;
	  dpsibarzrdx(i+(j-1)*nz,:) = (dhdx/g-h/g^2*dgdx)*(psibar(k)-psibar(k-1));
	  if 0
	  R = Rmajor(k) + h*(Rmajor(k-1)-Rmajor(k));
	  % A more accurate r
	  r = sqrt((rg(j)-R)^2+(zg(i)-z0)^2);
          psibarzr(i,j) = psibar(k) - (rminor(k)-r)/g*(psibar(k)-psibar(k-1));
	  dRdx = dRmajordp(k,:)*dpdx;
	  drdx = ((rg(j)-dRdx)+(zg(i)-dzcurdx))/r;
	  dpsibarzrdx(i+(j-1)*nz,:) = (rminor(k)-r)/g*(psibar(k)-psibar(k-1));
	  end
	end
      end
    end
  end
end
dRmajordx = dRmajordp*dpdx;
drminordx = drminordp*dpdx;

npla = sum(iplasma(:));
% Grab nearest grid cell if plasma is too tiny to cover a grid point
if npla == 0
  j = round((r0-rg(1))/dr+1);
  i = round((z0-zg(1))/dz+1);
  iplasma(i,j) = 1;
  psibarzr(i,j) = 1;
  npla = 1;
end
iknot = ones(nz,nr);
for i = 2:nkn
  iknot(psibarzr > psikn(i)) = i;
end
iknot(~iplasma) = nkn+1;
fpolzr = sign(rzero*bzero)*sqrt(2*(...
         f0(iknot) + ...
         f1(iknot).*psibarzr + ...
	 f2(iknot).*psibarzr.^2 + ...
	 f3(iknot).*psibarzr.^3 + ...
	 rzero^2*bzero^2/2));

preszr = p0(iknot) + ...
         p1(iknot).*psibarzr + ...
	 p2(iknot).*psibarzr.^2 + ...
	 p3(iknot).*psibarzr.^3;

pprimezr = (p1(iknot) + ...
            p2(iknot)*2.*psibarzr + ...
	    p3(iknot)*3.*psibarzr.^2)*twopi/(psibry-psimag);

ffprimzr = (f1(iknot) + ...
            f2(iknot)*2.*psibarzr + ...
	    f3(iknot)*3.*psibarzr.^2)*twopi/(psibry-psimag);

c1p = twopi/(psibry-psimag)*c1;
c2p = twopi/(psibry-psimag)*c2;
c3p = twopi/(psibry-psimag)*c3;
ps2 = 2*psibarzr(:)*ones(1,nkn+2);
ps3 = 3*psibarzr(:).^2*ones(1,nkn+2);
c1pzr = reshape(c1p(iknot,:),ngg,nkn+2);
c2pzr = reshape(c2p(iknot,:),ngg,nkn+2);
c3pzr = reshape(c3p(iknot,:),ngg,nkn+2);

dpprimezrdx = (p2(iknot(:))*2 + p3(iknot(:))*6.*psibarzr(:))*...
  ones(1,nx).*dpsibarzrdx*twopi/(psibry-psimag) - ...
  pprimezr(:)/(psibry-psimag)*(dpsibrydx-dpsimagdx);
dpprimezrdx(:,indsp) = dpprimezrdx(:,indsp) + c1pzr+c2pzr.*ps2+c3pzr.*ps3;

dffprimzrdx = (f2(iknot(:))*2 + f3(iknot(:))*6.*psibarzr(:))*...
  ones(1,nx).*dpsibarzrdx*twopi/(psibry-psimag) - ...
  ffprimzr(:)/(psibry-psimag)*(dpsibrydx-dpsimagdx);
dffprimzrdx(:,indsf) = dffprimzrdx(:,indsf) + c1pzr+c2pzr.*ps2+c3pzr.*ps3;

pcurrt = RA.*pprimezr+AR.*ffprimzr/mu0; % Toroidal current
dpcurrtdx = RA(:)*ones(1,nx).*dpprimezrdx+AR(:)*ones(1,nx).*dffprimzrdx/mu0;

sumjphi = sum(pcurrt(:));
sumrggjphi = sum(rgg(:).*pcurrt(:));
sumzggjphi = sum(zgg(:).*pcurrt(:));
dsumjphidx = sum(dpcurrtdx);

% Make the sum of pcurrt match cpasma
pcurrt(iplasma) = pcurrt(iplasma) + (cpasma-sumjphi)/npla;
dpcurrtdx(iplasma,:) = dpcurrtdx(iplasma,:) + ones(npla,1)*(dcpasmadx-dsumjphidx)/npla;

% Make the current centroid for pcurrt match rcur
pc2dr = ([zeros(nz,1) pcurrt(:,1:nr-1)]-[pcurrt(:,2:nr) zeros(nz,1)])/(2*dr);
rcur2 = rgg(:)'*pcurrt(:)/cpasma;
drcur2dx = rgg(:)'*dpcurrtdx/cpasma-rcur/cpasma*dcpasmadx;
pcurrt = pcurrt + (rcur-rcur2)*pc2dr;
dpcurrtdx = dpcurrtdx + pc2dr(:)*(drcurdx-drcur2dx);

% Make the current centroid for pcurrt match zcur
pc2dz = ([zeros(1,nr);pcurrt(1:nz-1,:)]-[pcurrt(2:nz,:);zeros(1,nr)])/(2*dz);
zcur2 = zgg(:)'*pcurrt(:)/cpasma;
dzcur2dx = zgg(:)'*dpcurrtdx/cpasma-zcur/cpasma*dcpasmadx;
pcurrt = pcurrt + (zcur-zcur2)*pc2dz;
dpcurrtdx = dpcurrtdx + pc2dz(:)*(dzcurdx-dzcur2dx);

rcurX = rgg(:)'*pcurrt(:)/cpasma;
zcurX = zgg(:)'*pcurrt(:)/cpasma;
drcurXdx = rgg(:)'*dpcurrtdx/cpasma-rcur/cpasma*dcpasmadx;
dzcurXdx = zgg(:)'*dpcurrtdx/cpasma-zcur/cpasma*dcpasmadx;

dpsizrdx = dpsizrpladpcurrt*dpcurrtdx;
dpsizrdx(:,indis) = dpsizrdx(:,indis) + [mpc mpv];

psitor = iplasma.*AR.*fpolzr; % Toroidal flux *in the plasma*

psizr_pla = reshape(dpsizrpladpcurrt*pcurrt(:),nz,nr);  

psizr = psizr_app + psizr_pla;


psiplaapp = psizr_app(:)'*pcurrt(:)/cpasma;
dpsiplaappdx = ...
  (pcurrt(:)'*[mpc mpv zeros(ngg,ns+1)] + ...
  psizr_app(:)'*dpcurrtdx)/cpasma - ...
  psiplaapp/cpasma*dcpasmadx;

% Flux at conductors
ys = mss*[ic;iv] + [mpc mpv]'*pcurrt(:);
% The response of flux at conductors
dysdx = [mpc mpv]'*dpcurrtdx;
dysdx(1:nc+nv,1:nc+nv) = dysdx(1:nc+nv,1:nc+nv)+mss;

gs_nulls
gscp_output
gscp_output_response
y = lae.y;

new_response_was_calculated = true;
if ~plasma_is_tiny
 % plasma_is_tiny = r0/a0 > 50;
end

% Remember last analyzed equilibrium
lae.plasma  = plasma;
lae.psiherr = psiherr;
lae.psizr   = psizr;
lae.psih    = psih;
lae.r0      = r0;
lae.z0      = z0;
lae.a0      = a0;
lae.ic      = ic;
lae.iv      = iv;
lae.sp      = sp;
lae.sf      = sf;
lae.li      = li;
lae.Wth     = Wth;
lae.betap   = betap;
lae.cpasma  = cpasma;
lae.rcur    = rcur;
lae.zcur    = zcur;
lae.psibry  = psibry;
lae.psimag  = psimag;
lae.rsurf   = r0;
lae.aminor  = a0;
lae.rbbbs   = rbbbs(:);
lae.zbbbs   = zbbbs(:);
lae.nbbbs   = nbbbs_max;
lae.rbdef   = rbdef;
lae.zbdef   = zbdef;
lae.rmaxis  = rmaxis;
lae.zmaxis  = zmaxis;
lae.ys      = ys;
