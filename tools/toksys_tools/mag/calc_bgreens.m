function [Br,Bz,dBrdr,dBrdz,dBzdr,dBzdz] = calc_bgreens(ra,rb,za,zb,k,Ek,Kk)
 %
%  SYNTAX: [Br,Bz,dBrdr,dBrdz,dBzdr,dBzdz] = calc_bgreens(ra,rb,za,zb,k,Ek,Kk)
%
%  PURPOSE:  Calculate Br, Bz green functions with the option of using already
%  computed elliptic integrals.
%	[Kk,Ek]=ellipke(m)
%       k=sqrt(m);
%  where  
%	m=((rb+ra).^2 - (rb-ra).^2)./((rb+ra).^2 + (za-zb).^2);
%
%  INPUT:
%	ra = radius of source element (m)
%	rb = radius of measurement point (m)
%	za = z position of source element (m)
%	zb = z position of measurement point (m)
%	k = parameter that is argument to elliptic integrals (optional)
%	Ek = calculated value of elliptic integral E(k) (optional)
%	Kk = calculated value of elliptic integral K(k)	 (optional)
%   All inputs must be same size (scalars, vectors, or matrices).
%   If k,Ek,Kk not specified, they are computed here.
%
%  OUTPUT:
%	Br = radial field at measurement point due to current in filament 
%		at source location
%	Bz = vertical field at measurement point due to current in filament 
%		at source location
%       dBrdr = derivative of Br w.r.t. rb at measurement point
%       dBrdz = derivative of Br w.r.t. zb at measurement point
%       dBzdr = derivative of Bz w.r.t. rb at measurement point
%       dBzdz = derivative of Bz w.r.t. zb at measurement point
%
%  RESTRICTIONS: This is a "filament-to-filament" calculation and is not correct
%  for distributed current sources.  The filament calculation is singular for 
%  "point b" equal to "point a".
 
%  METHOD:  All formulas can be found in R.A.Schill, General Relation for the
%   Vector Magnetic Field of a Circular Current Loop: A Closer Look, IEEE Trans.
%   Magnetics, vol.39, no.2, March 2003
%
%  WRITTEN BY:  Mike Walker 	ON 	8/22/97 (from IDL version)
%  (from Dave Humphreys' mindbf.pro and psical in efund6565.f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%help, ra, rb, za, zb, k, Ek, Kk, br, bz

% jal5-23-00: br blows up for rb=0 answer should be br=0 below is clug
  id= find(rb<eps);
  rb(id)= eps*ones(size(id));

if nargin < 7
   m=((rb+ra).^2 - (rb-ra).^2)./((rb+ra).^2 + (za-zb).^2);
   k=sqrt(m);
   [Kk,Ek]=ellipke(m);
end

twopir = 2.0*pi.*rb;

% Here are EFIT calculations (in psical.f).  I checked - these are exactly
% the analytically reduced form of Dave's calculations.
      delr = rb-ra;
      delz = zb-za;
      d = (rb+ra).^2 + delz.^2;

      G = (-Kk + (ra.^2+rb.^2+delz.^2)./(delr.^2+delz.^2).*Ek);
      H = (Kk + (ra.^2-rb.^2-delz.^2)./(delr.^2+delz.^2).*Ek);

      Br = 2e-7 * delz./(rb.*sqrt(d)).*G;
      Bz = 2e-7 * H./sqrt(d);

if(nargout > 2) 	% Calculations of derivatives from Schill paper

   dkdrb = k./(2*rb).*(1-(k.^2/2).*(1+rb./ra));
   dkdzb = -(k.^3.*delz)./(4*ra.*rb);

if 0	% testing
figure(1),clf
plot(rb,dkdrb)
hold on
plot(rb(2:end),diff(k)/mean(diff(rb)),'r--')
hold off
title('dkdrb')
end

   dKdk = (Ek./(1-k.^2) - Kk)./k;		% OK
   dEdk = (Ek - Kk)./k;			% OK

if 0	% testing
figure(1),clf
plot(k,dKdk)
hold on
plot(k(2:end),diff(Kk)./diff(k),'r--')
hold off
title('dKdk')

figure(2),clf
plot(k,dEdk)
hold on
plot(k(2:end),diff(Ek)./diff(k),'r--')
hold off
title('dEdk')
end

   Gnum = ra.^2 + rb.^2 + delz.^2;
   Hnum = ra.^2 - rb.^2 - delz.^2;
   denom = (delr.^2+delz.^2);

   dGdrb= -dKdk.*dkdrb + 2*((delr.^2+delz.^2).*rb-Gnum.*delr)./denom.^2.*Ek ...
	     + Gnum./denom.*dEdk.*dkdrb;
   dGdzb= -dKdk.*dkdzb +2*((delr.^2+delz.^2).*delz-Gnum.*delz)./denom.^2.*Ek ...
	     + Gnum./denom.*dEdk.*dkdzb;
   dHdrb= dKdk.*dkdrb - 2*((delr.^2+delz.^2).*rb+Hnum.*delr)./denom.^2.*Ek ...
	     + Hnum./denom.*dEdk.*dkdrb;
   dHdzb= dKdk.*dkdzb - 2*((delr.^2+delz.^2).*delz+Hnum.*delz)./denom.^2.*Ek ...
	     + Hnum./denom.*dEdk.*dkdzb;

if 0	% OK
figure(1),clf
plot(rb,dGdrb)
hold on
plot(rb(2:end),diff(G)./diff(rb),'r--')
hold off
title('dGdr')

figure(1),clf
plot(rb,dHdrb)
hold on
plot(rb(2:end),diff(H)./diff(rb),'r--')
hold off
title('dHdr')
end

if 0	% testing
figure(1),clf
plot(zb,dGdzb)
hold on
plot(zb(2:end),diff(G)./diff(zb),'r--')
hold off
title('dGdz')

figure(1),clf
plot(zb,dHdzb)
hold on
plot(zb(2:end),diff(H)./diff(zb),'r--')
hold off
title('dHdz')
end

   rtd = sqrt(d);

   dBrdr=2e-7*(delz./(rb.*rtd).*dGdrb - ...
		((rb.*(ra+rb)+d).*delz)./(rb.^2.*rtd.^3).*G);
   dBrdz = 2e-7*((ra+rb).^2./(rb.*rtd.^3).*G + (delz./(rb.*rtd)).*dGdzb);
   dBzdr = 2e-7*(dHdrb./rtd - (ra+rb)./(rtd.^3).*H);
   dBzdz = 2e-7*(dHdzb./rtd - delz./(rtd.^3).*H);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING TESTING TESTING TESTING TESTING TESTING TESTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if (0) 
%      dkdk=Ek/(k*(1.-k^2))-Kk/k
%      dedk=(Ek-Kk)/k
%      dmdk=8.0*pi*sqrt(ra*rb)*1.0e-7*((Ek-(1.-.5*k^2)*Kk)/k^2  ...
%          +((1.-.5*k^2)*dkdk-k*Kk-dedk)/k)
%      dkdz=(za-zb)*sqrt(((rb+ra)^2-(rb-ra)^2)/((rb+ra)^2+(za-zb)^2)^3)
%      br0 = -dmdk*dkdz/twopir
%      d = (rb+ra)^2 + (za-zb)^2
%      n = (rb+ra)^2 - (rb-ra)^2
%      dddr = 2.*(rb + ra)
%      dndr = 4.*ra
%      dkdr = (.5/k)*(d*dndr - n*dddr)/d^2
%      dmdr = 8.0*pi*0.5*sqrt(ra/rb)*1.0e-7*((1.-.5*k^2)*Kk-Ek)/k
%      bz0 = (dmdk*dkdr + dmdr)/twopir;
%end

%if (0) 
%print,'normalized norm(br-br0) = ',norm(br-br0)/norm(br)
%print,'normalized norm(bz-bz0) = ',norm(bz-bz0)/norm(bz)
%if (norm(br-br0)/norm(br) > 1e-3 |  ...
%    norm(bz-bz0)/norm(bz) > 1e-3) 
%   window,0
%   plot,br-br0
%   window,1
%   plot,br
%   oplot,br0
%   wait
%end
%end
