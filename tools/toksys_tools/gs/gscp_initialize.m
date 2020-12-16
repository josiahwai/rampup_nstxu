%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_initialize_circle
%
%  PURPOSE: Create an initial circular plasma
%
%  INPUTS: rbdef, zbdef, boundary-defining point
%          lim, index into rl, zl where bdef is in range lim:lim+1
%
%  OUTPUTS: Initialized equilibrium

%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  METHOD: 
	
%  NOTES:  
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander ON 2015-02-27
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plasma

  % Use rmaxis, zmaxis to calculate a0, could alternatively take aminor
  r0 = sum(RAcell(:))/sum(Acell(:));
  z0 = zcur;
  a0 = sqrt((rbdef-r0)^2+(zbdef-z0)^2);

  % ilim is index into Rlim, Zlim with bdef within ilim:ilim+1
  % Search for point on limiter closest to bdef
  dum = inf;
  for i = 1:nlim-1
    R = Rlim(i);
    Z = Zlim(i);
    dR = Rlim(i+1)-Rlim(i);
    dZ = Zlim(i+1)-Zlim(i);
    xl = -((R-rbdef)*dR+(Z-zbdef)*dZ)/(dR^2+dZ^2);
    rt = R + xl*dR;
    zt = Z + xl*dZ;
    s2 = (rt-rbdef)^2 + (zt-zbdef)^2;
    if s2 < dum & xl >= 0 & xl < 1
      dum = s2;
      ilim = i;
      xlim = xl;
      dlim = sqrt(dR^2+dZ^2);
      ulim = [dR dZ]/dlim;
      vlim = [turnin*ulim(:)]';
      rbdef = rt;
      zbdef = zt;
      dRlim = dR;
      dZlim = dZ;
      d2lim = dR^2+dZ^2;
      r0 = rbdef + vlim(1)*a0;
      z0 = zbdef + vlim(2)*a0;
    end
  end

  % Some extra numbers for the boundary
  rmin = r0 - a0;
  rmax = r0 + a0;
  zmin = z0 - a0;
  zmax = z0 + a0;

  % The moving coordinates where flux is measured
  rh = linspace(rmin,rmax,nrh);
  drh = (rmax-rmin)/(nrh-1);

  % Applied, plasma, and total flux at rh
  psihapp = gs_interp2(rg, zg, psizr_app, rh, z0);
  psihpla = cpasma*mfcircle2point(rh, 0, r0, 0, a0);
  psih = gs_interp2(rg, zg, psizr, rh, z0);;

else % There is no plasma
  r0 = nan;
  z0 = nan;
  a0 = 0;
  Atot = 0;
end
