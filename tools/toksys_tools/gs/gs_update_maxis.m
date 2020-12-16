%  USAGE:   gs_update_maxis
%
%  PURPOSE: Zoom in on magnetic axis
%
%  INPUTS: rmaxis, zmaxis, a best estimate of the position such as found
%            by a linear prediction of change since previous analysis
%          psizr, the flux at nz vertical positions * nr radial positions
%          rgg, zgg, dr, dz, nr, nz (grid variables)
%          For bicubic interpolation on the grid:
%            mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%            neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
%          The function isinpoly must have been called with limiter information:
%          isinpoly([],[],Rlim,Zlim)
%
%  OUTPUTS: maxis, a boolean true if axis still inside limiter after update
%           rmaxis, zmaxis, position of magnetic axis
%           wa, iia, weights and indices such that psimag = wa*psizr(iia)
%           drmaxisdpsi, weights such that drmaxis = drmaxisdpsi*dpsizr(iia)
%           dzmaxisdpsi, weights such that dzmaxis = dzmaxisdpsi*dpsizr(iia)
%	
%  METHOD: interpolation with bicubic Hermite splines, Newton-Rhapson to zoom
	
%  NOTES:  
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	4/12/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Search for magnetic axis with Newton-Rhapson
twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
j = 19;
while j > 0 & twopirbrzmax > 1e-10
  kr0 = min(nr-3,max(1,floor((rmaxis-rg(1))/dr)));
  kz1 = min(nz-2,max(2,ceil((zmaxis-zg(1))/dz)));
  k = kr0*nz+kz1;
  iia = k+neighbors'; % iia indexes 16 points around magnetic axis
  pp = psizr(iia);
  tr = (rmaxis-rgg(k))/dr;
  tz = (zmaxis-zgg(k))/dz;
  wa = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16); % weights to get flux at axis
  war = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16); % weights to get dpsi/dr
  waz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16); % weights to get dpsi/dz
  warr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16); % weights to get d2psi/dr2
  wazz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16); % weights to get d2psi/dr2
  warz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16); % weights to get d2psi/drdz
  ashift_Ba = -inv([warr*pp warz*pp; warz*pp wazz*pp]);
  cmaxis = ashift_Ba*[war*pp; waz*pp]; % Correction of position
  dum = 4*((cmaxis(1)/dr)^2+(cmaxis(2)/dz)^2);
  if dum > 1
    cmaxis = cmaxis/sqrt(dum);
  end
  rmaxis = rmaxis+cmaxis(1);
  zmaxis = zmaxis+cmaxis(2);
  j = j-1;
  twopirbrzmax = max(abs([war*pp waz*pp]));
end

drmaxisdpsi = ashift_Ba(1,1)*war+ashift_Ba(1,2)*waz;
dzmaxisdpsi = ashift_Ba(2,1)*war+ashift_Ba(2,2)*waz;

maxis = norm(cmaxis) < dr && ilimgg(k) > -1 && isinpoly(rmaxis,zmaxis);
