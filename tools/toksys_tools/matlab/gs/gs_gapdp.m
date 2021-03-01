%  USAGE:   gs_gapdp
%
%  PURPOSE: Find gaps between boundary and diagnostic points rdp, zdp
%           Gap distances measured along vector from axis to dp points
%
%  INPUTS: psibarzr, psibarzr = (psizr-psimag)/(psibry-psimag);
%          rdp, zdp, diagnostic points
%          rmaxis, zmaxis, position of axis
%          rg, zg, dr, dz, nr, nz (grid variables)
%          For cubic interpolation on the grid:
%            mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%            neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
%
%  OUTPUTS: gapdp, distance from dp to boundary in direction toward axis
%           rdpb, zdpb, = R, Z, of boundary points
%           idpb, indices to 16 grid points around each of rdpb, zdpb
%           wdpb, weights for 16 grid points around each of rdpb, zdpb
%           dgapdpdpsi, how gapdp respondes to (dpsizr(idpb) - dpsibry)
%	
%  METHOD: Newton-Rhapson method used on all contour-points in parallel
	
%  NOTES:  To-do: deduce rhodpb min, max after each iteration
%          to safe-guard against points shooting off?
%          Right now a speed limit of dr/3 per iteration is imposed
	
%
%  WRITTEN BY: Anders Welander ON 2015-02-23
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~plasma
  gapdp = nan(ndp,1);
  return
end

rhodp = sqrt((rdp-rmaxis).^2+(zdp-zmaxis).^2);
thetadp = angle(rdp-rmaxis + 1i*(zdp-zmaxis));
% thbbbs goes from -pi to +pi but not guaranteed to be monotonic
rhodpb = interp1([thbbbs(nbbbs-1)-2*pi; thbbbs(1:nbbbs)],...
                [rhobbbs(nbbbs-1); rhobbbs(1:nbbbs)], thetadp);

cosdpb = cos(thetadp);
sindpb = sin(thetadp);
rdpb = rmaxis+rhodpb.*cosdpb;
zdpb = zmaxis+rhodpb.*sindpb;

contour_iterations = 0;
drhomax = dr;
while drhomax > dr/1e9 & contour_iterations < 50
contour_iterations = contour_iterations+1;
kdpb = floor((rdpb-rg(1))/dr)*nz+floor((zdpb-zg(1))/dz)+1;
xdpb = (rdpb-rgg(kdpb))/dr;
wrc = [ones(ndp,1) xdpb(:) xdpb(:).^2 xdpb(:).^3]*mx;
wrc_r = [zeros(ndp,1) ones(ndp,1) 2*xdpb(:) 3*xdpb(:).^2]*mx/dr;

xdpb = (zdpb-zgg(kdpb))/dz;
wzc = [ones(ndp,1) xdpb(:) xdpb(:).^2 xdpb(:).^3]*mx;
wzc_z = [zeros(ndp,1) ones(ndp,1) 2*xdpb(:) 3*xdpb(:).^2]*mx/dz;

idpb = [kdpb(:)-(nz+1)   kdpb(:)-nz     kdpb(:)-(nz-1)   kdpb(:)-(nz-2) ...
        kdpb(:)-1        kdpb(:)        kdpb(:)+1        kdpb(:)+2      ...
	kdpb(:)+(nz-1)   kdpb(:)+nz     kdpb(:)+(nz+1)   kdpb(:)+(nz+2) ...
	kdpb(:)+(2*nz-1) kdpb(:)+(2*nz) kdpb(:)+(2*nz+1) kdpb(:)+(2*nz+2)];

wdpb = [wzc(:,1).*wrc(:,1) wzc(:,2).*wrc(:,1) wzc(:,3).*wrc(:,1) wzc(:,4).*wrc(:,1) ...
	wzc(:,1).*wrc(:,2) wzc(:,2).*wrc(:,2) wzc(:,3).*wrc(:,2) wzc(:,4).*wrc(:,2) ...
	wzc(:,1).*wrc(:,3) wzc(:,2).*wrc(:,3) wzc(:,3).*wrc(:,3) wzc(:,4).*wrc(:,3) ...
	wzc(:,1).*wrc(:,4) wzc(:,2).*wrc(:,4) wzc(:,3).*wrc(:,4) wzc(:,4).*wrc(:,4)];

wdpb_r = [wzc(:,1).*wrc_r(:,1) wzc(:,2).*wrc_r(:,1) wzc(:,3).*wrc_r(:,1) wzc(:,4).*wrc_r(:,1) ...
	  wzc(:,1).*wrc_r(:,2) wzc(:,2).*wrc_r(:,2) wzc(:,3).*wrc_r(:,2) wzc(:,4).*wrc_r(:,2) ...
	  wzc(:,1).*wrc_r(:,3) wzc(:,2).*wrc_r(:,3) wzc(:,3).*wrc_r(:,3) wzc(:,4).*wrc_r(:,3) ...
	  wzc(:,1).*wrc_r(:,4) wzc(:,2).*wrc_r(:,4) wzc(:,3).*wrc_r(:,4) wzc(:,4).*wrc_r(:,4)];

wdpb_z = [wzc_z(:,1).*wrc(:,1) wzc_z(:,2).*wrc(:,1) wzc_z(:,3).*wrc(:,1) wzc_z(:,4).*wrc(:,1) ...
	  wzc_z(:,1).*wrc(:,2) wzc_z(:,2).*wrc(:,2) wzc_z(:,3).*wrc(:,2) wzc_z(:,4).*wrc(:,2) ...
	  wzc_z(:,1).*wrc(:,3) wzc_z(:,2).*wrc(:,3) wzc_z(:,3).*wrc(:,3) wzc_z(:,4).*wrc(:,3) ...
	  wzc_z(:,1).*wrc(:,4) wzc_z(:,2).*wrc(:,4) wzc_z(:,3).*wrc(:,4) wzc_z(:,4).*wrc(:,4)];

psibardpb = reshape(sum(wdpb'.*psibarzr(idpb)')',ndp,1);
psibardpb_r = reshape(sum(wdpb_r'.*psibarzr(idpb)')',ndp,1);
psibardpb_z = reshape(sum(wdpb_z'.*psibarzr(idpb)')',ndp,1);
psibardpb_rho = cosdpb.*psibardpb_r + sindpb.*psibardpb_z;

drhodpb = (1-psibardpb) ./ psibardpb_rho;
drhodpb(drhodpb > dr/3) = dr/3;
drhodpb(drhodpb < -dr/3) = -dr/3;
rhodpb = rhodpb + drhodpb;
rdpb = rmaxis+rhodpb.*cosdpb;
zdpb = zmaxis+rhodpb.*sindpb;
dpout = ~isinpoly(rdpb,zdpb); % True for points outside limiter
rhodpb(dpout) = rhodpb(dpout) - drhodpb(dpout); % Cancel change for such points

drhomax = max(abs(drhodpb));

end

dum = mean(abs(psibardpb_rho));
for i = 1:ndp
  if psibardpb_rho(i) >= 0 & psibardpb_rho(i) < dum/10000
    psibardpb_rho(i) = dum/10000;
  elseif psibardpb_rho(i) <= 0 & psibardpb_rho(i) > -dum/10000
    psibardpb_rho(i) = -dum/10000;
  end
end

dgapdpdpsi = 1./psibardpb_rho/(psibry-psimag);
gapdp = rhodp-rhodpb;

% The gapdp response is made of a response along rhodp and axis displacement
% Axis displacement causes a displacement of dpb if different from dp
% Hence, the axis displacement only matters when gapdp ~= 0
% dpb moves by axis displacement as: 
%[sindpb*drmaxis,-cosdpb*dzmaxis]*(rhodp-rhodpb)/rhodp
% This changes the flux by: 
% (wdpb_r*sindpb*drmaxis-wdpb_z*cosdpb*dzmaxis)*(rhodp-rhodpb)/rhodp
