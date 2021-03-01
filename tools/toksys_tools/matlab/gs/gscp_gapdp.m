%  USAGE:   gscp_gapdp
%
%  PURPOSE: Return gapdp for circular plasma model
%           By definition gapdp points toward the magnetic axis but
%           here the circle center is used instead since it is easier
%
%  INPUTS: r0, z0, a0
%
%  OUTPUTS: gapdp, distance from dp to boundary in direction toward axis
%           rdpb, zdpb, = R, Z, of boundary points
%           idpb, indices to 16 grid points around each of rdpb, zdpb
%           wdpb, weights for 16 grid points around each of rdpb, zdpb
%           dgapdpdpsi, how gapdp respondes to (dpsizr(idpb) - dpsibry)
%	

%  METHOD: 
	
%  NOTES:  
	
%  VERSION %W% %G%
%
%  WRITTEN BY: Anders Welander ON 2015-03-04
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~plasma
  gapdp = nan(ndp,1);
  return
end

thetadp = angle(rdp-r0 + 1i*(zdp-z0));
rhodpb = a0;

cosdpb = cos(thetadp);
sindpb = sin(thetadp);
rdpb = r0+rhodpb.*cosdpb;
zdpb = z0+rhodpb.*sindpb;

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
rdpb = rdpb + drhodpb.*cosdpb;
zdpb = zdpb + drhodpb.*sindpb;

drhomax = max(abs(drhodpb));

dum = mean(abs(psibardpb_rho));
for i = 1:ndp
  if psibardpb_rho(i) >= 0 & psibardpb_rho(i) < dum/100
    psibardpb_rho(i) = dum/100;
  elseif psibardpb_rho(i) <= 0 & psibardpb_rho(i) > -dum/100
    psibardpb_rho(i) = -dum/100;
  end
end

dgapdpdpsi = -1./psibardpb_rho;
rhodp = sqrt((rdp-r0).^2+(zdp-z0).^2);
gapdp = rhodp-rhodpb;
