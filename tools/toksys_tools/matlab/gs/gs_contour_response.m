%  USAGE:   gs_contour_response
%
%  PURPOSE: Calculate response of contours rcont,zcont to x
%
%  INPUTS: outputs from gs_trace_contours & gs_response
%
%  OUTPUTS: drcontdx, dzcontdx, dthabdx
%	
%  METHOD: Some vector algebra and derivatives
	
%  NOTES:  
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	5/11/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~plasma
  return
end

rhoab = sqrt((rbdef-rmaxis)^2+(zbdef-zmaxis)^2);
dthabdx = (cos(thab)*(dzbdefdx-dzmaxisdx)-sin(thab)*(drbdefdx-drmaxisdx))/rhoab;

rhocont = sqrt((rcont-rmaxis).^2+(zcont-zmaxis).^2);

drcontdx = zeros(npola*ncont,nx);
dzcontdx = zeros(npola*ncont,nx);

for i = 1:npola*ncont
  
  if psibarcont_target(i) == 0

    % The response of rcont
    drcontdx(i,:) = drmaxisdx;

    % The response of zcont
    dzcontdx(i,:) = dzmaxisdx;

  else

    % The change of psibar at the original rcont, zcont 
    dumxp1 = ...
      (wcont(i,:)*dpsizrdx(icont(i,:),:) - ...
      (1-psibarcont(i))*dpsimagdx - ...
      psibarcont(i)*dpsibrydx)/(psibry-psimag);

    % Translation & rotation
    dumxr = drmaxisdx - rhocont(i)*sincont(i)*dthabdx;
    dumxz = dzmaxisdx + rhocont(i)*coscont(i)*dthabdx;

    % The change of psibar due to tranlation & rotation
    dumxp2 = psibarcont_r(i)*dumxr + psibarcont_z(i)*dumxz;

    % The change of rho for the translated & rotated vector
    dumrh = -(dumxp1+dumxp2)/psibarcont_rho(i);

    % The response of rcont
    drcontdx(i,:) = dumxr + coscont(i)*dumrh;

    % The response of zcont
    dzcontdx(i,:) = dumxz + sincont(i)*dumrh;
  
  end
  
end

for i = 1:ncont
  if psibarc(i) == 1
    drcontdx((i-1)*npola+1,:) = drbdefdx;
    dzcontdx((i-1)*npola+1,:) = dzbdefdx;
    drcontdx(i*npola,:) = drbdefdx;
    dzcontdx(i*npola,:) = dzbdefdx;
  end
end

