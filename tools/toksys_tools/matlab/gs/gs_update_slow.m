%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_update_slow
%
%  PURPOSE: Update equilibrium: psizr, ic, iv, sp, sf and output y
%           
%  INPUTS:  Must have run gs_update_fast
%	
%  OUTPUTS:  Updated:
%            psizr, flux on grid [Wb]
%            ic, coil currents [A]
%            iv, vessel currents [A]
%            sp, spline parameters for pressure versus normalized flux
%            sf, spline parameters for fpol^2/2-(rzero*bzero)^2/2
%            y, diagnostic outputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander ON 2015-03-08
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLASMA CREATION
if ~lae.plasma & plasma
  circle_model = use_circle_model; % Use circle_model if flag says so
  ic = lae.ic + dx(indic); % gs_find_dx was already run in gs_update_fast
  iv = lae.iv + dx(indiv); % Currents ic, iv used to calculate psizr_app
  gs_breakdown % Find a place for the first plasma blip
  if circle_model
    gscp_creation
    er = 0;
  else
    gs_creation
    gs_eq_analysis
    gs_response
    er = 0;
  end
  plasma_is_tiny = true;
  gs_update_fast % Rerun gs_update_fast with this new plasma
  if equilibrium_update_is_complete
    return
  end
end

% PLASMA EXTINCTION
if lae.plasma & ~plasma 
  ic = lae.ic + dx(indic); % gs_find_dx was already run in gs_update_fast
  iv = lae.iv + dx(indiv); % Currents ic, iv used to calculate psizr_app
  gs_eq_analysis
  gs_response
  err.ys(:) = 0;
  err.y(:) = 0;
  equilibrium_update_is_complete = true;
  return
end

% Change equilibrium as much as allowed by fdx
dxallowed = dx*fdx;

% Archive predictions with old equilibrium at this junction
junction.y  = lae.y  + dydx*dxallowed;
junction.ys = lae.ys + dysdx*dxallowed;
junction.ic = lae.ic + dxallowed(indic);
junction.iv = lae.iv + dxallowed(indiv);
junction.li = lae.li + dlidx*dxallowed;
junction.Wth = lae.Wth + dWthdx*dxallowed;
junction.betap = lae.betap + dbetapdx*dxallowed;
junction.cpasma = lae.cpasma + dcpasmadx*dxallowed;

% Analyze the equilibrium at the junction set by fdx
psizr(:) = lae.psizr(:) + dpsizrdx*dxallowed;
ic = lae.ic + dxallowed(indic);
iv = lae.iv + dxallowed(indiv);
sp = lae.sp + dxallowed(indsp);
sf = lae.sf + dxallowed(indsf);
rmaxis = lae.rmaxis + drmaxisdx*dxallowed; % Updated to improve success 
zmaxis = lae.zmaxis + dzmaxisdx*dxallowed; % rate for gs_update_maxis
if circle_model
  psih(:) = lae.psih(:) + dpdx(1:nrh,:)*dxallowed;
  r0 = lae.r0 + dr0dx*dxallowed;
  z0 = lae.z0 + dz0dx*dxallowed;
  a0 = lae.a0 + da0dx*dxallowed;
  gscp_analysis_response
else
  gs_eq_analysis
  gs_response
end

% Transition between circle and grid models as needed
if circle_model & Atot/Ag > ncp2gr
  circle_model = 0;
  gs_eq_analysis
  gs_response
elseif ~circle_model & Atot/Ag < ngr2cp & use_circle_model
  circle_model = 1;
  gscp_configure
  gscp_initialize
  gscp_analysis_response
  gscp_converge
end

% Make the converge error correction of lae in full right away
if circle_model
  dpsihbar = dpsihbardx(:,nx);
  dr0 = dr0dx(nx);
  dz0 = dz0dx(nx);
  da0 = da0dx(nx);
  nlcritnx = max(abs([dpsihbar/dpsibar_limit_for_linear; ...
    [dr0; dz0; da0]/min(a0/5,dr)]));
  % Stay inside grid
  if lae.r0+dr0 + lae.a0+da0 > rg(nr-1)
    nlcritnx = max(nlcritnx, abs((dr0+da0)/(rg(nr-1)-lae.r0-lae.a0)));
  end  
  if lae.r0+dr0 - lae.a0-da0 < rg(2)
    nlcritnx = max(nlcritnx, abs((dr0-da0)/(rg(2)-lae.r0+lae.a0)));
  end
else
  dpsizrx(:) = dpsizrdx(:,nx);
  dpsimagx = wa*dpsizrx(iia);
  dpsibryx = wb*dpsizrx(iib);
  azr = (1-lae.psibarzr)*dpsimagx;
  bzr = lae.psibarzr*dpsibryx;
  dpsibarzrx = (dpsizrx-azr-bzr)/(lae.psibry-lae.psimag);
  dpsibar = max(abs(dpsibarzrx(iplasma)));
  drmaxis = drmaxisdx(nx);
  dzmaxis = dzmaxisdx(nx);
  nlcritnx = max([abs(drmaxis/dr*2) abs(dzmaxis/dz*2) ...
    dpsibar/dpsibar_limit_for_linear]);  
end
% gnx is the largest fraction of dx(nx) that can ever be done
gnx = min(1,1/nlcritnx);
lae.cpasma(:) = lae.cpasma(:) - dcpasmadx(:,nx)*gnx;
lae.psizr(:) = lae.psizr(:) - dpsizrdx(:,nx)*gnx;
lae.betap(:) = lae.betap(:) - dbetapdx(:,nx)*gnx;
lae.Wth(:) = lae.Wth(:) - dWthdx(:,nx)*gnx;
lae.li(:) = lae.li(:) - dlidx(:,nx)*gnx;
lae.ys(:) = lae.ys(:) - dysdx(:,nx)*gnx;
lae.y(:) = lae.y(:) - dydx(:,nx)*gnx;
if circle_model
  lae.psih(:) = lae.psih(:) - dpdx(1:nrh,nx)*gnx;
  lae.r0 = lae.r0 - dr0dx(1,nx)*gnx;
  lae.z0 = lae.z0 - dz0dx(1,nx)*gnx;
  lae.a0 = lae.a0 - da0dx(1,nx)*gnx;
end
er = 0; % Now lae has no formal error

% If plasma quantities other than sp, sf are states then sp, sf
% need be corrected so that the new equilibrium matches the states
if constraints
  if constraints == 1
    dx2 = dxdxc(:,nc+nv+(1:3)) * [junction.cpasma-lae.cpasma; ...
      junction.li-lae.li; junction.betap-lae.betap];
  else % if constraints == 2
    dx2 = dxdxc(:,nc+nv+(1:3)) * [junction.cpasma-lae.cpasma; ...
      junction.li-lae.li; junction.Wth-lae.Wth];
  end
  % Check this dx2 for problems
  if circle_model
    dpsihbar = dpsihbardx*dx;
    dr0 = dr0dx*dx;
    dz0 = dz0dx*dx;
    da0 = da0dx*dx;
    nlcrit2 = max(abs([dpsihbar/dpsibar_limit_for_linear; ...
      [dr0; dz0; da0]/min(a0/5,dr)]));
  else
    dpsizrx(:) = dpsizrdx*dx2;
    dpsimagx = wa*dpsizrx(iia);
    dpsibryx = wb*dpsizrx(iib);
    azr = (1-lae.psibarzr)*dpsimagx;
    bzr = lae.psibarzr*dpsibryx;
    dpsibarzrx = (dpsizrx-azr-bzr)/(lae.psibry-lae.psimag);
    dpsibar = max(abs(dpsibarzrx(iplasma)));
    drmaxis = drmaxisdx*dx;
    dzmaxis = dzmaxisdx*dx;
    nlcrit2 = max([abs(drmaxis/dr*2) abs(dzmaxis/dz*2) ...
      dpsibar/dpsibar_limit_for_linear]);
  end
  fdx2 = min(1,1/nlcrit2);
  lae.y = lae.y + dydx*dx2;
  lae.ys = lae.ys + dysdx*dx2;
  lae.sp = lae.sp + dx2(indsp);
  lae.sf = lae.sf + dx2(indsf);
  lae.li = lae.li + dlidx*dx2;
  lae.betap = lae.betap + dbetapdx*dx2;
  lae.cpasma = lae.cpasma + dcpasmadx*dx2;
  lae.psizr(:) = lae.psizr(:) + dpsizrdx*dx2;
  if circle_model
    lae.psih(:) = lae.psih(:) + dpdx(1:nrh,:)*dx2;
    lae.r0 = lae.r0 + dr0dx*dx2;
    lae.z0 = lae.z0 + dz0dx*dx2;
    lae.a0 = lae.a0 + da0dx*dx2;
  end
else
  fdx2 = 1;
end

if fdx2 == 1
  % Now lae is the (formally error-free) equilibrium at the junction
  % Archive errors in predictions for junction with original equilibrium
  err.y(iy) = junction.y(iy) - lae.y(iy);
  err.ys = junction.ys - lae.ys;
else
  % No sense phasing out errors that may be huge
  err.y(:) = 0;
  err.ys(:) = 0;
end

% Try again to update the equilibrium
gs_update_fast

if ~equilibrium_update_is_complete
  % More than one update is needed to get to targets
  if circle_model
    gscp_converge
  else
    gs_converge
  end
  % At this point it is no longer meaningful to phase out errors
  % because prediction with original equilibrium can be hugely wrong
  % when evaluated at the target state 
  err.y(:) = 0;
  err.ys(:) = 0;
  gs_update_fast
end

