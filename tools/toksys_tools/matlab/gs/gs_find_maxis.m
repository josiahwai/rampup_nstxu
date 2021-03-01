%  USAGE:   gs_find_maxis
%
%  PURPOSE: Find magnetic axis (from scratch)
%
%  INPUTS: psizr, the flux at nz vertical positions * nr radial positions
%          ilimgg, flags for grid, -1 = outside vessel, >-1 = inside
%          psimag, psibry, assumed fluxes (only sign(psibry-psimag) matters)
%          The function isinpoly must have been called with limiter information:
%          isinpoly([],[],Rlim,Zlim)
%
%  OUTPUTS: rmaxis, zmaxis, position of magnetic axis (nan if no axis found)
%           wa, iia, weights and indices such that psimag = wa*psizr(iia)
%           drmaxisdpsi, weights such that drmaxis = drmaxisdpsi*dpsizr(iia)
%           dzmaxisdpsi, weights such that dzmaxis = dzmaxisdpsi*dpsizr(iia)
%           plasma_current_has_flipped_sign, boolean true if assumed 
%             sign(psibry-psimag) was wrong
%	
%  METHOD: Search for local minimum in normalized flux, 
%            and call gs_update_maxis to zoom in on it
%          If no minimum is found then search for local maximum and set flag:
%            plasma_current_has_flipped_sign if such an axis is found
	
%  NOTES:  
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	4/12/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmaxis = nan;
zmaxis = nan;

plasma_current_has_flipped_sign = false;
for j = 1:ngg
  if ilimgg(j) > -1 & j > nz & j <= ngg-nz & ...
     (psimag-psibry)*psizr(j-nz) < (psimag-psibry)*psizr(j) & ...
     (psimag-psibry)*psizr(j+nz) < (psimag-psibry)*psizr(j) & ...
     (psimag-psibry)*psizr(j- 1) < (psimag-psibry)*psizr(j) & ...
     (psimag-psibry)*psizr(j+ 1) < (psimag-psibry)*psizr(j)
    rmaxis = rgg(j);
    zmaxis = zgg(j);
  end
end

if isnan(rmaxis) % In this case search for a magnetic axis for opposite current
  for j = 1:ngg
    if ilimgg(j) > -1 & j > nz & j <= ngg-nz & ...
       (psimag-psibry)*psizr(j-nz) > (psimag-psibry)*psizr(j) & ...
       (psimag-psibry)*psizr(j+nz) > (psimag-psibry)*psizr(j) & ...
       (psimag-psibry)*psizr(j- 1) > (psimag-psibry)*psizr(j) & ...
       (psimag-psibry)*psizr(j+ 1) > (psimag-psibry)*psizr(j)
      rmaxis = rgg(j);
      zmaxis = zgg(j);
      plasma_current_has_flipped_sign = true;
    end
  end
end

if isnan(rmaxis)
  cpasma = 0;
else
  gs_update_maxis
end
