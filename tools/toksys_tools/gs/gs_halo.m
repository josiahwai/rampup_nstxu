%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_halo
%
%  PURPOSE: Calculate halo current on grid
%
%  INPUTS:  halo, index to a halo model
%           Several equilibrium variables (made by gs_eq_analysis)
%
%  OUTPUTS: ihzr, halo current on grid
	
%  VERSION @(#)gs_halo.m	1.1 02/23/15
%
%  WRITTEN BY:  Anders Welander  ON	2/14/15
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if halo == 1
  reszr = Bp2zr;
  resmean = mean(reszr(ilimgg > -1));
  reszr(reszr < 0.01*resmean) = 0.01*resmean;
  ihzr = geohzr./reszr.*(Ag-Acell);
  ihzr = ihalo/sum(ihzr(:))*ihzr;
end
