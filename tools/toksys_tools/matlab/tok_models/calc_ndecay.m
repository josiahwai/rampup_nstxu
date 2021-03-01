
   function ndecay = calc_ndecay(psivac,rg,zg,r0,z0,rdec,zdec);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  ndecay = calc_ndecay(psivac,rg,zg,r0,z0,rdec,zdec);
%
%  PURPOSE:  Calculate decay index at point (or vector of points) from
%		vacuum flux on grid.
%
%  INPUTS:
%	psivac = flux on grid (either vector or nz x nr array) (Wb)
%	rg = grid radial position vector [m]
%	zg = grid vertical position vector [m]
%       r0 = radial position of reference point for Bz0 (scalar) [m]
%       z0 = vertical position of reference point for Bz0 (scalar) [m]
%	rdec = radial position(s) to calc decay index at (scalar/vector) [m]
%	zdec = vertical position(s) to calc decay index at (scalar/vector) [m]
%
%  OUTPUTS:
%	ndecay = decay index (indices if rdec, zdec are vectors)
%
%  RESTRICTIONS:
%       r0,z0 must be scalar. rdec,zdec can be vectors.
%	If only r0,z0 are  given (no rdec,zdec), will use that point as
%	both Bz0 reference and location to calculate ndecay. But if
%	rdec,zdec are given as well, r0,z0 will be location of Bz0 reference,
%	and rdec,zdec will be location(s) at which to calc ndecay. 
%
%  METHOD:  
%	see rzrig.m, calc_decind_victor.m

%  WRITTEN BY:  Dave Humphreys 	ON	9/27/07
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prelims:
   mu0 = 0.4*pi;

% Derived Values:
    nr = length(rg);
    nz = length(zg);
    rgg = ones(nz,1)*rg';
    rggv = rgg(:);
    twopi = 2*pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct fields from flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  psivac = reshape(psivac,nz,nr); 
  [dfgdr, dfgdz] = gradient(psivac,rg,zg);
  br = - dfgdz./(twopi*rgg);
  bz = + dfgdr./(twopi*rgg);

  bz0 =  interp2(rg,zg,bz,r0,z0);
  [tmp, dbrdz] = gradient(br,rg,zg);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calc decay index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ndig = - rggv.*dbrdz(:)./bz0;
  if (exist('rdec')~=1 | exist('zdec')~=1), rdec=r0; zdec=z0; end
  ndecay = interp2(rg,zg,reshape(ndig,nz,nr),rdec,zdec);


