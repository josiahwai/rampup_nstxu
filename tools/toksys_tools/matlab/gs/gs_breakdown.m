%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_breakdown (breakdown meaning formation of embryonic GS plasma)
%
%  PURPOSE: Decide initial plasma position
%
%  INPUTS: ic, iv, coil and vessel currents
%          cpasma, initial plasma current (only its sign is needed)
%
%  OUTPUTS:  plasma, flag if a plasma now exists
%            rmaxis, zmaxis, the position of breakdown
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	7/15/14
%
%  MODIFICATION HISTORY: Latest version 2015-02-18		
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Breakdown can occur in several places but for now only 1 is chosen

% User may specify vicinity of breakdown with init.rmaxis, init.zmaxis
if isnan(rmaxis) | isnan(zmaxis) | ~isinpoly(rmaxis,zmaxis)
  % Otherwise the inboard midplane is chosen:
  zmaxis = (max(zl)+min(zl))/2;
  rmaxis = max(rl);
  for i = 1:nl-1
    if (zl(i)-zmaxis)*(zl(i+1)-zmaxis) <= 0 & zl(i) ~= zl(i+1)
      dum = ((zl(i+1)-zmaxis)*rl(i)-(zl(i)-zmaxis)*rl(i+1))/(zl(i+1)-zl(i));
      if rmaxis > dum
        rmaxis = dum;
      end
    end
  end
end

% Calculate poloidal flux on the grid
psizr_app = reshape(mpc*ic+mpv*iv,nz,nr); % from coils and vessel
if halo % ih is current on open field lines [Amperes/cell]
  psizr_halo = reshape(dpsizrpladpcurrt*ih(:),nz,nr);
else
  psizr_halo = zeros(nz,nr);
end
psizr = psizr_app + psizr_halo;

% Now find a place near rmaxis, zmaxis where breakdown can occur

% Nulls are possible places. These breakdowns give rise to diverted plasmas.
gs_nulls
% Possible points are stored in bdef.candidates

% Breakdown can also occur at the limiter where the radial field is zero
% and the vertical field pushes the plasma into the wall. This makes a circle.
gs_touches

dum = inf;
for i = 1:touches.count
  if touches.inicursign(i)*cpasma > 0
    if dum > (touches.r(i)-rmaxis)^2 + (touches.z(i)-zmaxis)^2
      dum = (touches.r(i)-rmaxis)^2 + (touches.z(i)-zmaxis)^2;
      plasma = 1;
      rbdef = touches.r(i);
      zbdef = touches.z(i);
      psibry = touches.psi(i);
      lim = touches.limiter(i);
    end
  end
end
