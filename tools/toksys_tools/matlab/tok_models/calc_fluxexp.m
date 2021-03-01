function [fluxexp, r2, z2] = calc_fluxexp(eq, rr, zz, dd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   fluxexp = calc_fluxexp(eq, r, z)
%           [fluxexp, r2, d2] = calc_fluxexp(eq, r, z, d)
%
%  PURPOSE: Calculate flux expansion at points r,z for equilibrium eq
%
%  INPUTS: eq, TokSys description of an equilibrium containing fields:
%              rg, zg = grid coordinates
%              psizr = flux at points rg, zg
%              rbbbs, zbbbs, nbbbs = nbbbs plasma boundary coordinates
%              psibry = boundary flux
%          r,z, coordinates of points where flux expansion is calculated
%          d, distance in outboard midplane to outside boundary (default 0) 
%
%  OUTPUTS: fluxexp, if d==0,flux expansion = abs(pr)./sqrt(yr.^2+yz.^2)
%    where pr = d(psi)/dR at maximum major radius for the boundary
%    and sqrt(yr.^2+yz.^2) is magnitude of flux gradient at points r,z
%    if d > 0, fluxexp = D/d, where D is perpendicular distance from r,z
%    to point r2, z2 with same flux difference as distance d in midplane
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Anders Welander ON 2015-04-27
%
%  MODIFICATION HISTORY: ASW 20171129 allow arbitrary r,z,d dimensions				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4 | isempty(dd)
  dd = 0;
end

% Find where boundary is at maximum major radius
[rmax, i] = max(eq.rbbbs(1:eq.nbbbs));
zmid = eq.zbbbs(i);
% Solve the equations:
% pz + prz*x + pzz*y = 0
% p + pr*x + pz*y = psibry
% for x and y, and repeat a few times to zoom in on max r of boundary
for i = 1:4
  [p, pr, pz, prr, prz, pzz] = gs_interp2(eq.rg, eq.zg, eq.psizr, rmax, zmid);
  xy = inv([prz pzz; pr pz])*[-pz;eq.psibry-p];
  rmax = rmax+xy(1);
  zmid = zmid+xy(2);
end

% Outboard midplane flux and gradient
[p, pr] = gs_interp2(eq.rg, eq.zg, eq.psizr, rmax, zmid);

% Allocate memory for fluxexp
fluxexp = rr+zz+dd;

for k = 1:numel(fluxexp)
  if k <= numel(rr)
    r = rr(k);
  end
  if k <= numel(zz)
    z = zz(k);
  end
  if k <= numel(dd)
    d = dd(k);
  end
  [y, yr, yz] = gs_interp2(eq.rg, eq.zg, eq.psizr, r, z);
  if d == 0
    fluxexp(k) = abs(pr)./sqrt(yr.^2+yz.^2);
    r2 = nan;
    z2 = nan;
  else % Find point r2, z2 where flux equals y+p2-p
    dr = eq.rg(2)-eq.rg(1);
    dz = eq.zg(2)-eq.zg(1);
    p2 = gs_interp2(eq.rg, eq.zg, eq.psizr, rmax+d, zmid);
    % y+p2-p = y2 + yr2*yr*t + yz2*yz*t
    r2 = r;
    z2 = z;
    y2 = y;
    yr2 = yr;
    yz2 = yz;
    for it = 1:5
      t = (y-y2+p2-p)./(yr.*yr2+yz.*yz2);
      Dr = abs(yr.*t);
      Dz = abs(yz.*t);
      % Don't take full Newton-Rhapson step if larger than half a grid cell
      for i = 1:length(t(:))
	t(i) = min([1 dr/Dr(i)/2 dz/Dz(i)/2])*t(i);
      end
      if it == 1
	t(yr.*t<0) = -t(yr.*t<0); % Go to larger radii, will add tracing here
      end
      r2 = r2 + yr.*t;
      z2 = z2 + yz.*t;
      [y2, yr2, yz2] = gs_interp2(eq.rg, eq.zg, eq.psizr, r2, z2);
    end
    D = sqrt((r2-r).^2+(z2-z).^2);
    fluxexp(k) = D/d;
  end
end



