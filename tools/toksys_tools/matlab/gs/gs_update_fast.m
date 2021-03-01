%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_update_fast
%
%  PURPOSE: Update equilibrium: psizr, ic, iv, sp, sf and output y
%           
%  INPUTS:  The gs workspace, and:
%           xc or xs, structure with inputs 
%                ic: coil currents such that psizr_app = mpc*ic
%                iv: vessel currents such that psizr_app = mpv*iv
%                ip: total plasma current
%                li: normalized inductance
%             betap: poloidal beta
%            plotit: Special flag for plotting equilibrium (default 0)
%	
%  OUTPUTS:  Updated:
%            psizr, flux on grid [Wb]
%            ic, coil currents [A]
%            iv, vessel currents [A]
%            sp, spline parameters for pressure versus normalized flux
%            sf, spline parameters for fpol^2/2-(rzero*bzero)^2/2
%            y, diagnostic outputs
%            g, largest fraction of dx(nx) that can be done
%            equilibrium_update_is_complete, if false run gs_update_slow
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  WRITTEN BY:  Anders Welander ON 2015-03-08
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get dx needed from last-analyzed-equilibrium (lae) to reach new state
gs_find_dx

if isfield(index_in_y,'xs') & exist('xs','var') & ...
  length(index_in_y.xs) == length(xs)
  lae.y(index_in_y.xs) = xs;
end

% Determine if the change can be done all at once
if halo
  nlcriterion = inf; % Open field line currents always done by gs_update_slow
elseif ~plasma & lae.plasma
  nlcriterion = inf; % Plasma extinction done by gs_update_slow
elseif ~plasma
  nlcriterion = 0;
elseif ~lae.plasma
  nlcriterion = inf; % Plasma creation done by gs_update_slow
elseif plasma_is_tiny
  nlcriterion = 0;
  if daminordx*dx > 0
    nlcriterion = inf;
    plasma_is_tiny = 0;
  end
elseif circle_model % Using circular plasma model
  dpsihbar = dpsihbardx*dx;
  dr0 = dr0dx*dx;
  dz0 = dz0dx*dx;
  da0 = da0dx*dx;
  nlcriterion = max(abs([dpsihbar/dpsibar_limit_for_linear; ...
    [dr0; dz0; da0]/min(a0/5,dr)]));
  % Stay inside grid
  if lae.r0+dr0 + lae.a0+da0 > rg(nr-1)
    nlcriterion = max(nlcriterion, abs((dr0+da0)/(rg(nr-1)-lae.r0-lae.a0)));
  end  
  if lae.r0+dr0 - lae.a0-da0 < rg(2)
    nlcriterion = max(nlcriterion, abs((dr0-da0)/(rg(2)-lae.r0+lae.a0)));
  end
else
  dpsizrx(:) = dpsizrdx(:,1:nx-1)*dx(1:nx-1);
  dpsimagx = wa*dpsizrx(iia);
  dpsibryx = wb*dpsizrx(iib);
  azr = (1-lae.psibarzr)*dpsimagx;
  bzr = lae.psibarzr*dpsibryx;
  dpsibarzrx = (dpsizrx-azr-bzr)/(lae.psibry-lae.psimag);
  dpsibar = max(abs(dpsibarzrx(iplasma)));
  drmaxis = drmaxisdx*dx;
  dzmaxis = dzmaxisdx*dx;
  nlcriterion = max([abs(drmaxis/dr*2) abs(dzmaxis/dz*2) ...
    dpsibar/dpsibar_limit_for_linear]);
end
% fdx is the fraction of dx(1:nx-1) that can be done
fdx = min(1,1/nlcriterion);
new_aminor = lae.aminor + daminordx*dx;
new_rsurf = lae.rsurf + drsurfdx*dx;
new_plasma_is_tiny = new_aminor/new_rsurf < 1/50;
if fdx > 0 & plasma_is_tiny & new_plasma_is_tiny
  fdx = 1;
end
if fdx == 1 % Simply update the outputs
  if phaseouterrors
    % err is previous minus present equilibrium for state used to make present
    y  = lae.y  + dydx*dx  + (1-nlcriterion)*err.y; % Subtract err gradually 
    ys = lae.ys + dysdx*dx + (1-nlcriterion)*err.ys;
  else
    y  = lae.y  + dydx*dx;
    ys = lae.ys + dysdx*dx;
  end
  if plasma_is_tiny & new_plasma_is_tiny
    dcpasma = dcpasmadx*dx*fdx;
    if dcpasma/lae.cpasma < 0
      rbdef = lae.rbdef + drbdefdx*dx;
      zbdef = lae.zbdef + dzbdefdx*dx;
      cpasma = lae.cpasma+dcpasma;
      a0 = cpasma/lae.cpasma*lae.aminor;
      r0 = rbdef + a0*vlim(1);
      z0 = zbdef + a0*vlim(2);
      v = linspace(0,2*pi,nbbbs_max)';
      if isfield(index_in_y,'rbbbs')
        y(index_in_y.rbbbs) = r0 + a0*cos(v);
      end
      if isfield(index_in_y,'zbbbs')
        y(index_in_y.zbbbs) = z0 + a0*sin(v);
      end
      if isfield(index_in_y,'nbbbs')
        y(index_in_y.nbbbs) = nbbbs_max;
      end
    end
  end
  equilibrium_update_is_complete = 1;
else
  equilibrium_update_is_complete = 0;
end

if isfield(index_in_y,'time')
  y(index_in_y.time) = time;
end
