%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_converge
%
%  PURPOSE: Iterate to a small flux error
%
%  INPUTS: The work space for gs codes and:
%          min_iterations (default = 1)
%          max_iterations (default = 9)
%
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	10/23/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('min_iterations','var')
  min_iterations = 1;
end
if ~exist('max_iterations','var')
  max_iterations = 19;
end

% Counter for number of iterations
iteration_counter = 0;

done = false;

while ~done

  gs_eq_analysis
  gs_response
  gs_find_dx
  
  dpsizrx(iplasma) = dpsizrdx(iplasma,1:nx-1)*dx(1:nx-1);
  dpsimagx = wa*dpsizrx(iia);
  dpsibryx = wb*dpsizrx(iib);
  azr = (1-lae.psibarzr)*dpsimagx;
  bzr = lae.psibarzr*dpsibryx;
  dpsibarzrx = (dpsizrx-azr-bzr)/(lae.psibry-lae.psimag);
  dpsibar = max(abs(dpsibarzrx(iplasma)));
  
  costfun = max(abs(psizr_err(:)/(psibry-psimag)));
  if evolve_option == 0 & exist('xc','var') & isfield(xc,'cpasma')
    if constraints == 1 | constraints == 2
      costfun = max(costfun,abs(xc.cpasma-cpasma)/abs(xc.cpasma));
    end
  end
  
  dpsizr(:) = dpsizrdx*dx;
  drbbbs = drbdx*dx;
  dzbbbs = dzbdx*dx;
  drmaxis = drmaxisdx*dx;
  dzmaxis = dzmaxisdx*dx;

  g = max(max(abs(drbbbs(1:nbbbs)/dr)), max(abs(dzbbbs(1:nbbbs)/dz)));
  g = max(g,abs(drmaxis/dr));
  g = max(g,abs(dzmaxis/dz));
  %f = min([1,  0.2/max(abs(dpsizr(:)/(psibry-psimag))),  2/g]);
  f = min([1,  0.05/max(abs(dpsizr(:)/(psibry-psimag))),  0.5/g]);
  f = min([1,  0.1/max(abs(dpsizr(:)/(psibry-psimag))),  0.25/g]);
  
  % Deal with non-linear change of current
  fd = psibry-psimag;
  dfd = dpsibrydx*dx-dpsimagdx*dx;
  if 0*abs(dfd*f/fd) > 0.2
    fd_app = wb*psizr_app(iib)-wa*psizr_app(iia);
    fd_pla = wb*psizr_pla(iib)-wa*psizr_pla(iia);
    dfdis = dpsibrydx(indis)*dx(indis)-dpsimagdx(indis)*dx(indis);
    dfdsp = dpsibrydx(indsp)*dx(indsp)-dpsimagdx(indsp)*dx(indsp);
    dfdsf = dpsibrydx(indsf)*dx(indsf)-dpsimagdx(indsf)*dx(indsf);
    cpasma_target = lae.cpasma+dcpasmadx*dx; % Where we want to go
    % The current scales with [sp sf]/(psibry-psimag)
    % How change sp, sf if we assume bdef and maxis don't move?
    % cpasma/cpasma_target = [sp sf]/fd*(fd_app+dfdis+fd_pla/h)/(x*[sp sf])
    h = cpasma/cpasma_target;
    % x*h = (fd_app+dfdis+h*fd_pla)/fd
    x = (fd_app+dfdis+fd_pla/h)/fd/h;
    dx(indsp) = x*sp;
    dx(indsf) = x*sf;
    psizr(:) = psizr_app(:) + dpsizrdx(:,indis)*dx(indis) + psizr_pla(:)/h;
    ic = ic+dx(indic);
    iv = iv+dx(indiv);
    sp = x*sp;
    sf = x*sf;
    er = 0;
  else  
    psizr = psizr + f*dpsizr;
    ic    = ic    + f*dx(indic);
    iv    = iv    + f*dx(indiv);
    sp    = sp    + f*dx(indsp);
    sf    = sf    + f*dx(indsf);
    er    = er    + f*dx(inder);
  end
  
  % Plot convergence progress if desired
  
  iteration_counter = iteration_counter+1;
  if plotit > 1
    gs_plot_progress
  end
  drawnow % For simulink
  
  done = costfun < dpsibar_limit_for_linear | ...
    iteration_counter == max_iterations;
  done = done & iteration_counter >= min_iterations;

end

% Plot end result if desired
if plotit == 1
  gs_plot_progress
end

% This is the flux at conductors for this equilibrium
ys = lae.ys + dysdx*dx;

% This is the output for this equilibrium
y = lae.y + dydx*dx;
