%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_update_eq
%
%  PURPOSE: Update equilibrium: psizr, ic, iv, sp, sf and output y
%           Update response when dpsibar_limit_for_linear is reached
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  VERSION @(#)gs_update_eq.m	1.5 02/23/15
%
%  WRITTEN BY:  Anders Welander  ON	10/23/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The equilibrium is updated to new state vector

% Get dx needed from last-analyzed-equilibrium (lae) to reach new state
plasma_lae = any(lae.sp) | any(lae.sf);
gs_find_dx

% Determine if the change can be done all at once
if ~plasma_lae & plasma
  nlcriterion = inf;
elseif ~plasma
  nlcriterion = 0;
elseif circle % Using simple circular plasma model
  drmaxisx = drmaxisdx(1:nx-1)*dx(1:nx-1);
  dzmaxisx = dzmaxisdx(1:nx-1)*dx(1:nx-1);
  nlcriterion = 8*max(abs(drmaxisx)/dr,abs(dzmaxisx)/dz);
else
  dpsizrx(:) = dpsizrdx(:,1:nx-1)*dx(1:nx-1);
  dpsimagx = wa*dpsizrx(iia);
  dpsibryx = wb*dpsizrx(iib);
  azr = (1-lae.psibarzr)*dpsimagx;
  bzr = lae.psibarzr*dpsibryx;
  dpsibarzrx = (dpsizrx-azr-bzr)/(lae.psibry-lae.psimag);
  dpsibar = max(abs(dpsibarzrx(iplasma)));
  nlcriterion = dpsibar/dpsibar_limit_for_linear;
end
% f is the fraction of dx that can be done
f = min(1,1/nlcriterion);

if f == 1 % Simply update the outputs
  % Subtract flux error correction gradually, in full just before new analysis
  dx(nx) = -nlcriterion;
  % err is previous minus present equilibrium for state used to make present
  % Subtract err gradually 
  y  = lae.y  + dydx*dx  + (1-nlcriterion)*err.y;
  ys = lae.ys + dysdx*dx + (1-nlcriterion)*err.ys;
  return
end

% Here f < 1, change the equilibrium as much as allowed by nlcriterion
dx = f*dx;
dx(nx) = -1; % Time to subtract all of the flux error
% Remember previous equilibrium analyzed at state for new
previous.y  = lae.y  + dydx*dx;
previous.ys = lae.ys + dysdx*dx;
previous.li = lae.li + dlidx*dx;
previous.Wth = lae.Wth + dWthdx*dx;
previous.betap = lae.betap + dbetapdx*dx;
previous.cpasma = lae.cpasma + dcpasmadx*dx;

% Analyze the equilibrium obtained with approved changes to the state
psizr(:) = lae.psizr(:) + dpsizrdx*dx;
ic = lae.ic + dx(indic);
iv = lae.iv + dx(indiv);
sp = lae.sp + dx(indsp);
sf = lae.sf + dx(indsf);
rmaxis = lae.rmaxis + drmaxisdx*dx; % Should be updated to improve
zmaxis = lae.zmaxis + dzmaxisdx*dx; % success rate for gs_update_maxis
gs_eq_analysis
gs_response

% If the equilibrium is far from converged then give up on smooth outputs
if abs(dcpasmadx(nx)/cpasma) > dpsibar_limit_for_linear
  err.y(:) = 0;
  err.ys(:) = 0;
  gs_converge
  return
end

% If plasma quantities other than sp, sf are states then sp, sf
% need be corrected so this new equilibrium matches the states
if evolve_option == 1 | evolve_option == 99 % cpasma, li, betap are states
  dx = dxdxc(:,nc+nv+(1:3)) * [previous.cpasma-lae.cpasma; ...
    previous.li-lae.li; previous.betap-lae.betap];
  lae.y = lae.y + dydx*dx;
  lae.ys = lae.ys + dysdx*dx;
  lae.sp = lae.sp + dx(indsp);
  lae.sf = lae.sf + dx(indsf);
  lae.li = lae.li + dlidx*dx;
  lae.betap = lae.betap + dbetapdx*dx;
  lae.cpasma = lae.cpasma + dcpasmadx*dx;
  lae.psizr(:) = lae.psizr(:) + dpsizrdx*dx;
end
    
% Now archive error of previous equilibrium at new state
err.y  = previous.y  - lae.y ;
err.ys = previous.ys - lae.ys;


% Now try again to update the equilibrium

% Get dx needed from newly last-analyzed-equilibrium (lae) to reach new state
gs_find_dx

sp = lae.sp + dx(indsp);
sf = lae.sf + dx(indsf);
plasma = any(sp) | any(sf);

% Determine if the change can be done all at once
if ~plasma
  nlcriterion = 0;
elseif circle % Using simple circular plasma model
  drmaxisx = drmaxisdx(1:nx-1)*dx(1:nx-1);
  dzmaxisx = dzmaxisdx(1:nx-1)*dx(1:nx-1);
  nlcriterion = 8*max(abs(drmaxisx)/dr,abs(dzmaxisx)/dz);
else
  dpsizrx(:) = dpsizrdx(:,1:nx-1)*dx(1:nx-1);
  dpsimagx = wa*dpsizrx(iia);
  dpsibryx = wb*dpsizrx(iib);
  azr = (1-lae.psibarzr)*dpsimagx;
  bzr = lae.psibarzr*dpsibryx;
  dpsibarzrx = (dpsizrx-azr-bzr)/(lae.psibry-lae.psimag);
  dpsibar = max(abs(dpsibarzrx(iplasma)));
  nlcriterion = dpsibar/dpsibar_limit_for_linear;
end
% f is the fraction of dx that can be done
f = min(1,1/nlcriterion);

if f == 1 % Simply update the outputs
  % Subtract flux error correction gradually, in full just before new analysis
  dx(nx) = -nlcriterion;
  % err is previous minus present equilibrium for state used to make present
  % Subtract err gradually 
  y  = lae.y  + dydx*dx  + (1-nlcriterion)*err.y;
  ys = lae.ys + dysdx*dx + (1-nlcriterion)*err.ys;
  return
end


% At this point a new equilibrium analysis failed to get close enough
% to the desired new state vector. Give up on smoothness of outputs
err.y(:) = 0;
err.ys(:) = 0;
gs_converge

return










% Loop until targets are reached
done = false;
iteration_counter = 0;

while ~done

  gs_find_dx
  g = 0;
  
  % Total changes due to inputs since last analyzed equilibrium
  dpsizrx(:) = dpsizrdx(:,1:nx-1)*dx(1:nx-1);
  dpsimagx = wa*dpsizrx(iia);
  dpsibryx = wb*dpsizrx(iib);

  if nbbbs > 1
  
    % Decide how much change dx can be done before running gs_eq_analysis

    azr = (1-lae.psibarzr)*dpsimagx;
    bzr = lae.psibarzr*dpsibryx;
    dpsibarzrx = (dpsizrx-azr-bzr)/(lae.psibry-lae.psimag);

    dpsibar = max(abs(dpsibarzrx(iplasma)));

    dx(nx) = -dpsibar/dpsibar_limit_for_linear;
    % All of the flux error effect is now already subtracted from y, ys
    dx(nx) = 0;
    % Amount of the error in y according to last gs_eq_analysis to keep
    g = 1-min(1,dpsibar/dpsibar_limit_for_linear);

    % The amount (f) of dpsizrx that will be applied is from 0 to 100%
    f = min(1,dpsibar_limit_for_linear/dpsibar);

  else % No plasma
    f = 1;
  end 
    
  dpsizr(:) = f*(dpsizrx(:) + dpsizrdx(:,nx)*dx(nx));
   
  psizr = lae.psizr + dpsizr;
  ic = lae.ic + f*dx(indic);
  iv = lae.iv + f*dx(indiv);
  sp = lae.sp + f*dx(indsp);
  sf = lae.sf + f*dx(indsf);
    
  done = f == 1 | iteration_counter > 1;
  iteration_counter = iteration_counter+1;
  
  % Update the outputs
  y  = lae.y + dydx*dx + err.y*g;
  ys  = lae.ys + dysdx*dx + err.ys*g;

  if f < 1
  
    % The equilibrium analysis will reveal errors in these quantities
    if iteration_counter == 1
      previous.y  = lae.y + dydx*dx*f;
      previous.ys = lae.ys + dysdx*dx*f;
      previous.cpasma = lae.cpasma + dcpasmadx*dx*f;
      previous.li = lae.li + dlidx*dx*f;
      previous.betap = lae.betap + dbetapdx*dx*f;
    end
  
    if f < 1
      gs_converge
      er = 0;
    else
      gs_eq_analysis
      gs_response
      % Subtract all of the flux error right away from y, cpasma, li, betap
      lae.y = lae.y - dydx(:,nx);
      lae.ys = lae.ys - dysdx(:,nx);
      lae.cpasma = lae.cpasma - dcpasmadx(1,nx);
      lae.li = lae.li - dlidx(1,nx);
      lae.betap = lae.betap - dbetapdx(1,nx);
      lae.psizr(:) = lae.psizr(:) - dpsizrdx(:,nx);
      er = 0;
    end
    
    if iteration_counter == 1
      if evolve_option == 1 || evolve_option == 99
	% If cpasma, li, betap are states in xs then these values
	% can not change for the new equilibrium
	% Correct sp, sf, to preserve cpasma, li, betap
        dx = dxdxc(:,nc+nv+(1:3)) * [previous.cpasma-lae.cpasma; ...
          previous.li-lae.li; previous.betap-lae.betap];
	lae.sp = lae.sp + dx(indsp);
	lae.sf = lae.sf + dx(indsf);
        lae.cpasma = lae.cpasma + dcpasmadx*dx;
        lae.li = lae.li + dlidx*dx;
        lae.betap = lae.betap + dbetapdx*dx;
        lae.y = lae.y + dydx*dx;
        lae.ys = lae.ys + dysdx*dx;
        lae.psizr(:) = lae.psizr(:) + dpsizrdx*dx;
      end
    end
        
    % Calculate the errors to remove them gradually rather than all at once
    err.y  = previous.y  - lae.y ;
    err.ys = previous.ys - lae.ys;

  end
      
end % End of while loop with condition ~done

if evolve_option == 1 || evolve_option == 99
  if isfield(index_in_y,'cpasma')
    err.y(index_in_y.cpasma) = 0;
  end
  if isfield(index_in_y,'li')
    err.y(index_in_y.li) = 0;
  end
  if isfield(index_in_y,'betap')
    err.y(index_in_y.betap) = 0;
  end
end

gs_find_dx
% Update the outputs
y  = lae.y + dydx*dx;
ys  = lae.ys + dysdx*dx;
