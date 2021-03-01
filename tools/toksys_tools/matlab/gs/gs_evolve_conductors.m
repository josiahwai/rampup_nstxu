%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_evolve_conductors
%
%  PURPOSE: Calculate the time derivative of conductor currents
%
%  INPUTS:  new_response_was_calculated, flag to update Amat, Bmat
%           constraints, flag for constraints on plasma profiles
%           lae.ys, flux at conductors from plasma for last analyzed equilibrium
%           dysdx, response of flux at conductors from plasma
%           dx, change of state vector since last analyzed equilibrium
%           f, fraction of dx to apply (normally 1)
%           mss, mutual inductances between conductors when no plasma
%           ic, iv, coil and vessel current vectors
%           u1, voltages applied to conductors
%
%  OUTPUTS: x1, state of conductor system
%           x1dot, time derivative of x1
%           v1, voltage induced at conductors due to x1dot
%           y1, total flux at conductors including flux error corrections
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	11/15/13
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if new_response_was_calculated
    
  % Inductance matrix including plasma response
  lstar = dysdx*dxdxc(:,1:nc+nv);
  
  % Taking the inverse to solve for xdot
  lstari = inv(lstar);
  
  % x1dot = Amat*x1 + Bmat*u1
  A1 = -lstari*rss;
  B1 = lstari;
  
end

% State vector for the conductor system
x1 = [ic; iv];

% Time derivative of these states
x1dot = A1*x1 + B1*u1;

% Resistive voltage on conductors
vres = ress.*x1;
