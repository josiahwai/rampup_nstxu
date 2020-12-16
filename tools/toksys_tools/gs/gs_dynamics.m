%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_dynamics
%
%  PURPOSE: Calculate the time derivative xsdot
%
%  INPUTS:  new_response_was_calculated, flag to update Amat, Bmat
%           evolve_option, flag for how to evolve xs
%           dysdx, response of flux at conductors
%           dpsipladx, response of plasma flux
%
%  OUTPUTS: xsdot, time derivative of xs
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	4/2/14
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The different values of evolve_option:
% = 1 Scalar Ip with li, betap (or Wth) disturbances
% = 2 1-D with resistivity and net-heating profiles
% = 3 2-D with resistivity that varies on flux surfaces
% = 9 Like 1 except Ip is also evolved externally, same as  li, betap/Wth

if evolve_option == 1
  if time < plares.times(1)
    Rpla = plares.values(1);
  elseif time > plares.times(end)
    Rpla = plares.values(end);
  else
    Rpla = interp1(plares.times,plares.values,time);
  end
  if Rpla ~= Rpla_previous | new_response_was_calculated
    calculate_AB = true;
  else
    calculate_AB = false;    
  end
  Rpla_previous = Rpla;
else
  calculate_AB = new_response_was_calculated;
end
if plasma & calculate_helical_voltage
  gs_trace_contours
  gs_contour_response
  gs_contour_profiles
  gs_profile_response
  gs_helical_voltage
end
if calculate_AB | isempty(Amat)  
  if evolve_option == 1 % Mimic scalar Ip evolution done with linear models    
    if plasma
      Rhat(nci+nv+1,nci+nv+1) = Rpla;
      % Inductances for conductors and plasma before circuits are connected
      L_star = [dysdx; dpsipladx; zeros(2,nx)]*dxdxc(:,1:nc+nv+3);
      L_star(end-1,end-1) = 1; % The equation dli/dt = dlit/dt
      L_star(end-0,end-0) = 1; % The equation dbetap/dt = dbetapt/dt
    else % No plasma    
      L_star = [dysdx(:,1:nc+nv+3); zeros(3,nc+nv+3)];
      L_star(end-0,end-0) = 1;
      L_star(end-1,end-1) = 1;
      L_star(end-2,end-2) = 1;
    end
    
    lstar = Pxx'*(L_star+Lextra)*Pxx; % Circuits connected  
    lstari = inv(lstar);
    Amat = -lstari*Rhat;
    Bmat = lstari*Vhat;
    
  elseif evolve_option == 2 % 1D or 2D resistivity
    if plasma
      
      % Pj projects ncont voltage equations on nkn+2 equations
      lai = inv([dysdx(:,1:nc+nv+nkn+2); zeros(nkn+2,nc+nv) eye(nkn+2)]);
      dxdsf = [-lai*[dysdx(:,indsf); zeros(nkn+2)]; eye(nkn+2)];
      Pj = [pinv(dvindcdxdot(1:ncont-1,1:nx-1)*dxdsf) zeros(nkn+2,1)];
      
      L_star = [dysdx(:,1:nx-1); ...                      % Conductor eqs
        [zeros(nkn+2,nc+nv) eye(nkn+2) zeros(nkn+2)]; ... % pressure eqs
        Pj*dvindcdxdot(:,1:nx-1)];                        % current eqs
      
      % Parallel resistance not in Amat for following reasons:
      % 1. Exact same xsdot is obtained by -Bmat*vresc
      % 2. No need to calculate offset voltage dvrescFIXEDdxs*xs
      % 3. gamma for Amat depends only weakly on this resistance
      % 4. jbs done outside so its effect on gamma is missing anyway
      % rxx(indsf,:) = Pj*dvrescFIXEDdx(:,1:nx-1); % SO SKIPPING THIS
      % Equations for j:  V_CD-vresc = Lc*xsdot
    else % No plasma
      L_star = [dysdx(:,1:nx-1); [zeros(2*nkn+4,nc+nv) eye(2*nkn+4)]];
      Pj = zeros(nkn+2,ncont);
    end
    
    lstar = Pxx'*(L_star+Lextra)*Pxx; % Circuits connected
    lstari = inv(lstar);
    Amat = -lstari*(Pxx'*(rxx+Rextra)*Pxx);
    Bmat = lstari*[eye(nci+nv+nkn+2,nci+nv+nkn+2+ncont); ...
                   zeros(nkn+2,nci+nv+nkn+2) Pj];
    
  elseif evolve_option == 9
  % Evolve ic, iv, input time-derivatives of current and pressure
    if constraints
      % profiles constrained to a cpasma, li, betap or Wth   
      if plasma
	L_star = [dysdx; zeros(3,nx)]*dxdxc(:,1:nc+nv+3);
      else
	L_star = [dysdx(:,1:nc+nv+3); zeros(3,nc+nv+3)];
      end
      L_star(end-2,end-2) = 1; % The equation dcpasma/dt = dcpasma/dt
      L_star(end-1,end-1) = 1; % The equation dli/dt = dlit/dt
      L_star(end-0,end-0) = 1; % The equation dbetap/dt = dbetapt/dt
    else
      L_star = [dysdx; [zeros(2*nkn+4,nc+nv) eye(2*nkn+4) zeros(2*nkn+4,1)]];
    end
    % Inductances with circuits connected
    lstar = Pxx'*(L_star+Lextra)*Pxx;

    lstari = inv(lstar);
    Amat = -lstari*Rhat;
    Bmat = lstari;
  else
    error('gs_dynamics: Invalid evolve_option')
  end

  stabilized = false;
  if stabilize | ...
     isfield(index_in_y,'gamma') | ...
     isfield(index_in_y,'drcurdv') | ...
     isfield(index_in_y,'dzcurdv')
    [vecs,vals] = eigsort(Amat);            
    gamma = vals(1);
    if plasma
      if constraints
        drcurdv = drcurdx*dxdxc(:,1:nxc-1)*Pxx*vecs(:,1);
        dzcurdv = dzcurdx*dxdxc(:,1:nxc-1)*Pxx*vecs(:,1);
      else
        drcurdv = drcurdx(1:nx-1)*Pxx*vecs(:,1);
        dzcurdv = dzcurdx(1:nx-1)*Pxx*vecs(:,1);
      end
    else
      drcurdv = NaN;
      dzcurdv = NaN;
    end
    s = sign(dzcurdv);
    if isfield(index_in_y,'gamma')
      y(index_in_y.gamma) = gamma;
      lae.y(index_in_y.gamma) = gamma;
    end
    if isfield(index_in_y,'drcurdv')
      y(index_in_y.drcurdv) = s*drcurdv;
      lae.y(index_in_y.drcurdv) = s*drcurdv;
    end
    if isfield(index_in_y,'dzcurdv')
      y(index_in_y.dzcurdv) = s*dzcurdv;
      lae.y(index_in_y.dzcurdv) = s*dzcurdv;
    end
    if stabilize
      if vals(1) > 0 && all(isreal(vecs(:,1)))	
%         Amat = Amat - Amat*vecs(:,1)*vecs(:,1)';
        Amat = real(negateFirstEigval(Amat));
%         Amat = negateFirstEigval(Amat);        
        stabilized = true;
      end
    end
  end    
  if isfield(index_in_y,'dgapdpdvps') & plasma
    if constraints
      lae.y(index_in_y.dgapdpdvps) = ...
        dgapdpdx*dxdxc(:,1:end-1)*Pxx*Bmat(:,1:nps);  
    else
      lae.y(index_in_y.dgapdpdvps) = ...
        dgapdpdx(:,1:nx-1)*Pxx*Bmat(:,1:nps);  
    end
  end
  if isfield(index_in_y,'dgapdpdvcd') & plasma
    if constraints
      lae.y(index_in_y.dgapdpdvcd) = ...
        dgapdpdx*dxdxc(:,1:end-1)*Pxx*Bmat(:,nci+nv+(1:ncd));  
    else
      lae.y(index_in_y.dgapdpdvcd) = ...
        dgapdpdx(:,1:nx-1)*Pxx*Bmat(:,nci+nv+(1:ncd));
    end
  end
  if isfield(index_in_y,'dgapdpdxs') & plasma
    if constraints
      lae.y(index_in_y.dgapdpdxs) = ...
        dgapdpdx*dxdxc(:,1:end-1)*Pxx;  
    else
      lae.y(index_in_y.dgapdpdxs) = ...
        dgapdpdx(:,1:nx-1)*Pxx;
    end
  end
  if isfield(index_in_y,'dzcurdotdvps') & plasma
    if constraints
      lae.y(index_in_y.dzcurdotdvps) = ...
        dzcurdx*dxdxc(:,1:end-1)*Pxx*Bmat(:,1:nps);  
    else
      lae.y(index_in_y.dzcurdotdvps) = ...
        dzcurdx(1,1:nx-1)*Pxx*Bmat(:,1:nps);  
    end
  end
  if isfield(index_in_y,'dzcurdotdvcd') & plasma
    if constraints
      lae.y(index_in_y.dzcurdotdvcd) = ...
        dzcurdx*dxdxc(:,1:end-1)*Pxx*Bmat(:,nci+nv+(1:ncd));  
    else
      lae.y(index_in_y.dzcurdotdvcd) = ...
        dzcurdx(1,1:nx-1)*Pxx*Bmat(:,nci+nv+(1:ncd));
    end
  end
  if isfield(index_in_y,'dcpasmadotdvps') & plasma
    if constraints
      lae.y(index_in_y.dcpasmadotdvps) = ...
        dcpasmadx*dxdxc(:,1:end-1)*Pxx*Bmat(:,1:nps);  
    else
      lae.y(index_in_y.dcpasmadotdvps) = ...
        dcpasmadx(1,1:nx-1)*Pxx*Bmat(:,1:nps);  
    end
  end
  if isfield(index_in_y,'dcpasmadotdvcd') & plasma
    if constraints
      lae.y(index_in_y.dcpasmadotdvcd) = ...
        dzcurdx*dxdxc(:,1:end-1)*Pxx*Bmat(:,nci+nv+(1:ncd));  
    else
      lae.y(index_in_y.dcpasmadotdvcd) = ...
        dcpasmadx(1,1:nx-1)*Pxx*Bmat(:,nci+nv+(1:ncd));
    end
  end
end

if isempty(xs) % Need to return an (initial) state vector xs0
  if evolve_option == 1
    if constraints == 1
      xs0 = [picci*ic; iv; cpasma; li; betap];
    elseif constraints == 2
      xs0 = [picci*ic; iv; cpasma; li; Wth];
    end
  elseif evolve_option == 2
    xs0 = [picci*ic; iv; sp; sf];
  elseif evolve_option == 9
    if constraints == 1
      xs0 = [picci*ic; iv; cpasma; li; betap];
    elseif constraints == 2
      xs0 = [picci*ic; iv; cpasma; li; Wth];
    end
  else  
    error('gs_dynamics: Invalid evolve_option')
  end
end
