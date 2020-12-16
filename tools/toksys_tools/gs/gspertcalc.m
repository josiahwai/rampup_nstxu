function [dpsizr,djphi,dbr,dbz,dpsizrerror] = gspertcalc(djphidpsizr,mps,rg,zg,mpp,method,options)
%
%  USAGE: [dpsizr,djphi,dbr,dbz,dpsizrerror] = gspertcalc(djphidpsizr,dpsizr_app,rg,zg,mpp,method,options)
%
%  PURPOSE: Calculate plasma response to applied flux for a given djphi/dpsizr
%
%  INPUTS:  djphidpsizr, how jphi responds to local flux perturbation [MA/m2/Wb]
%            dpsizr_app, perturbation(s) of applied flux [Wb] (size [nz,nr] or [nz*nr,:], e.g. mpc)
%                rg, zg, vectors of R, Z coordinates for the grid [m]
%                   mpp, (optional) supply (nz*nr,nr) mutuals [H] if available.
%                method, (optional, otherwise set automatically) calculation method (1, 2 or 3)
%                        1 = Exact solution optimized for small grids and few dpsizr_app
%                        2 = Exact solution optimized for larger grids and several dpsizr_app
%                        3 = Find solution for coarser grid and remove the error that comes from
%                            interpolating back to original grid by iterations that remove the
%                            error and associated plasma response. This method may fail if the
%                            coarse grid can't well represent the gradients in djphidpsizr
%                options, (with method 3), structure with optional fields:
%                         xr, rg for coarse grid is rg(1:xr:nr), defaults for xr, xz depend on gradients
%                         xz, zg for coarse grid is zg(1:xz:nz), in djphidpsizr and available memory
%                         maxerror (default = 1e-12), maximum acceptable value of dpsizrerror
%                         maxiter (default = 9), maximum iterations (even if dpsizrerror > maxerror)
%                         idoplot (default = 0), plot the remaining dpsizr errors after each iteration
%
%  OUTPUTS: dpsizr, perturbed flux [Wb] (= gscalc(djphi,rg,zg)+dpsizr_app)
%            djphi, perturbed current density [MA] (= djphidpsizr.*dpsizr)
%              dbr, perturbed radial field [T]
%              dbz, perturbed vertical field [T]
%      dpsizrerror, max(max(abs(gscalc(djphi,rg,zg)+dpsizr_app-dpsizr))), only recommended for method 3 
%
%  METHOD: Finite differences are solved for interior and mutuals are used for edge of grid

%	
%  VERSION @(#)gspertcalc.m	1.7 09/27/12
%
%  WRITTEN BY:  Anders Welander  ON	3/21/12
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mu0 = 4e-7*pi;

  nr = length(rg);
  nz = length(zg);
  ngg = nr*nz;
  rg = rg(:)';
  rgg = ones(nz,1)*rg;
  dr = (rg(nr)-rg(1))/(nr-1);
  dz = (zg(nz)-zg(1))/(nz-1);
  
  djphidpsizr = reshape(djphidpsizr,nz,nr); % Allow shape to be ngg*1 or 1*ngg, (reshape takes 36 microseconds on THOR)
    
  [n1,n2] = size(mps); % Will use n1,n2 to shape output same as mps a.k.a. dpsizr_app
  if n1*n2 == ngg % If true, we only have 1 applied flux
    mps = mps(:); % Collapse it to ngg*1
    nfluxes = 1;
  else
    nfluxes = n2;
  end
  
  if ~exist('mpp','var')
    mpp = [];
  end

  if size(mpp,1) ~= nz*nr | size(mpp,2) ~= nr
    calculate_mpp = 1;
    if ~isempty(mpp) % User supplied an mpp but it can't be used so warn
      disp('WARNING gspertcalc: Supplied mpp is wrong size. A new will be made.')
    end
    % Variables to be used with fine for calculation of mutual inductances
    zgg = zg*ones(1,nr);
    ri = rgg(:)-dr/2;
    ro = rgg(:)+dr/2;
    zl = zgg(:)-dz/2;
    zu = zgg(:)+dz/2;
    cur = zeros(ngg,1)+mu0;
    rp = zeros(ngg,1); % The rp value will be set before each call to fine
    zp = zeros(ngg,1)+zg(1);
  else
    calculate_mpp = 0;
  end  
  
  % Specify the calculation method
  if ~exist('method','var') | isempty(method)
    if ngg <=65*129
      method = 1;
    else
      method = 2;
    end
  end

  if method == 1 | method == 2
    if nargout > 4
      dpsizr_app = mps; % Needed to calculate dpsizr_err at the end of the code
    end
  end
  
  % Indices
  ii = (2:nz-1)'*ones(1,nr-2)+ones(nz-2,1)*(1:nr-2)*nz; ii = ii(:); % Indices to the interior of the grid
  iedges = [2:nz-1 ngg+(2-nz:-1) nz+1:nz:ngg-nz 2*nz:nz:ngg-nz]';
  icorners = [1 nz ngg-nz+1 ngg]';
  isides = [iedges; icorners];
  iplasma = find(djphidpsizr);
  jplasma = find(djphidpsizr(ii));
  if length(iplasma) > length(jplasma);
    djphidpsizr(1,:) = 0;
    djphidpsizr(nz,:) = 0;
    djphidpsizr(:,1) = 0;
    djphidpsizr(:,nr) = 0;
    iplasma = find(djphidpsizr);
  end
  np = length(iplasma);
  
  if method == 1 | method == 2  
    % The right-hand-side should be delstar for the interior, change mps(ii,j) to delstar(mps(ii,j))
    for j = 1:nfluxes
      mps(ii,j) = mps(ii,j)*(-2/dr^2-2/dz^2)+mps(ii-1,j)/dz^2+mps(ii+1,j)/dz^2+...
        mps(ii-nz,j).*(1/dr^2+1/2./rgg(ii)/dr)+mps(ii+nz,j).*(1/dr^2-1/2./rgg(ii)/dr);
    end
  end
  
  
  % All unknowns are changes of the total flux on the grid, dpsizr = dpsizr_app+dpsizr_pla
  % Equations are for flux change at the edges of the grid and change in delstar for interior
    
  if method == 1 % Method that uses a lot of memory
    
    if calculate_mpp
      mpp = zeros(ngg,nr);
      for k = 1:nr
        rp(:) = rg(k);
        [hr, hz, fl] = fine(ri,ro,zl,zu,cur,rp,zp);
        mpp(:,k) = fl;
      end
    end
    
    % Solving both finite difference equations and mutuals equations for the edge all at once
    T = sparse(ii,ii,   2*pi*1e6*mu0*rgg(ii).*djphidpsizr(ii)-2/dr^2-2/dz^2, ngg, ngg)+... % Equations for change in delstar: 
	sparse(ii,ii-nz,1/dr^2+1/2./rgg(ii)/dr,                              ngg, ngg)+... % delstar (dpsi_pla+dpsi_app) becomes:
	sparse(ii,ii+nz,1/dr^2-1/2./rgg(ii)/dr,                              ngg, ngg)+... %  (delstar+mu0*R*djdpsi) dpsi = delstar(dpsi_app)
	sparse(ii,ii-1, 1/dz^2,                                              ngg, ngg)+...
	sparse(ii,ii+1, 1/dz^2,                                              ngg, ngg)+...
	sparse(isides, isides, 1,                                            ngg, ngg)+... % The 1 in (1-mpp*di/dpsi)
	sparse(isides*ones(1,np), ones(length(isides),1)*iplasma', 1,        ngg, ngg);    % Allocate memory for edge equations

    % Edge equations for flux change, dpsi = dpsi_pla+dpsi_app becomes (1-mpp*di/dpsi)*dpsi = dpsi_app
    T(1:nz:ngg,iplasma) = -1e6*dr*dz*ones(nr,1)*djphidpsizr(iplasma)'.*mpp(iplasma,:)';
    T(nz:nz:ngg,iplasma) = -1e6*dr*dz*ones(nr,1)*djphidpsizr(iplasma)'.*mpp(iplasma+mod(nz-iplasma,nz)-mod(iplasma-1,nz),:)';
    for j = 2:nz-1
      izshift = (1+abs(j-(1:nz))-(1:nz))'*ones(1,nr);
      T(j,iplasma) = -1e6*dr*dz*djphidpsizr(iplasma)'.*mpp(iplasma+izshift(iplasma),1)';
      T(j+(nr-1)*nz,iplasma) = -1e6*dr*dz*djphidpsizr(iplasma)'.*mpp(iplasma+izshift(iplasma),end)';
    end

    % Here is where the magic happens
    dpsizr = T\mps;
  
  elseif method == 2 % Use less memory by solving FD first in terms of edges and then solve edges
    
    if calculate_mpp
      rp(:) = rg(1);
      [hr, hz, fl] = fine(ri,ro,zl,zu,cur,rp,zp);
      mp1 = fl;
      rp(:) = rg(nr);
      [hr, hz, fl] = fine(ri,ro,zl,zu,cur,rp,zp);    
      mpnr = fl;    
    end
    
    % Solving for interior in terms of edges
    mr = nr-2;
    mz = nz-2;
    mgg = mr*mz; % Number of grid points inside the edges
    jj = (1:mz)'*ones(1,mr)+ones(mz,1)*(0:mr-1)*mz; % Indices into (nz-2)*(nr-2) points of interior grid

    % Finite difference equations for interior points & only interior points as unknowns
    T = sparse(jj, jj,  2*pi*1e6*mu0*rgg(ii).*djphidpsizr(ii)-2/dr^2-2/dz^2,   mgg, mgg)+... % Equations for change in delstar: 
	sparse(jj(:,2:mr), jj(:,1:mr-1), 1/dr^2+1/2./rgg(ii(jj(:,2:mr)))/dr,   mgg, mgg)+... % delstar (dpsi_pla+dpsi_app) becomes:
	sparse(jj(:,1:mr-1), jj(:,2:mr), 1/dr^2-1/2./rgg(ii(jj(:,1:mr-1)))/dr, mgg, mgg)+... % (delstar+mu0*R*djdpsi) dpsi = delstar(dpsi_app)
	sparse(jj(2:mz,:), jj(1:mz-1,:), 1/dz^2,                               mgg, mgg)+...
	sparse(jj(1:mz-1,:), jj(2:mz,:), 1/dz^2,                               mgg, mgg);

    % The unknowns at the edges of the grid are on the right-hand-side in the finite difference equations
    R = sparse(1:mz,         1:mz,           -1/dr^2-1/2/rg(2)/dr,    mgg, 2*mz+2*mr)+... % Inner edge
	sparse(mgg+(1-mz:0), mz+(1:mz),      -1/dr^2+1/2/rg(nr-1)/dr, mgg, 2*mz+2*mr)+... % Outer edge
	sparse(1:mz:mgg,     2*mz+(1:mr),    -1/dz^2,                 mgg, 2*mz+2*mr)+... % Lower edge
	sparse(mz:mz:mgg,    2*mz+mr+(1:mr), -1/dz^2,                 mgg, 2*mz+2*mr);    % Upper edge

    dpsiinddelsmps = T\mps(ii,:); % Derivative of dpsizr(ii) w.r.t. delstar(mps(ii))
    dpsiindpsiedge = T\R;         % Derivative of dpsizr(ii) w.r.t. dpsizr(iedges)

    M = zeros(2*mz+2*mr); % Will hold equations for flux change at the edge in terms of flux changes at the edge
    C = zeros(2*mz+2*mr,nfluxes); % Small correction of edge flux if delstar(mps)~=0
    djphidpsiedge = djphidpsizr(iplasma)*ones(1,2*mz+2*mr).*dpsiindpsiedge(jplasma,:); % - delstar(dpsizr)/mu0/R
    djphiddelsmps = djphidpsizr(iplasma)*ones(1, nfluxes ).*dpsiinddelsmps(jplasma,:); % + delstar(dpsizr_app)/mu0/R
    for j = 2:nz-1
      izshift = (1+abs(j-(1:nz))-(1:nz))'*ones(1,nr);
      if calculate_mpp
	M(j-1,:) = -1e6*dr*dz*mp1(iplasma+izshift(iplasma))'*djphidpsiedge;
	C(j-1,:) = -1e6*dr*dz*mp1(iplasma+izshift(iplasma))'*djphiddelsmps;
	M(j-1+mz,:) = -1e6*dr*dz*mpnr(iplasma+izshift(iplasma))'*djphidpsiedge;
	C(j-1+mz,:) = -1e6*dr*dz*mpnr(iplasma+izshift(iplasma))'*djphiddelsmps;
      else
	M(j-1,:) = -1e6*dr*dz*mpp(iplasma+izshift(iplasma),1)'*djphidpsiedge;
	C(j-1,:) = -1e6*dr*dz*mpp(iplasma+izshift(iplasma),1)'*djphiddelsmps;
	M(j-1+mz,:) = -1e6*dr*dz*mpp(iplasma+izshift(iplasma),nr)'*djphidpsiedge;
	C(j-1+mz,:) = -1e6*dr*dz*mpp(iplasma+izshift(iplasma),nr)'*djphiddelsmps;
      end
    end
    izshift = (1+abs(nz-(1:nz))-(1:nz))'*ones(1,nr);
    for j = 2:nr-1
      if calculate_mpp
        rp(:) = rg(j);
        [hr, hz, fl] = fine(ri,ro,zl,zu,cur,rp,zp);
	M(j-1+2*mz,:) = -1e6*dr*dz*fl(iplasma)'*djphidpsiedge;
	C(j-1+2*mz,:) = -1e6*dr*dz*fl(iplasma)'*djphiddelsmps;
	M(j-1+2*mz+mr,:) = -1e6*dr*dz*fl(iplasma+izshift(iplasma))'*djphidpsiedge;
	C(j-1+2*mz+mr,:) = -1e6*dr*dz*fl(iplasma+izshift(iplasma))'*djphiddelsmps;
      else
	M(j-1+2*mz,:) = -1e6*dr*dz*mpp(iplasma,j)'*djphidpsiedge;
	C(j-1+2*mz,:) = -1e6*dr*dz*mpp(iplasma,j)'*djphiddelsmps;
	M(j-1+2*mz+mr,:) = -1e6*dr*dz*mpp(iplasma+izshift(iplasma),j)'*djphidpsiedge;
	C(j-1+2*mz+mr,:) = -1e6*dr*dz*mpp(iplasma+izshift(iplasma),j)'*djphiddelsmps;
      end
    end
    for j = 1:2*(nz-2)+2*(nr-2)
      M(j,j) = M(j,j)+1;
    end
    dpsie = inv(M)*(mps(iedges,:)-C);

    dpsizr = zeros(ngg,nfluxes);
    dpsizr(ii,:) = dpsiindpsiedge*dpsie+dpsiinddelsmps;
    dpsizr(iedges,:) = dpsie;

    % Finally do the corners
    djphi = zeros(ngg,nfluxes);
    djphi(iplasma,:) = djphidpsizr(iplasma)*ones(1,nfluxes).*dpsizr(iplasma,:);
    if calculate_mpp
      dpsizr(1,:) = mps(1,:)+mp1(iplasma)'*djphi(iplasma,:)*dr*dz*1e6;
      dpsizr(nz,:) = mps(nz,:)+mp1(iplasma+izshift(iplasma))'*djphi(iplasma,:)*dr*dz*1e6;
      dpsizr(ngg-nz+1,:) = mps(ngg-nz+1,:)+mpnr(iplasma)'*djphi(iplasma,:)*dr*dz*1e6;
      dpsizr(ngg,:) = mps(ngg,:)+mpnr(iplasma+izshift(iplasma))'*djphi(iplasma,:)*dr*dz*1e6;
    else
      dpsizr(1,:) = mps(1,:)+mpp(iplasma,1)'*djphi(iplasma,:)*dr*dz*1e6;
      dpsizr(nz,:) = mps(nz,:)+mpp(iplasma+izshift(iplasma),1)'*djphi(iplasma,:)*dr*dz*1e6;
      dpsizr(ngg-nz+1,:) = mps(ngg-nz+1,:)+mpp(iplasma,nr)'*djphi(iplasma,:)*dr*dz*1e6;
      dpsizr(ngg,:) = mps(ngg,:)+mpp(iplasma+izshift(iplasma),nr)'*djphi(iplasma,:)*dr*dz*1e6;
    end
    
  elseif method == 3 % Use even less memory by solving for a coarser grid and iterating to exact solution
    
    if ~exist('options','var')
      options.idoplot = 0;
    end
    if ~isstruct(options)
      if length(options) == 1
        options.idoplot = options;
      elseif length(options) == 2
        xr = options(1);
	options.xz = options(2);
	options.xr = xr;
      else
        options.idoplot = 0;
      end
    end
    
    if isfield(options,'maxerror')
      maxerror = options.maxerror;
    else
      maxerror = 1e-12;
    end
    
    if isfield(options,'maxiter')
      maxiter = options.maxiter;
    else
      maxiter = 9;
    end
    
    if isfield(options,'idoplot')
      idoplot = options.idoplot;
    else
      idoplot = 0;
    end
    
    % xr and xz are how much coarser the grid will be
    if isfield(options,'xr') % More clever defaults are claimed in the help section but not yet implemented!
      xr = options.xr;
    else
      xr = max(min(8,(nr-1)/16),1); % Never coarser than 17 and reduce by 8 at most
    end
    if isfield(options,'xz')
      xz = options.xz;
    else
      xz = max(min(8,(nz-1)/16),1); % Never coarser than 17 and reduce by 8 at most
    end
    ir = 1:xr:nr;
    iz = 1:xz:nz;
    nxr = length(ir);
    nxz = length(iz);
    igg = iz'*ones(1,nxr)+ones(nxz,1)*(ir-1)*nz;
    if ~exist('zgg','var')
      zgg = zg*ones(1,nr);
    end
    
    % Do recursive calls to this function
    if calculate_mpp
      [dpsizrc,djphic] = gspertcalc(djphidpsizr(iz,ir),mps(igg,:),rg(ir),zg(iz));
    else
      [dpsizrc,djphic] = gspertcalc(djphidpsizr(iz,ir),mps(igg,:),rg(ir),zg(iz),mpp(igg,ir));
    end
    
    dpsizr = zeros(ngg,nfluxes);
    djphi = zeros(ngg,nfluxes);
    for j = 1:nfluxes
      dpsizr_plac = reshape(dpsizrc(:,j)-mps(igg,j),nxz,nxr);
      dpsizr_pla2 = interp2(rg(ir),zg(iz),dpsizr_plac,rgg,zgg,'spline');
      dpsizr(:,j) = dpsizr_pla2(:)+mps(:,j);
      djphi(:,j) = djphidpsizr(:).*dpsizr(:,j);
    end
        
    dpsizrerror = inf;
    iteration = 0;
    
    while dpsizrerror > maxerror & iteration < maxiter
      iteration = iteration+1;
      if calculate_mpp
	dpsizr_err = dpsizr-gscalc(djphi,rg,zg)-mps;
	[dpsizrc,djphic] = gspertcalc(djphidpsizr(iz,ir),dpsizr_err(igg,:),rg(ir),zg(iz));
      else
	dpsizr_err = dpsizr-gscalc(djphi,rg,zg,mpp)-mps;
	[dpsizrc,djphic] = gspertcalc(djphidpsizr(iz,ir),dpsizr_err(igg,:),rg(ir),zg(iz),mpp(igg,ir));
      end
      for j = 1:nfluxes
        dpsizr_plac = reshape(dpsizrc(:,j)-dpsizr_err(igg,j),nxz,nxr);
        dpsizr_pla2 = interp2(rg(ir),zg(iz),dpsizr_plac,rgg,zgg,'spline');
	dpsizr(:,j) = dpsizr(:,j)-dpsizr_pla2(:)-dpsizr_err(:,j);
        djphi(:,j) = djphidpsizr(:).*dpsizr(:,j);
	if idoplot
	  plot(dpsizr_err(:,j))
	  title(['Errors in dpsizr for flux ' num2str(j) ', iteration ' num2str(iteration)])
	  a = axis; a(2) = size(dpsizr_err,1); axis(a);
	  try % Sometimes drawnow results in an error such as: "timeout waiting for window to show up"
	    drawnow
	  catch
	    drawproblem = lasterror;
	    disp(['gspertcalc: ' drawproblem.message ' Continuing the calculations...'])
	  end
        end
      end
      dpsizrerror = max(max(abs(dpsizr_err)));
    end
        
  end % End of choice of method
  
  
  % All is done. Pack it up.
  
  if nargout > 1 % return djphi
    if ~exist('djphi','var')
      djphi = djphidpsizr(:)*ones(1,nfluxes).*dpsizr; 
    end
    if nfluxes == 1
      djphi = reshape(djphi,n1,n2);
    end
  end
  
  if nargout > 2 % return dbr
    for j = 1:nfluxes
      y = reshape(dpsizr(:,j),nz,nr);
      a = [zeros(1,nr); (y(1:end-2,:)-y(3:end,:))/4/pi/dz./rgg(2:end-1,:); zeros(1,nr)];
      a(1,:) = 5*a(2,:)-7*a(3,:)+3*a(4,:); % Taylor expansion to order 2 around br(3,:)
      a(nz,:) = 5*a(nz-1,:)-7*a(nz-2,:)+3*a(nz-3,:);
      dbr(:,j) = a(:);
    end
    if nfluxes == 1
      dbr = reshape(dbr,n1,n2);
    end
  end

  if nargout > 3 % return dbz
    for j = 1:nfluxes
      y = reshape(dpsizr(:,j),nz,nr);
      a = [zeros(nz,1) (y(:,3:end)-y(:,1:end-2))/4/pi/dr./rgg(:,2:end-1) zeros(nz,1)];
      a(:,1) = 5*a(:,2)-7*a(:,3)+3*a(:,4); % Taylor expansion to order 2 around bz(:,3)
      a(:,nr) = 5*a(:,nr-1)-7*a(:,nr-2)+3*a(:,nr-3);
      dbz(:,j) = a(:);
    end
    if nfluxes == 1
      dbz = reshape(dbz,n1,n2);
    end
  end
 
  if nfluxes == 1
    dpsizr = reshape(dpsizr,n1,n2);
  end
  
  if nargout > 4 & method ~= 3 % return dpsizrerror (was already calculated if method == 3)
    if nfluxes == 1
      dpsizr_app = reshape(dpsizr_app,n1,n2);
    end
    if ~exist('dpsizrerror','var')
      if calculate_mpp
        dpsizrerror = max(max(abs(gscalc(djphi,rg,zg)+dpsizr_app-dpsizr)));
      else
        dpsizrerror = max(max(abs(gscalc(djphi,rg,zg,mpp)+dpsizr_app-dpsizr)));
      end
    end
  end
  
