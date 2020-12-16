%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_xsdot
%
%  PURPOSE: Calculate xsdot
%           d
%  INPUTS:  Amat, Bmat, xs, u
%	
%  OUTPUTS:  xsdot, time derivative of state vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  WRITTEN BY:  Anders Welander ON 2015-03-08
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate time derivative of xs

xsdotA = Amat*xs;
if evolve_option == 2 & plasma
  % The resistivity can be handled both internally (here)
  % and/or outside by calculating voltages applied to each surface
  % Resistivity handled here can come in either once in config.eta_vs_rhot or
  % each time in u.eta_vs_rhot or u.etacont
  if exist('u','var') & isfield(u,'etacont') % User sent eta for contours
    etacont = u.etacont;
    etac = (etacont(1:npola-1,:)+etacont(2:npola,:))/2;
  elseif exist('u','var') &  isfield(u,'eta_vs_rhot')
    eta_vs_rhot = u.eta_vs_rhot;
    rhoe = linspace(0,1,length(eta_vs_rhot))';
    eta = spline(rhoe, eta_vs_rhot, rhot);
    etacont = ones(npola,1)*eta';
    etac = (etacont(1:npola-1,:)+etacont(2:npola,:))/2;
  else
    etac = lae.etac; % Resistivity between contour points
  end
  
  % Parallel surface averaged voltages
  vresc0 = sum(dvrescFdetac.*etac)'; % Voltage with new etac but lae currents
  vresc = vresc0 + dvrescFdx*dx; % Add change due to new currents, dx
  if isfield(index_in_y,'vres')
    y(index_in_y.vres) = vresc;
  end
  
  xsdotA = xsdotA - Bmat(:,nci+nv+2+nkn+(1:ncont))*vresc;
end

if isfield(u,'vps')
  nvps = length(u.vps(:));
  xsdotB1 = Bmat(:,1:nvps)*u.vps(:);
else
  xsdotB1 = zeros(nxs,1);
end

if isfield(u,'vcd')
  nvcd = length(u.vcd(:));
  xsdotB2 = Bmat(:,nci+nv+(1:nvcd))*u.vcd(:);
else
  xsdotB2 = zeros(nxs,1);
end
if isfield(index_in_y,'xsdotB2')
  y(index_in_y.xsdotB2) = xsdotB2;
end
xsdotB = xsdotB1 + xsdotB2;

xsdot = xsdotA + xsdotB;

if isfield(index_in_y,'vplas')
  if constraints
    y(index_in_y.vplas) = dydx(index_in_y.vplas,:)*dxdxc*...
      [Pxx; zeros(1,nxs)]*[xsdot(1:nci+nv); 0; xsdot(nci+nv+2:nxs)];
  else
    y(index_in_y.vplas) = dydx(index_in_y.vplas,:)*...
      [Pxx; zeros(1,nxs)]*[xsdot(1:nci+nv); zeros(nxs-nci-nv,1)];
  end
end

if isfield(index_in_y,'gapdpdot')
  if constraints
    y(index_in_y.gapdpdot) = dgapdpdx*dxdxc*[Pxx*xsdot; 0];
  else
    y(index_in_y.gapdpdot) = dgapdpdx*[Pxx*xsdot; 0];
  end
end

if isfield(index_in_y,'gapdpdotA')
  if plasma
    if constraints
      y(index_in_y.gapdpdotA) = dgapdpdx*dxdxc*[Pxx*xsdotA; 0];
    else
      y(index_in_y.gapdpdotA) = dgapdpdx*[Pxx*xsdotA; 0];
    end
  else
    y(index_in_y.gapdpdotA) = nan;
  end
end

if isfield(index_in_y,'gapdpdotB')
  if constraints
    y(index_in_y.gapdpdotB) = dgapdpdx*dxdxc*[Pxx*xsdotB; 0];
  else
    y(index_in_y.gapdpdotB) = dgapdpdx*[Pxx*xsdotB; 0];
  end
end

if isfield(index_in_y,'gapdpdotB2')
  if constraints
    y(index_in_y.gapdpdotB2) = dgapdpdx*dxdxc*[Pxx*xsdotB2; 0];
  else
    y(index_in_y.gapdpdotB2) = dgapdpdx*[Pxx*xsdotB2; 0];
  end
end

if isfield(index_in_y,'zcurdot')
  if constraints
    y(index_in_y.zcurdot) = dzcurdx*dxdxc*[Pxx*xsdot; 0];
  else
    y(index_in_y.zcurdot) = dzcurdx*[Pxx*xsdot; 0];
  end
end

if isfield(index_in_y,'zcurdotA')
  if plasma
    if constraints
      y(index_in_y.zcurdotA) = dzcurdx*dxdxc*[Pxx*xsdotA; 0];
    else
      y(index_in_y.zcurdotA) = dzcurdx*[Pxx*xsdotA; 0];
    end
  else
    y(index_in_y.zcurdotA) = nan;
  end
end

if isfield(index_in_y,'zcurdotB')
  if constraints
    y(index_in_y.zcurdotB) = dzcurdx*dxdxc*[Pxx*xsdotB; 0];
  else
    y(index_in_y.zcurdotB) = dzcurdx*[Pxx*xsdotB; 0];
  end
end

if isfield(index_in_y,'zcurdotB2')
  if constraints
    y(index_in_y.zcurdotB2) = dzcurdx*dxdxc*[Pxx*xsdotB2; 0];
  else
    y(index_in_y.zcurdotB2) = dzcurdx*[Pxx*xsdotB2; 0];
  end
end

if isfield(index_in_y,'cpasmadot')
  if constraints
    y(index_in_y.cpasmadot) = dcpasmadx*dxdxc*[Pxx*xsdot; 0];
  else
    y(index_in_y.cpasmadot) = dcpasmadx*[Pxx*xsdot; 0];
  end
end

if isfield(index_in_y,'cpasmadotA')
  if plasma
    if constraints
      y(index_in_y.cpasmadotA) = dcpasmadx*dxdxc*[Pxx*xsdotA; 0];
    else
      y(index_in_y.cpasmadotA) = dcpasmadx*[Pxx*xsdotA; 0];
    end
  else
    y(index_in_y.cpasmadotA) = nan;
  end
end

if isfield(index_in_y,'cpasmadotB')
  if constraints
    y(index_in_y.cpasmadotB) = dcpasmadx*dxdxc*[Pxx*xsdotB; 0];
  else
    y(index_in_y.cpasmadotB) = dcpasmadx*[Pxx*xsdotB; 0];
  end
end

if isfield(index_in_y,'cpasmadotB2')
  if constraints
    y(index_in_y.cpasmadotB2) = dcpasmadx*dxdxc*[Pxx*xsdotB2; 0];
  else
    y(index_in_y.cpasmadotB2) = dcpasmadx*[Pxx*xsdotB2; 0];
  end
end

if isfield(index_in_y,'vind')
  if constraints
    lae.y(index_in_y.vind) = dvindcdxdot*dxdxc*[Pxx*xsdot; 0];
  else
    lae.y(index_in_y.vind) = dvindcdxdot(:,1:nx-1)*Pxx*xsdot;
%    psipladot = dydx(index_in_y.psipla,1:end-1)*Pxx*xsdot
  end
end
