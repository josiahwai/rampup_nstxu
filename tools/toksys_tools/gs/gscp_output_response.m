%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gscp_output_response
%
%  PURPOSE: Populate dydx with contents controlled by index_in_y
%
%  INPUTS: index_in_y and workspace after gscp_analysis_response
%
%  OUTPUTS:  dydx
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander ON 2015-03-04
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(index_in_y,'psizr')
  dydx(index_in_y.psizr,:) = dpsizrdx;
end
if isfield(index_in_y,'fl')
  dfldx = [mlc mlv]*eye(nc+nv,nx)+mpl'*dpcurrtdx;
  dydx(index_in_y.fl,:) = dfldx;
end
if isfield(index_in_y,'lv')
  dlvdx = [mhc mhv]*eye(nc+nv,nx)+mph'*dpcurrtdx;
  dydx(index_in_y.lv,:) = dlvdx;
end
if isfield(index_in_y,'bp')
  dbpdx = [gbc gbv]*eye(nc+nv,nx)+gpb'*dpcurrtdx;
  dydx(index_in_y.bp,:) = dbpdx;
end
if isfield(index_in_y,'rog')
  drogdx = rldata(:,1:nc+nv)*eye(nc+nv,nx)+rldata(:,nc+nv+1)*dcpasmadx;
  dydx(index_in_y.rog,:) = drogdx;
end
if isfield(index_in_y,'psidp')
  dpsidpdx = [mdc mdv]*eye(nc+nv,nx)+mpd'*dpcurrtdx;
  dydx(index_in_y.psidp,:) = dpsidpdx;
end
if isfield(index_in_y,'brdp')
  dbrdpdx = [grdc grdv]*eye(nc+nv,nx)+grpd'*dpcurrtdx;
  dydx(index_in_y.brdp,:) = dbrdpdx;
end
if isfield(index_in_y,'bzdp')
  dbzdpdx = [gzdc gzdv]*eye(nc+nv,nx)+gzpd'*dpcurrtdx;
  dydx(index_in_y.bzdp,:) = dbzdpdx;
end
if isfield(index_in_y,'btdp')
  dbtdpdx = zeros(ndp,nx);  % MUST ADD CODE HERE!
  dydx(index_in_y.btdp,:) = dbtdpdx;
end
if plasma & (...
    isfield(index_in_y,'gapdp') | ...
    isfield(index_in_y,'gapdpdot') | ...
    isfield(index_in_y,'gapdpdotA') | ...
    isfield(index_in_y,'gapdpdotB') | ...
    isfield(index_in_y,'dgapdpdvps'))
  dr0dx = dpdx(nrh+1,:);
  dz0dx = dpdx(nrh+2,:);
  da0dx = dpdx(nrh+3,:);
  for i = 1:ndp
    dgapdpdx(i,:) = -((rdp(i)-r0)*dr0dx+(zdp(i)-z0)*dz0dx)/rhodp(i) - da0dx;
  end
end
if isfield(index_in_y,'gapdp') & plasma
  dydx(index_in_y.gapdp,:) = dgapdpdx;
end
if isfield(index_in_y,'gaps')
  dgapsdx = zeros(ngaps,nx);  % MUST ADD CODE HERE!
  dydx(index_in_y.gaps,:) = dgapsdx;
end
if isfield(index_in_y,'ys')
  dydx(index_in_y.ys,:) = dysdx;
end
if isfield(index_in_y,'rx')
  drxoutputdx = zeros(nxpoints,nx); % nan's can't be used in Simulink so using zeros
  n = min(nxpoints,nulls.count);
  for i = 1:n
    ii = nulls.k(kxpoints(i))+neighbors;
    drxoutputdx(i,:) = nulls.drdpsi(kxpoints(i),:)*dpsizrdx(ii,:);
  end
  dydx(index_in_y.rx,:) = drxoutputdx;
end
if isfield(index_in_y,'zx')
  dzxoutputdx = zeros(nxpoints,nx);
  n = min(nxpoints,nulls.count);
  for i = 1:n
    ii = nulls.k(kxpoints(i))+neighbors;
    dzxoutputdx(i,:) = nulls.dzdpsi(kxpoints(i),:)*dpsizrdx(ii,:);
  end
  dydx(index_in_y.zx,:) = dzxoutputdx;
end
if isfield(index_in_y,'psix')
  dpsixoutputdx = zeros(nxpoints,nx);
  n = min(nxpoints,nulls.count);
  for i = 1:n
    ii = nulls.k(kxpoints(i))+neighbors;
    dpsixoutputdx(i,:) = nulls.w(kxpoints(i),:)*dpsizrdx(ii,:);
  end
  dydx(index_in_y.psix,:) = dpsixoutputdx;
end
if isfield(index_in_y,'rcur')
  dydx(index_in_y.rcur,:) = drcurdx;
end
if isfield(index_in_y,'zcur')
  dydx(index_in_y.zcur,:) = dzcurdx;
end
if isfield(index_in_y,'cpasma')
  dydx(index_in_y.cpasma,:) = dcpasmadx;
end
if isfield(index_in_y,'aminor')
  dydx(index_in_y.aminor,:) = da0dx;
end
if isfield(index_in_y,'rbbbs')
  dydx(index_in_y.rbbbs,:) = drbdx;
end
if isfield(index_in_y,'zbbbs')
  dydx(index_in_y.zbbbs,:) = dzbdx;
end
if isfield(index_in_y,'rmaxis')
  dydx(index_in_y.rmaxis,:) = drmaxisdx;
  if plasma_is_tiny
    dydx(index_in_y.rmaxis,:) = 0;
  end
end
if isfield(index_in_y,'zmaxis')
  dydx(index_in_y.zmaxis,:) = dzmaxisdx;
end
if isfield(index_in_y,'psimag')
  dydx(index_in_y.psimag,:) = dpsimagdx;
end
if isfield(index_in_y,'rbdef')
  dydx(index_in_y.rbdef,:) = drbdefdx;
end
if isfield(index_in_y,'zbdef')
  dydx(index_in_y.zbdef,:) = dzbdefdx;
end
if isfield(index_in_y,'psibry')
  dydx(index_in_y.psibry,:) = dpsibrydx;
end
if isfield(index_in_y,'li')
  dydx(index_in_y.li,:) = dlidx;
end
if isfield(index_in_y,'betap')
  dydx(index_in_y.betap,:) = dbetapdx;
end
if isfield(index_in_y,'Wth')
  dydx(index_in_y.Wth,:) = dWthdx;
end
if isfield(index_in_y,'psipla')
  dydx(index_in_y.psipla,:) = dpsipladx;
end
if isfield(index_in_y,'psiplaapp')
  dydx(index_in_y.psiplaapp,:) = dpsiplaappdx;
end
if isfield(index_in_y,'vplas') % This one multiplies by xdot, not dx
  dydx(index_in_y.vplas,:) = dpsipladx;
  % KLUGE
  dydx(index_in_y.vplas,:) = 0;
  dydx(index_in_y.vplas,1:nc) = wa*mpc(iia,:);
end
if isfield(index_in_y,'bp2flx')
  dydx(index_in_y.bp2flx,:) = dbp2flxdx;
end
if isfield(index_in_y,'Lpla') % This is an output made by gs_response
  lae.y(index_in_y.Lpla) = dpsipladx*dxdxc(:,nc+nv+1);
end
if isfield(index_in_y,'Ltot')
  dydx(index_in_y.Ltot,:) = dLtotdx;
end
if isfield(index_in_y,'Atot')
  dydx(index_in_y.Atot,:) = dAtotdx;
end
if isfield(index_in_y,'Vtot')
  dydx(index_in_y.Vtot,:) = dVtotdx;
  if plasma_is_tiny
    dydx(index_in_y.Vtot,:) = 0;
  end
end
if isfield(index_in_y,'rhot')
  dydx(index_in_y.rhot,:) = drhotdx;
end
if isfield(index_in_y,'jtav')
  dydx(index_in_y.jtav,:) = djtavdx;
end
if isfield(index_in_y,'Vres')
  dydx(index_in_y.Vres,:) = dVresdx;
end
for i = 1:length(user_signal)
  if isfield(index_in_y,user_signal(i).name)
    n = length(user_signal(i).y0);
    dydx(getfield(index_in_y,user_signal(i).name),:) = ...
      [user_signal(i).dydic user_signal(i).dydiv zeros(n,nx-nc-nv)] + ...
      user_signal(i).dydpcurrt*dpcurrtdx;
  end
end
