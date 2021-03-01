function [xresp,x,bresp,b,xtarget] = calc_x(xtarget,tok_data_struct,resp,eq,idoplot)
%
%  USAGE: [xresp,x,bresp,b,xtarget] = calc_x(xtarget,tok_data_struct,resp,eq,idoplot)
%
%  PURPOSE: Find x-points and responses thereof *in vicinity of* specified target points
%           and also return Br, Bz and their responses *at* specified target points.
%
%  INPUTS: xtarget: specification on the form [r z] where r and z are column vectors
%                   of approximate or target x-point locations.
%                   If xtarget is empty, the default is bottom and top of boundary.
%  tok_data_struct: toksys description of tokamak
%             resp: (optional) the output from a response model (gspert or rzrig)
%               eq: equilibrium structure
%          idoplot: flag to plot limiter, boundary, xtarget points and x-points
%
%  OUTPUTS: xresp: response structure with fields:
%                  drxdis: rx response to conductors [m/A]
%                  dzxdis: zx response to conductors [m/A]
%                  drxdip,*dli,*dbetap, etc: exogenous responses
%               x: the x-points on form [rx zx] found near xtarget points
%           bresp: response structure with fields:
%                  dbrdis: br response at points xtarget to conductors [T/A]
%                  dbzdis: bz response at points xtarget to conductors [T/A]
%                  dbrdip,*dli,*dbetap, etc: exogenous responses
%               b: The fields at xtarget points on the form [br bz]
%         xtarget: Target location for x-points (same as input if input is not empty)

%
%  RESTRICTIONS: x-points must be within the grid

%
%  METHOD: 
%	
%  VERSION @(#)calc_x.m	1.2 01/31/12
%
%  WRITTEN BY:  Anders Welander  ON	1/31/12
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if nargin<4
    help calc_x
    disp('*********************************************************')
    disp('Supply at least 4 arguments! (xtarget, resp can be empty)')
    disp('*********************************************************')
    return
  end
  
  if ~exist('idoplot','var')
    idoplot = 0;
  end
 
  % Unpack tok_data_struct and eq
  struct_to_ws(tok_data_struct);
  dr = mean(diff(rg)); dz = mean(diff(zg));
  if size(limdata,2)==2, limdata=limdata'; end
  rlim = limdata(2,:); zlim = limdata(1,:);
  if min(rlim)<0, rlim = limdata(1,:); zlim = limdata(2,:); end
  if rlim(1)~=rlim(end) | zlim(1)~=zlim(end)
    rlim(end+1) = rlim(1); zlim(end+1) = zlim(1);
  end
  nlim = length(rlim);
  
  mismatch = exp(max(abs(log([eq.dr/dr eq.dz/dz eq.rg(1)/rg(1) eq.zg(1)/zg(1) eq.nw/nr eq.nh/nz]))))*100-100;
  if mismatch > 1e-4
    if isempty(resp)
      disp(['Warning in calc_x: The grids in tok_data_struct and eq differ by ' num2str(mismatch) '%!'])
    else
      disp(['Error in calc_x: The grids in tok_data_struct and eq differ by ' num2str(mismatch) '%!'])
      disp('xresp can not be calculated.')
      pause
      return
    end
  end
  struct_to_ws(eq);
  
  if isempty(xtarget) % Use default xtarget, i.e. the top and bottom points of the boundary
    [zbot, k] = min(zbbbs);
    xtarget(1,:) = [rbbbs(k) zbbbs(k)];
    [ztop, k] = max(zbbbs);
    xtarget(2,:) = [rbbbs(k) zbbbs(k)];
  end
 
  nx = size(xtarget,1);
  xresp = []; x = [];
  bresp = []; b = [];

  % Find weights in grid points to calculate value at a point using cubic Hermite spline (ref. wikipedia)
  mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;
  neighbors = reshape([-1-nz -1 -1+nz -1+2*nz;-nz 0 nz 2*nz;1-nz 1 1+nz 1+2*nz;2-nz 2 2+nz 2+2*nz],1,16);
  
  % Find x-points in vicinity of xtarget points
  for ix = 1:nx
    c_previous = zg(end)-zg(1); cx = c_previous*.99; j = 9;
    rx = xtarget(ix,1); zx = xtarget(ix,2);
    % Try zooming in on x-point with Newton Rhapson. Also calculate Br, Bz at xtarget points
    while j>0 & norm(cx)<norm(c_previous)
      % Find indices and weights for grid points around the x-point
      kr0 = min(nr-3,max(1,floor((rx-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
      kz1 = min(nz-2,max(2,ceil((zx-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
      k = kr0*nz+kz1;
      iix = k+neighbors; % iix indexes 16 points around x point
      pp = psizr(iix); j = j-1; c_previous = cx;
      tr = (rx-rgg(k))/dr; tz = (zx-zgg(k))/dz;
      wx = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wxr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      wxz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wxrr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
      wxzz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wxrz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      % Calculate b at xtarget points
      if size(b,1) < ix % This is true on the first pass
        b(ix,1) = -pp*wxz'/rx/2/pi;
        b(ix,2) = +pp*wxr'/rx/2/pi;
      end
      xshift_Bx = -inv([pp*wxrr' pp*wxrz';pp*wxrz' pp*wxzz']);
      cx = xshift_Bx*[pp*wxr'; pp*wxz'];
      rx = rx+cx(1); zx = zx+cx(2);
    end
    if norm(cx)<dr % The test that we converged on an x-point
      x(ix,1) = rx;  x(ix,2) = zx;
    else
      x(ix,1) = NaN; x(ix,2) = NaN;
    end
  end
  
  if idoplot
    plot(xtarget(:,1),xtarget(:,2),'go','markerSize',8)
    hold on
    plot(x(:,1),x(:,2),'rx','markerSize',12,'linew',2)
    plot(rlim,zlim,'k','linew',2)
    plot(rbbbs(1:nbbbs),zbbbs(1:nbbbs),'b')
    axis('image')
    legend('xtarget','x-points');
    hold off
  end
  
  if isempty(resp)
    return
  end
  
  % Parse the resp structure
  if isfield(resp,'dcphidis') % true for gspert
  elseif isfield(resp,'drdis') & isfield(resp,'dcdr') % true for rzrig
    resp.dcphidis = resp.dcdr(:)*resp.drdis+resp.dcdz(:)*resp.dzdis;
    resp.units.dcphidis = 'A/A';
    resp.dcphidip = resp.dcdr(:)*resp.drdip+resp.dcdz(:)*resp.dzdip;
    resp.units.dcphidip = 'A/A';
    resp.dcphidli = resp.dcdr(:)*resp.drdli;
    resp.units.dcphidli = 'A/1';
    resp.dcphidbetap = resp.dcdr(:)*resp.drdbetap+resp.dcdz(:)*resp.dzdbetap;
    resp.units.dcphidbetap = 'A/1';
  else
    wait('error calc_gap: can not interpret structure resp')
    return
  end
  % Since exogenous variables may vary, use field names to extract all that apply
  f = fieldnames(resp); % f holds names of fields in resp
  ind_resp_in_f = []; % will become indices in f containing dcphid*
  for j = 1:length(f)
    if strfind(char(f(j)),'dcphid') == 1
      ind_resp_in_f(end+1) = j;
    end
  end
  
  % Now calculate x-point responses using dcphid* for plasma response and tok_data_struct for mutuals
  for j = 1:nx
    rx = x(j,1); zx = x(j,2);
    % Get mutuals from grid to the 16 points that surround this point
    kr0 = min(nr-3,max(1,floor((rx-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((zx-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    k = kr0*nz+kz1;
    ii = k+neighbors; % ii indexes 16 points around xtarget location
    m = zeros(16,nz*nr); % will hold mutuals to the 16 points around x(j,:)
    for jj = 1:16
      z_idx = mod(ii(jj)-1,nz)+1; % z index of output element ii(jj)
      r_idx = ceil(ii(jj)/nz);    % r index of output element ii(jj)
      for jr=1:nr
        k1 = (r_idx-1)*nz;
	k2 = (jr-1)*nz;
	for jz=1:nz
	  m(jj,k2+jz) = mpp(k1+abs(z_idx-jz)+1,jr);
	end
      end
    end
    % We now have dpsizr(ii)/dis = m*resp.dcphidis+[mpc(ii,:) mpv(ii,:)], etc
    tr = (rx-rgg(k))/dr; tz = (zx-zgg(k))/dz;
    wx = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wxr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    wxz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wxrr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
    wxzz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wxrz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    pp = psizr(ii);
    xshift_Bx = -inv([pp*wxrr' pp*wxrz';pp*wxrz' pp*wxzz']);
    % We now have dpsi(rx,zx)/dis = w*(m*dpsizr(ii))/dis, etc
    for jj = 1:length(ind_resp_in_f)
      respobj = char(f(ind_resp_in_f(jj)));
      if strmatch(respobj,'dcphidis')
        dpsidis = [(m*resp.dcphidis+[mpc(ii,:) mpv(ii,:)])];
        cx = xshift_Bx*[wxr*dpsidis; wxz*dpsidis];
        xresp.drxdis(j,:) = cx(1,:);
        xresp.dzxdis(j,:) = cx(2,:);
      elseif eval(['prod(size(resp.' respobj '))']) == nr*nz % respobj may be on form (nz,nr) so collapse
        eval(['dpsidx = m*resp.' respobj '(:);'])
	cx = xshift_Bx*[wxr*dpsidx; wxz*dpsidx];
        eval(['xresp.drxd' respobj(7:end) '(j,:) = cx(1,:);'])
        eval(['xresp.dzxd' respobj(7:end) '(j,:) = cx(2,:);'])
      else
        eval(['dpsidx = m*resp.' respobj ';'])
	cx = xshift_Bx*[wxr*dpsidx; wxz*dpsidx];
        eval(['xresp.drxd' respobj(7:end) '(j,:) = cx(1,:);'])
        eval(['xresp.dzxd' respobj(7:end) '(j,:) = cx(2,:);'])
      end
      if isnan(rx)
        eval(['xresp.drxd' respobj(7:end) '(j,:) = NaN;'])
        eval(['xresp.dzxd' respobj(7:end) '(j,:) = NaN;'])
      end
    end
  end
  for jj = 1:length(ind_resp_in_f)
    respobj = char(f(ind_resp_in_f(jj)));
    eval(['xresp.descriptions.drxd' respobj(7:end) ' = ''rx response to ' respobj(7:end) ''';'])
    eval(['xresp.descriptions.dzxd' respobj(7:end) ' = ''zx response to ' respobj(7:end) ''';'])
    s = eval(['resp.units.dcphid' respobj(7:end)]);
    eval(['xresp.units.drxd' respobj(7:end) ' = [''m/A * '' s]; '])
    eval(['xresp.units.dzxd' respobj(7:end) ' = [''m/A * '' s]; '])
  end
  
  % Now calculate Br, Bz responses at xtarget using dcphid* for plasma response and tok_data_struct for mutuals
  for j = 1:nx
    rx = xtarget(j,1); zx = xtarget(j,2);
    % Get mutuals from grid to the 16 points that surround this point
    kr0 = min(nr-3,max(1,floor((rx-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((zx-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    k = kr0*nz+kz1;
    ii = k+neighbors; % ii indexes 16 points around xtarget location
    m = zeros(16,nz*nr); % will hold mutuals to the 16 points around xtarget(j,:)
    for jj = 1:16
      z_idx = mod(ii(jj)-1,nz)+1; % z index of output element ii(jj)
      r_idx = ceil(ii(jj)/nz);    % r index of output element ii(jj)
      for jr=1:nr
        k1 = (r_idx-1)*nz;
	k2 = (jr-1)*nz;
	for jz=1:nz
	  m(jj,k2+jz) = mpp(k1+abs(z_idx-jz)+1,jr);
	end
      end
    end
    % We now have dpsizr(ii)/dis = m*resp.dcphidis+[mpc(ii,:) mpv(ii,:)], etc
    tr = (rx-rgg(k))/dr; tz = (zx-zgg(k))/dz;
    w = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    wr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    wz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    % We now have dpsi(rx,zx)/dis = w*(m*dpsizr(ii))/dis, etc
    for jj = 1:length(ind_resp_in_f)
      respobj = char(f(ind_resp_in_f(jj)));
      if strmatch(respobj,'dcphidis')
        dpsidis = [(m*resp.dcphidis+[mpc(ii,:) mpv(ii,:)])];
        bresp.dbrdis(j,:) = -wxz*dpsidis/rx/2/pi;
        bresp.dbzdis(j,:) = +wxr*dpsidis/rx/2/pi;
      elseif eval(['prod(size(resp.' respobj '))']) == nr*nz % respobj may be on form (nz,nr) so collapse
        eval(['dpsidx = m*resp.' respobj '(:);'])
        eval(['bresp.dbrd' respobj(7:end) '(j,:) = -wxz*dpsidx/rx/2/pi;'])
        eval(['bresp.dbzd' respobj(7:end) '(j,:) = +wxr*dpsidx/rx/2/pi;'])
      else
        eval(['dpsidx = m*resp.' respobj ';'])
        eval(['bresp.dbrd' respobj(7:end) '(j,:) = -wxz*dpsidx/rx/2/pi;'])
        eval(['bresp.dbzd' respobj(7:end) '(j,:) = +wxr*dpsidx/rx/2/pi;'])
      end
      if isnan(rx) % Not likely to happen for xtarget
        eval(['bresp.dbrd' respobj(7:end) '(j,:) = NaN;'])
        eval(['bresp.dbzd' respobj(7:end) '(j,:) = NaN;'])
      end
    end
  end
  for jj = 1:length(ind_resp_in_f)
    respobj = char(f(ind_resp_in_f(jj)));
    eval(['bresp.descriptions.dbrd' respobj(7:end) ' = ''br response to ' respobj(7:end) ''';'])
    eval(['bresp.descriptions.dbzd' respobj(7:end) ' = ''bz response to ' respobj(7:end) ''';'])
    s = eval(['resp.units.dcphid' respobj(7:end)]);
    eval(['bresp.units.dbrd' respobj(7:end) ' = [''T/A * '' s]; '])
    eval(['bresp.units.dbzd' respobj(7:end) ' = [''T/A * '' s]; '])
  end
  
  
