function [rb,zb,rx,zx,ra,za,r0,z0,ilimited,psimag,psibry]=trace_boundary(psizr,rg,zg,limdata)
%
%  USAGE:   [rb,zb,rx,zx,ra,za,r0,z0,ilimited,psimag,psibry]=trace_boundary(psizr,rg,zg,limdata)
%
%  PURPOSE: Find boundary
%
%  INPUTS:  psizr
%
%  OUTPUTS: Coordinates of boundary: rb, zb
%           Coordinates of x or touch point: r0,z0
%           Coordinates of axis: ra, za
%
%  RESTRICTIONS:
%
%  METHOD: 
%
	
%  VERSION @(#)trace_boundary.m	1.1 01/12/11
%
%  WRITTEN BY:  Anders Welander  ON	3/31/10
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Prelims
  mu0 = .4e-6*pi;
  twopi = 2*pi;
  nr = length(rg); nz = length(zg);
  dr = rg(2)-rg(1); dz = zg(2)-zg(1);
  rb=0; zb=0; rx=0; zx=0; ra=0; za=0; r0=0; z0=0; ilimited=0; % Always return values to avoid error
  
  if exist('limdata','var')
    if size(limdata,2)<size(limdata,1), limdata=limdata'; end
    nlim = size(limdata,2);
    rlm = interp1(1:nlim,limdata(2,:),1:.01:nlim);
    zlm = interp1(1:nlim,limdata(1,:),1:.01:nlim);
    psilim = interp2(rg,zg,psizr,rlm,zlm,'spline');
  end
  
  % Determine whether the axis is at max or min flux
  [psimax, imax] = max(psizr); zmax = zg(imax); [psimax, imax] = max(psimax); rmax = rg(imax); zmax=zmax(imax);
  [psimin, imin] = min(psizr); zmin = zg(imin); [psimin, imin] = min(psimin); rmin = rg(imin); zmin=zmin(imin);
  % Pick the extremum closest to the middle of the grid.
  if (rmax-mean(rg))^2+(zmax-mean(zg))^2 < (rmin-mean(rg))^2+(zmin-mean(zg))^2
    ra = rmax; za = zmax; % First stab at ra, za
  else
    ra = rmin; za = zmin;
    psizr = -psizr; % This way the flux always has max at axis
  end
 
  % Find weights in grid points to calculate value in a point using cubic Hermite spline (ref. wikipedia)
  mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;
  neighbors = reshape([-1-nz -1 -1+nz -1+2*nz;-nz 0 nz 2*nz;1-nz 1 1+nz 1+2*nz;2-nz 2 2+nz 2+2*nz],1,16);
  
  % Magnetic axis
  for j=1:50
    iz = max(find(zg<za)); ir = max(find(rg<ra)); ia = iz+nz*(ir-1);
    if length(ia)==0
      disp('trace_boundary: Could not find axis.')
      return
    end
    iia = ia+neighbors; % iia indexes 16 points around magnetic axis
    if min(iia)<1 | max(iia)>nr*nz
      disp('trace_boundary: Could not find axis.')
      return
    end
    pp = psizr(iia);
    tr = (ra-rg(ir))/dr; tz = (za-zg(iz))/dz;
    wa = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    war = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    waz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    warr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
    wazz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    warz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
    ashift_Ba = -inv([pp*warr' pp*warz';pp*warz' pp*wazz']);
    c = ashift_Ba*[pp*war'; pp*waz']; ra = ra+c(1)/5; za = za+c(2)/5;
  end
  psimag = pp*wa';

  % Find x-point
  na = 2*(nr-1); psix = -1e99; rho = linspace(0,sqrt((zg(end)-zg(1))^2+(rg(end)-rg(1))^2)/2,100);
  for ia = 1:na
    r = ra+rho*cos(ia*2*pi/na); z = za+rho*sin(ia*2*pi/na);
    p = interp2(rg,zg,psizr,r,z,'spline');
    j = 5;
    while p(j+1)<p(j) & j<99, j=j+1; end
    if p(j) > psix
      psix = p(j);
      rx = r(j); zx=z(j);
    end
  end

  % Now we are close
  for j=1:30
    if rx<rg(3) | rx>rg(end-3) | zx<zg(3) | zx>zg(end-3)
      if exist('limdata','var')
        ilimited = 1; j=99;
      else
        disp('trace_boundary: Could not find x-point.')
        return
      end
    end
    if ~ilimited
      iz = max(find(zg<zx)); ir = max(find(rg<rx)); ix = iz+nz*ir;
      iix = ix+neighbors;
      pp = psizr(iix);
      tr = (rx-rg(ir))/dr; tz = (zx-zg(iz))/dz;
      wx = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wxr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      wxz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wxrr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
      wxzz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wxrz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      xshift_Bx = -inv([pp*wxrr' pp*wxrz';pp*wxrz' pp*wxzz']);
      c = xshift_Bx*[pp*wxr'; pp*wxz']; rx = rx+c(1)/5; zx = zx+c(2)/5;
    end
  end
  if ilimited
    [psibry, ilim] = max(psilim);
    r0 = rlm(ilim); z0 = zlm(ilim);
  else
    psibry = pp*wx';
    w0 = wx; r0 = rx; z0 = zx;
  end

  th0 = angle(r0-ra+sqrt(-1)*(z0-za)); dth = 2*pi/na;
  % Trace boundary
  rb([1 na+1]) = r0; zb([1 na+1]) = z0;
  for ia = 2:na
    th = th0+dth*(ia-1);
    r = ra+rho*cos(th); z = za+rho*sin(th);
    p = interp2(rg,zg,psizr,r,z,'spline');
    j = 5;
    while p(j+1)<p(j) & j<99 & p(j)>psibry, j=j+1; end
    if p(j) > psibry
      if j<99 % In this case we are close to the other x-point
        rx2 = r(j); zx2 = z(j);
	rx2a = r(min(100,j+3)); zx2a = z(min(100,j+3));
	rx2b = r(max(1,j-3)); zx2b = z(max(1,j-3));
	for k=1:9
	  iz = max(find(zg<zx2)); ir = max(find(rg<rx2)); ix = iz+nz*ir;
	  if length(ix)==0 | iz<3 | iz>nz-3 | ir<3 | ir>nr-3 % Found no x-point so set point like normal
	    rx2 = spline(p(3:j),r(3:j),psibry);
	    zx2 = spline(p(3:j),z(3:j),psibry);
	    k = 9;
	  else
	    iix = ix+neighbors;
	    pp = psizr(iix);
	    tr = (rx2-rg(ir))/dr; tz = (zx2-zg(iz))/dz;
	    wx = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	    wxr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	    wxz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	    wxrr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
	    wxzz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	    wxrz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	    xshift_Bx = -inv([pp*wxrr' pp*wxrz';pp*wxrz' pp*wxzz']);
	    c = xshift_Bx*[pp*wxr'; pp*wxz']; rx2 = rx2+c(1); zx2 = zx2+c(2);
	  end
	end
	c = sort([rx2a rx2 rx2b]); rb(ia) = c(2);
	c = sort([zx2a zx2 zx2b]); zb(ia) = c(2);
      else
        disp('trace_boundary: error in tracing boundary. Could not find flux value below psibry.')
      end
    else
      rb(ia) = spline(p(3:j),r(3:j),psibry); zb(ia) = spline(p(3:j),z(3:j),psibry);
    end
  end
  if ~ilimited & exist('limdata','var')
    k = find(zlm<max(zb) & zlm>min(zb));
    if max(psilim(k))>psibry
      ilimited = 1;
      [psibry, ilim] = max(psilim(k));
      r0 = rlm(k(ilim)); z0 = zlm(k(ilim));
      th0 = angle(r0-ra+sqrt(-1)*(z0-za)); dth = 2*pi/na;
      % Trace boundary again because it was limited
      rb([1 na+1]) = r0; zb([1 na+1]) = z0;
      for ia = 2:na
	th = th0+dth*(ia-1);
	r = ra+rho*cos(th); z = za+rho*sin(th);
	p = interp2(rg,zg,psizr,r,z,'spline');
	j = 5;
	while p(j+1)<p(j) & j<99 & p(j)>psibry, j=j+1; end
	if p(j) > psibry
	  if j<99 % In this case we are close to the other x-point
            rx2 = r(j); zx2 = z(j);
	    rx2a = r(min(100,j+3)); zx2a = z(min(100,j+3));
	    rx2b = r(max(1,j-3)); zx2b = z(max(1,j-3));
	    for k=1:9
	      iz = max(find(zg<zx2)); ir = max(find(rg<rx2)); ix = iz+nz*ir;
	      iix = ix+neighbors;
	      pp = psizr(iix);
	      tr = (rx2-rg(ir))/dr; tz = (zx2-zg(iz))/dz;
	      wx = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	      wxr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	      wxz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	      wxrr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
	      wxzz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	      wxrz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	      xshift_Bx = -inv([pp*wxrr' pp*wxrz';pp*wxrz' pp*wxzz']);
	      c = xshift_Bx*[pp*wxr'; pp*wxz']; rx2 = rx2+c(1); zx2 = zx2+c(2);
	    end
	    c = sort([rx2a rx2 rx2b]); rb(ia) = c(2);
	    c = sort([zx2a zx2 zx2b]); zb(ia) = c(2);
	  else
            disp('trace_boundary: error in tracing boundary. Could not find flux value below psibry.')
	  end
	else
	  rb(ia) = spline(p(3:j),r(3:j),psibry); zb(ia) = spline(p(3:j),z(3:j),psibry);
	end
      end
    end
  end
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
