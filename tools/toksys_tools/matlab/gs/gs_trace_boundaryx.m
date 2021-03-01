%  USAGE:   gs_trace_boundaryx
%
%  PURPOSE: Trace the plasma boundary, this x-version also finds the x-point
%
%  INPUTS: psizr, the flux at nz vertical positions * nr radial positions
%          psibry, the flux at the boundary (see gs_find_bdef)
%          ilimgg, flags for grid, -1 = outside vessel, >-1 = inside
%          nbbbs_max, maximum number of boundary points
%          zbot, ztop, min and max z position of plasma
%          zbotg = (floor((zbot-zg(1))/dz)+0.5)*dz+zg(1);
%          ztopg = (floor((ztop-zg(1))/dz)+0.5)*dz+zg(1);
%          rbdef, zbdef, position that defines the boundary (touch or x point)
%          ix1, flag which is true if lower x-point defines the boundary
%              if ix1 is true then these derivates of the flux are needed:
%              psix1_rr, psix1_rz, psix1_zz (see gs_find_bdef)
%          ix2, flag which is true if upper x-point defines the boundary
%              if ix2 is true then these derivates of the flux are needed:
%              psix2_rr, psix2_rz, psix2_zz (see gs_find_bdef)
%          rmaxis, zmaxis, position of axis (used to order points w.r.t theta)
%          rgg, zgg, dr, dz, nr, nz (grid variables)
%          R13 = (1+1i*sqrt(3))^2/4; % for solving cubic equations
%          For cubic interpolation on the grid:
%            mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%            neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
%
%  OUTPUTS: rbbbs, zbbbs, nbbbs, R, Z of boundary and number of points
%           rhobbbs, thbbbs, distance to magnetic axis and poloidal angle
%	
%  METHOD:  The grid is searched for points that are on either side of the boundary
%           The exact locations are then found by solving 1-dim cubic equations
%           The boundary defining point is added and also extra points if needed
%           to make all poloidal angles between adjacent points < dthbbbs_threshold 
	
%  NOTES:  The points that are added in the event dtheta > dthbbbs_threshold
%          are dthbbbs_max from previous points rather than evenly distributed
%          in the gap between regular points. The reason is to ensure smoothness
%          in the transition from needing an extra point to not needing it, 
%          when the plasma changes.
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	3/12/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


trace_status = -1;

xh = nan(ngg,3); % Solutions between horizontal neighbors
xv = nan(ngg,3); % Solutions between vertical neighbors
fh = logical(ones(ngg,3)); % Flags for okay to visit
fv = logical(ones(ngg,3)); % Flags for okay to visit

gbbbs(:) = 0;
j = floor((rmaxis-rg(1))/dr)+1;
i = floor((zmaxis-zg(1))/dz)+1;
k = i+(j-1)*nz;
xh(k,:) = gs_solve_hermite_cubic(psizr(i,j-1:j+2), psibry);
while sum(imag(xh(k,:)) == 0 & xh(k,:) >=0 & xh(k,:) < 1) == 0 & j > 2
  j = j-1;
  k = i+(j-1)*nz;
  xh(k,:) = gs_solve_hermite_cubic(psizr(i,j-1:j+2), psibry);
end
if sum(imag(xh(k,:)) == 0 & xh(k,:) >=0 & xh(k,:) < 1) == 0
  trace_status = 1;
  return
else
  n = 1;
  gbbbs(1) = nz;
  ibbbs(1) = k;
  for p = 1:3
    if sum(fh(k,:)) == 3 & imag(xh(k,p)) ==0 & xh(k,p) >= 0 & xh(k,p) < 1
      fh(k,p) = 0;
      xx(1) = xh(k,p);
    end
  end
end

gradmin = inf; % Will become lowest flux gradient found
ix = -1; % Index in bbbs to boundary-defining x-point

while trace_status < 0
  k = i+(j-1)*nz;
  if isnan(xh(k,1))
    xh(k,:) = gs_solve_hermite_cubic(psizr(i,j-1:j+2), psibry);
  end
  if isnan(xv(k,1))
    xv(k,:) = gs_solve_hermite_cubic(psizr(i-1:i+2,j), psibry);
  end
  if isnan(xh(k+1,1))
    xh(k+1,:) = gs_solve_hermite_cubic(psizr(i+1,j-1:j+2), psibry);
  end
  if isnan(xv(k+nz,1))
    xv(k+nz,:) = gs_solve_hermite_cubic(psizr(i-1:i+2,j+1), psibry);
  end
  f12 = [fh(k,:) fv(k,:) fh(k+1,:) fv(k+nz,:)];
  r12 = ~[imag(xh(k,:)) imag(xv(k,:)) imag(xh(k+1,:)) imag(xv(k+nz,:))];
  k12 = [xh(k,:) >= 0 & xh(k,:) < 1, ...
         xv(k,:) >= 0 & xv(k,:) < 1, ...
	 xh(k+1,:) >= 0 & xh(k+1,:) < 1, ...
         xv(k+nz,:) >= 0 & xv(k+nz,:) < 1];

  if sum(r12 & k12 & f12) == 0
    trace_status = 0;
  elseif sum(r12 & k12 & f12) == 1 % There is only one way to go
    n = n+1;
    for p = 1:3
      if r12(p) & k12(p) & f12(p)
	fh(k,p) = 0;
        xx(n) = xh(k,p);
        gbbbs(n) = nz;
        ibbbs(n) = k;
	i = i-1;
      end
      if r12(p+3) & k12(p+3) & f12(p+3)
	fv(k,p) = 0;
        xx(n) = xv(k,p);
        gbbbs(n) = 1;
        ibbbs(n) = k;
	j = j-1;
      end
      if r12(p+6) & k12(p+6) & f12(p+6)
	i = i+1;
	fh(k+1,p) = 0;
        xx(n) = xh(k+1,p);
        gbbbs(n) = nz;
        ibbbs(n) = k+1;
      end
      if r12(p+9) & k12(p+9) & f12(p+9)
	j = j+1;
	fv(k+nz,p) = 0;
        xx(n) = xv(k+nz,p);
        gbbbs(n) = 1;
        ibbbs(n) = k+nz;
      end
    end
  else % More than one way to go
    if gbbbs(n) == 1
      ir = (ibbbs(n)-k)/nz;
      iz = xx(n);
    else
      ir = xx(n);
      iz = ibbbs(n)-k;
    end
    [ir, iz, m, vr, vz] = gs_trace_step_v(ir, iz, 0, 1, 0, 1, psizr(i-1:i+2,j-1:j+2), psibry);
    for p = 1:length(vr)
      tr = vr(p);
      tz = vz(p);
      wr =  reshape(([1 tz   tz^2   tz^3]*mx)'*[0 1  2*tr   3*tr^2]*mx,1,16);
      wz =  reshape(([0 1  2*tz   3*tz^2]*mx)'*[1 tr   tr^2   tr^3]*mx,1,16);
      dum = (wr*psizr(k+neighbors'))^2+(wz*psizr(k+neighbors'))^2;
      if dum < gradmin
        gradmin = dum;
	iix = k+neighbors';
	wxr = wr;
	wxz = wz;
	rx = rg(j)+tr*dr;
	zx = zg(i)+tz*dz;
	ix = n+1;
      end
    end
    n = n+1; % Update number of contour points
    if m == -4 % gs_trace_step stepped inward
      j = j-1;
      [~, p] = min(abs(xv(k,:)-iz));
      if fv(k,p)
	fv(k,p) = 0;
	ibbbs(n) = k;
	gbbbs(n) = 1;
	xx(n) = iz;
      else
	trace_status = 0;
      end
    elseif m == -1 % gs_trace_step stepped down
      i = i-1;
      [~, p] = min(abs(xh(k,:)-ir));
      if fh(k,p)
	fh(k,p) = 0;
	ibbbs(n) = k;
	gbbbs(n) = nz;
	xx(n) = ir;
      else
	trace_status = 0;
      end
    elseif m == +1 % gs_trace_step stepped up
      i = i+1;
      [~, p] = min(abs(xh(k+1,:)-ir));
      if fh(k+1,p)
	fh(k+1,p) = 0;
	ibbbs(n) = k+1;
	gbbbs(n) = nz;
	xx(n) = ir;
      else
	trace_status = 0;
      end
    elseif m == +4 % gs_trace_step stepped outward
      j = j+1;
      [~, p] = min(abs(xv(k+nz,:)-iz));
      if fv(k+nz,p)
	fv(k+nz,p) = 0;
	ibbbs(n) = k+nz;
	gbbbs(n) = 1;
	xx(n) = iz;
      else
	trace_status = 0;
      end
    elseif m == 0
      trace_status = 2;
    end
  end
  if trace_status == -1 & i < 2 | i > nz-2 | j < 2 | j > nr-2
    trace_status = 1;
  end
  if i+(j-1)*nz == ibbbs(1)
    trace_status = 0;
  end
end

if trace_status == 0
  
  jr = gbbbs == nz;
  jz = gbbbs == 1;
  rbbbs(jr) = rgg(ibbbs(jr))+xx(jr)*dr;
  rbbbs(jz) = rgg(ibbbs(jz));
  zbbbs(jr) = zgg(ibbbs(jr));
  zbbbs(jz) = zgg(ibbbs(jz))+xx(jz)*dz;
  thbbbs(1:n)  = angle(rbbbs(1:n)-rmaxis + 1i*zbbbs(1:n)-zmaxis*1i);

  % Order the points according to theta from -pi to +pi
  [~, i] = min(thbbbs(1:n));
  if sum(diff(thbbbs(1:n)) > 0) > sum(diff(thbbbs(1:n)) < 0)
    thbbbs(1:n)  = thbbbs([i:n, 1:i-1]);
    gbbbs(1:n) = gbbbs([i:n, 1:i-1]);
    ibbbs(1:n) = ibbbs([i:n, 1:i-1]);
    rbbbs(1:n) = rbbbs([i:n, 1:i-1]);
    zbbbs(1:n) = zbbbs([i:n, 1:i-1]);
    if ix >= 1 & ix <= n
      ix = find([i:n, 1:i-1]==ix);
    end
  else
    thbbbs(1:n)  = thbbbs([i:-1:1, n:-1:i+1]);
    gbbbs(1:n) = gbbbs([i:-1:1, n:-1:i+1]);
    ibbbs(1:n) = ibbbs([i:-1:1, n:-1:i+1]);
    rbbbs(1:n) = rbbbs([i:-1:1, n:-1:i+1]);
    zbbbs(1:n) = zbbbs([i:-1:1, n:-1:i+1]);
    if ix >= 1 & ix <= n
      ix = find([i:-1:1, n:-1:i+1]==ix);
    end
  end
  jr = gbbbs == nz;
  jz = gbbbs == 1;

  id4(jr,:) = [ibbbs(jr)-nz ibbbs(jr) ibbbs(jr)+nz ibbbs(jr)+2*nz];
  id4(jz,:) = [ibbbs(jz)- 1 ibbbs(jz) ibbbs(jz)+ 1 ibbbs(jz)+2   ];
  id8(jr,:) = [ibbbs(jr)-1-nz ibbbs(jr)-1 ibbbs(jr)-1+nz ibbbs(jr)-1+2*nz ...
               ibbbs(jr)+1-nz ibbbs(jr)+1 ibbbs(jr)+1+nz ibbbs(jr)+1+2*nz];
  id8(jz,:) = [ibbbs(jz)-nz-1 ibbbs(jz)-nz ibbbs(jz)-nz+1 ibbbs(jz)-nz+2 ...
               ibbbs(jz)+nz-1 ibbbs(jz)+nz ibbbs(jz)+nz+1 ibbbs(jz)+nz+2];

  wbbbs(1:n,:) = [ones(n,1) xx(1:n) xx(1:n).^2 xx(1:n).^3]*mx;
  wds(1:n,:) = [zeros(n,1) ones(n,1) 2*xx(1:n) 3*xx(1:n).^2]*mx;

  dpsid4(1:n) = sum(wds(1:n,:)'.*psizr(id4(1:n,:))')';
  dpsid8(1:n) = sum([-wbbbs(1:n,:) wbbbs(1:n,:)]'.*psizr(id8(1:n,:))')';

  dpsibbbsdr(jr) = dpsid4(jr)/dr;
  dpsibbbsdr(jz) = dpsid8(jz)/(2*dr);
  dpsibbbsdz(jr) = dpsid8(jr)/(2*dz);
  dpsibbbsdz(jz) = dpsid4(jz)/dz;

  rhobbbs(1:n) = sqrt((rbbbs(1:n)-rmaxis).^2+(zbbbs(1:n)-zmaxis).^2);

  nbbbs = n+1;
  gbbbs(nbbbs) = gbbbs(1);
  ibbbs(nbbbs) = ibbbs(1);
  rbbbs(nbbbs) = rbbbs(1);
  zbbbs(nbbbs) = zbbbs(1);
  wbbbs(nbbbs,:) = wbbbs(1,:);
  thbbbs(nbbbs) = thbbbs(1)+2*pi;
  rhobbbs(nbbbs) = rhobbbs(1);
  dpsibbbsdr(nbbbs) = dpsibbbsdr(1);
  dpsibbbsdz(nbbbs) = dpsibbbsdz(1);
  
  if ix == -1
    [~,i] = min(dpsibbbsdr(1:n).^2+dpsibbbsdz(1:n).^2);
    rx = rbbbs(i);
    zx = zbbbs(i);
    ir = (rx-rg(1))/dr+1;
    iz = (zx-zg(1))/dz+1;
    kr = floor(ir);
    kz = floor(iz);
    tr = ir-kr;
    tz = iz-kz;
    iix = kz+nz*(kr-1)+neighbors';
    if ~exist('drzdgp','var')
      drzdgp = zeros(2);
    end
    if ~exist('wr','var')
      wr = zeros(1,16);
      wz = zeros(1,16);
    end
    dzxdpsi = (drzdgp(2,1)*wr+drzdgp(2,2)*wz)*dz;
  end

  k = iix(6);
  tr = (rx-rgg(k))/dr;
  tz = (zx-zgg(k))/dz;
  pp = psizr(iix);
  m = 0;
  while m < 9
    m = m+1;
    wx =   reshape(([1 tz   tz^2   tz^3]*mx)'*[1 tr   tr^2   tr^3]*mx,1,16);
    wr =  reshape(([1 tz   tz^2   tz^3]*mx)'*[0 1  2*tr   3*tr^2]*mx,1,16);
    wz =  reshape(([0 1  2*tz   3*tz^2]*mx)'*[1 tr   tr^2   tr^3]*mx,1,16);
    wrr = reshape(([1 tz   tz^2   tz^3]*mx)'*[0 0  2      6*tr  ]*mx,1,16);
    wrz = reshape(([0 1  2*tz   3*tz^2]*mx)'*[0 1  2*tr   3*tr^2]*mx,1,16);
    wzz = reshape(([0 0    2    6*tz  ]*mx)'*[1 tr   tr^2   tr^3]*mx,1,16);
    psix_r = wr*pp;
    psix_z = wz*pp;
    psix_rr = wrr*pp;
    psix_rz = wrz*pp;
    psix_zz = wzz*pp;
    drzdgp = -inv([psix_rr psix_rz; psix_rz psix_zz]);
    cx = drzdgp*[psix_r; psix_z];
    if norm(cx) > 0.25
      cx = 0.25*cx/norm(cx);
    end
    tr = tr+cx(1);
    tz = tz+cx(2);
  end
  rx = rgg(k)+tr*dr;
  zx = zgg(k)+tz*dz;
  psix = wx*pp;
  drxdpsi = (drzdgp(1,1)*wr+drzdgp(1,2)*wz)*dr;
  dzxdpsi = (drzdgp(2,1)*wr+drzdgp(2,2)*wz)*dz;
  thx = angle(rx-rmaxis+1i*zx-1i*zmaxis);
  rhox = sqrt((rx-rmaxis)^2+(zx-zmaxis)^2);
  if ix == -1 % bbbs element with index i had weakest field
    if i > 1
      i1 = i-1;
    else
      i1 = nbbbs-1;
    end
    if (thbbbs(i1)-thx)*(thbbbs(i)-thx) < 0
      ix = i;
    else
      ix = i+1;
    end
  end

  % Squeeze in this x-point in bbbs
  nbbbs = nbbbs+1;
  gbbbs(ix+1:nbbbs) = gbbbs(ix:nbbbs-1);
  ibbbs(ix+1:nbbbs) = ibbbs(ix:nbbbs-1);
  rbbbs(ix+1:nbbbs) = rbbbs(ix:nbbbs-1);
  zbbbs(ix+1:nbbbs) = zbbbs(ix:nbbbs-1);
  wbbbs(ix+1:nbbbs,:) = wbbbs(ix:nbbbs-1,:);
  thbbbs(ix+1:nbbbs) = thbbbs(ix:nbbbs-1);
  rhobbbs(ix+1:nbbbs) = rhobbbs(ix:nbbbs-1);
  dpsibbbsdr(ix+1:nbbbs) = dpsibbbsdr(ix:nbbbs-1);
  dpsibbbsdz(ix+1:nbbbs) = dpsibbbsdz(ix:nbbbs-1);

  gbbbs(ix) = 0;
  ibbbs(ix) = iix(6);
  rbbbs(ix) = rx;
  zbbbs(ix) = zx;
  thbbbs(ix) = thx;
  rhobbbs(ix) = rhox;
  dpsibbbsdr(ix) = 0;
  dpsibbbsdz(ix) = 0;
  
  rbdef = rx;
  zbdef = zx;
  wb = wx;
  iib = iix;
  drbdefdpsi = drxdpsi;
  dzbdefdpsi = dzxdpsi;

  iplasma = logical(zeros(nz,nr));
  for i = 1:nz
    k = 0;
    for j = 1:nr
      iplasma(i,j) = k;
      if mod(sum(fh(i+nz*(j-1),:)==0),2)
	k = ~k;
      end
    end
  end

end


