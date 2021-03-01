%  USAGE:   gs_trace_boundary
%
%  PURPOSE: Trace plasma boundary
%
%  INPUTS: rbdef, zbdef, position that defines the boundary (touch or x point)
%          psizr, the flux at nz vertical positions * nr radial positions
%          psibry, the flux at the boundary (see gs_find_bdef)
%          rmaxis, zmaxis, position of axis (used to order points w.r.t theta)
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

pstart = 0;
while trace_status ~= 0 & pstart < 12

  pstart = pstart+1;
  trace_status = -1;

  gbbbs(:) = 0;
  rbbbs(1) = rbdef;
  zbbbs(1) = zbdef;
  n = 1;
  j = max(2,min(nr-2,floor((rbdef-rg(1))/dr)+1));
  i = max(2,min(nz-2,floor((zbdef-zg(1))/dz)+1));
  k = i+(j-1)*nz;
  ibbbs(1) = k;
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

  fh = logical(ones(ngg,3)); % Flags for okay to visit
  fv = logical(ones(ngg,3)); % Flags for okay to visit

  n = 2;
  
  if lim
    % Trace from touch point to edge of cell
    pstart = 12;
    ir = (rbdef-rg(j))/dr;
    iz = (zbdef-zg(i))/dz;
    ir0 = ir;
    iz0 = iz;
    wr0 = [1 ir ir^2 ir^3]*mx;
    wz0 = mx'*[1 iz iz^2 iz^3]';
    wr1 = [0 1 2*ir 3*ir^2]*mx/dr;
    wz1 = mx'*[0 1 2*iz 3*iz^2]'/dz;
    y44 = psizr(i-1:i+2,j-1:j+2);
    psi_r = reshape(wz0*wr1,1,16)*y44(:);
    psi_z = reshape(wz1*wr0,1,16)*y44(:);
    done_tracing = false;
    counter = 0;
    irect = 1;
    if abs(psi_r) > abs(psi_z)
      % Make horizontal cut in the cell, field line is kind of vertical
      while ~done_tracing
	counter = counter+1;
	if counter > 4
	  trace_status = 1;
	  done_tracing = true;
	end
	if irect == 1
	  [ir, iz, m] = gs_trace_step(ir, iz, 0, 1, 0, iz0, y44, psibry);
	  if m == 1
            irect = 2;
	  else
            done_tracing = true;
	  end
	else % irect = 2
	  [ir, iz, m] = gs_trace_step(ir, iz, 0, 1, iz0, 1, y44, psibry);
	  if m == -1
            irect = 1;
	  else
            done_tracing = true;
	  end
	end
      end
    else
      % Make vertical cut in the cell, field line is kind of horizontal
      while ~done_tracing
	counter = counter+1;
	if counter > 4
	  trace_status = 1;
	  done_tracing = true;
	end
	if irect == 1
	  [ir, iz, m] = gs_trace_step(ir, iz, 0, ir0, 0, 1, y44, psibry);
	  if m == 4
            irect = 2;
	  else
            done_tracing = true;
	  end
	else % irect = 2
	  [ir, iz, m] = gs_trace_step(ir, iz, ir0, 1, 0, 1, y44, psibry);
	  if m == -4
            irect = 1;
	  else
            done_tracing = true;
	  end
	end
      end
    end
    if m == -1
      [~, p] = min(abs(xh(k,:)-ir));
    elseif m == -4
      [~, p] = min(abs(xv(k,:)-iz));
      p = p+3;
    elseif m == +1
      [~, p] = min(abs(xh(k+1,:)-ir));
      p = p+6;
    else
      [~, p] = min(abs(xv(k+nz,:)-iz));
      p = p+9;
    end
    
  else % Testing x-point, assume all exits with psizr = psibry are accessible
    p = pstart;
    ir = (rbdef-rg(j))/dr;
    iz = (zbdef-zg(i))/dz;
    if p < 4
      jstart = xh(k,p);
      istart = 0;
    elseif p < 7
      jstart = 0;
      istart = xv(k,p-3);
    elseif p < 10
      jstart = xh(k+1,p-6);
      istart = 1;
    else
      jstart = 1;
      istart = xv(k+nz,p-9);
    end
    rzhit = [ir+[-1 1]*1e-4 iz+[-1 1]*1e-4];
    if rzhit(1) < 0
      rzhit(1:2) = [0 2]*1e-4;
    end
    if rzhit(3) < 0
      rzhit(3:4) = [0 2]*1e-4;
    end
    if rzhit(2) > 1
      rzhit(1:2) = 1-[2 0]*1e-4;
    end
    if rzhit(4) > 1
      rzhit(3:4) = 1-[2 0]*1e-4;
    end
    hit = 0;
    if imag(jstart) == 0 & jstart >= 0 & jstart <= 1 & ...
       imag(istart) == 0 & istart >= 0 & istart <= 1
      [ir, iz, m, hit] = gs_trace_step(...
        jstart, istart, 0, 1, 0, 1, psizr(i-1:i+2,j-1:j+2), psibry, ...
        xh(k,:), xv(k,:), xh(k+1,:), xv(k+nz,:), rzhit);
    end
    if ~hit
      trace_status = 1;
    end
  end

  
  if p < 4
    fh(k,p) = 0;
    xx(n) = xh(k,p);
    gbbbs(n) = nz;
    ibbbs(n) = k;
    i = i-1;
    r = rg(j) + xx(2)*dr;
    z = zg(i+1);
  elseif p < 7
    fv(k,p-3) = 0;
    xx(n) = xv(k,p-3);
    gbbbs(n) = 1;
    ibbbs(n) = k;
    j = j-1;
    r = rg(j+1);
    z = zg(i) + xx(2)*dz;
  elseif p < 10
    i = i+1;
    fh(k+1,p-6) = 0;
    xx(n) = xh(k+1,p-6);
    gbbbs(n) = nz;
    ibbbs(n) = k+1;
    r = rg(j) + xx(2)*dr;
    z = zg(i-1);
  else
    j = j+1;
    fv(k+nz,p-9) = 0;
    xx(n) = xv(k+nz,p-9);
    gbbbs(n) = 1;
    ibbbs(n) = k+nz;
    r = rg(j-1);
    z = zg(i) + xx(2)*dz;
  end

  if imag(xx(2)) == 0 & xx(2) >= 0 & xx(2) <= 1 & ...
     i > 1 & j > 1 & i < nz-1 & j < nr-1
  
    % Needed since the grid coincides with limiter element in DIII-D
    if ((r-rbdef)/dr)^2 + ((z-zbdef)/dz)^2 < 1e-8
      n = 1;
    end

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
	[ir, iz, m] = gs_trace_step(...
	  ir, iz, 0, 1, 0, 1, psizr(i-1:i+2,j-1:j+2), psibry, ...
	  xh(k,:), xv(k,:), xh(k+1,:), xv(k+nz,:));
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
      if trace_status == -1
        if i < 2 | i > nz-2 | j < 2 | j > nr-2
	  trace_status = 1;
	elseif ilimgg(i,j) == -1
	  trace_status = 1;
	end
      end
      if i+(j-1)*nz == ibbbs(1)
	trace_status = 0;
      end
    end % End of the while loop:    while trace_status < 0

    if trace_status == 0
    
      jr = gbbbs == nz;
      jz = gbbbs == 1;
      rbbbs(jr) = rgg(ibbbs(jr))+xx(jr)*dr;
      rbbbs(jz) = rgg(ibbbs(jz));
      zbbbs(jr) = zgg(ibbbs(jr));
      zbbbs(jz) = zgg(ibbbs(jz))+xx(jz)*dz;
      thbbbs(1:n)  = angle(rbbbs(1:n)-rmaxis + 1i*zbbbs(1:n)-zmaxis*1i);
      
      if ~any(thbbbs(1:n) > 0 & thbbbs(1:n) < +pi) | ...
         ~any(thbbbs(1:n) < 0 & thbbbs(1:n) > -pi)
        trace_status = 1;
      end
    
    end

    if trace_status == 0
   
      % Order the points according to theta from -pi to +pi
      [~, i] = min(thbbbs(1:n));
      if sum(diff(thbbbs(1:n)) > 0) > sum(diff(thbbbs(1:n)) < 0)
	thbbbs(1:n)  = thbbbs([i:n, 1:i-1]);
	gbbbs(1:n) = gbbbs([i:n, 1:i-1]);
	ibbbs(1:n) = ibbbs([i:n, 1:i-1]);
	xx(1:n) = xx([i:n, 1:i-1]);
      else
	thbbbs(1:n)  = thbbbs([i:-1:1, n:-1:i+1]);
	gbbbs(1:n) = gbbbs([i:-1:1, n:-1:i+1]);
	ibbbs(1:n) = ibbbs([i:-1:1, n:-1:i+1]);
	xx(1:n) = xx([i:-1:1, n:-1:i+1]);
      end

      jr = gbbbs == nz;
      jz = gbbbs == 1;
      rbbbs(jr) = rgg(ibbbs(jr))+xx(jr)*dr;
      rbbbs(jz) = rgg(ibbbs(jz));
      zbbbs(jr) = zgg(ibbbs(jr));
      zbbbs(jz) = zgg(ibbbs(jz))+xx(jz)*dz;
      rbbbs(gbbbs==0) = rbdef;
      zbbbs(gbbbs==0) = zbdef;

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
      
      % Squeeze in x-point wannabes
      for j = 1:xwannabe.count
        r = xwannabe.r(j);
        z = xwannabe.z(j);
	drdpsi = xwannabe.drdpsi(j,:);
	dzdpsi = xwannabe.dzdpsi(j,:);
	rhox = sqrt((r-rmaxis)^2+(z-zmaxis)^2);      
        thx  = angle(r-rmaxis +1i*z-zmaxis*1i);
	[dum,i] = min((thbbbs(2:nbbbs)-thx).*(thbbbs(1:nbbbs-1)-thx));
	i = i+1;
	nbbbs = nbbbs+1;
	gbbbs(i+1:nbbbs) = gbbbs(i:nbbbs-1);
	ibbbs(i+1:nbbbs) = ibbbs(i:nbbbs-1);
	rbbbs(i+1:nbbbs) = rbbbs(i:nbbbs-1);
	zbbbs(i+1:nbbbs) = zbbbs(i:nbbbs-1);
	wbbbs(i+1:nbbbs,:) = wbbbs(i:nbbbs-1,:);
	thbbbs(i+1:nbbbs) = thbbbs(i:nbbbs-1);
	rhobbbs(i+1:nbbbs) = rhobbbs(i:nbbbs-1);
	dpsibbbsdr(i+1:nbbbs) = dpsibbbsdr(i:nbbbs-1);
	dpsibbbsdz(i+1:nbbbs) = dpsibbbsdz(i:nbbbs-1);
	gbbbs(i) = j/10;
	thbbbs(i) = thx;
	dth = thbbbs(i+1)-thbbbs(i-1);
	rhobbbs(i) = (thbbbs(i+1)-thx)/dth*rhobbbs(i-1) + ...
	             (thx-thbbbs(i-1))/dth*rhobbbs(i+1);
	sinth = sin(thx);
	costh = cos(thx);
	rbbbs(i) = rmaxis + rhobbbs(i)*costh;
	zbbbs(i) = zmaxis + rhobbbs(i)*sinth;
	p = 0;
	drhob = realmax;
	while abs(drhob) > 1e-12 & p < 9
	  p = p+1;
	  kr0 = floor((rbbbs(i)-rg(1))/dr);
	  kz1 = floor((zbbbs(i)-zg(1))/dz)+1;
	  if kr0 < 1 | kr0 > nr-3 | kz1 < 2 | kz1 > nz-2
	    trace_status = 1;
	    break
	  end
	  k = kr0*nz+kz1;
	  pp = psizr(k+neighbors');
	  tr = (rbbbs(i)-rgg(k))/dr;
	  tz = (zbbbs(i)-zgg(k))/dz;
	  wr0 = [1 tr tr^2 tr^3]*mx;
	  wz0 = mx'*[1 tz tz^2 tz^3]';
	  wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
	  wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
	  dpsibbbsdr(i) = reshape(wz0*wr1,1,16)*pp;
	  dpsibbbsdz(i) = reshape(wz1*wr0,1,16)*pp;
	  dum = (dpsibbbsdr(i)*costh+dpsibbbsdz(i)*sinth);
	  drhob = (psibry-reshape(wz0*wr0,1,16)*pp)/dum;
	  if abs(drhob) > dr/2^p
	    drhob = drhob/abs(drhob)*dr/2^p;
	  end
	  rhobbbs(i) = rhobbbs(i) + drhob;
	  rbbbs(i) = rmaxis + rhobbbs(i)*costh;
	  zbbbs(i) = zmaxis + rhobbbs(i)*sinth;
	end
	ibbbs(i) = k;
	dum16 = (drdpsi*sinth-dzdpsi*costh)/dum*rhobbbs(i)/rhox;
	drbdpsix(j,:) =  dum16*dpsibbbsdz(i);
	dzbdpsix(j,:) = -dum16*dpsibbbsdr(i);
	dum16 = (drmaxisdpsi*sinth-dzmaxisdpsi*costh)/dum*(rhox-rhobbbs(i))/rhox;
	drbdpsia(j,:) =  dum16*dpsibbbsdz(i);
	dzbdpsia(j,:) = -dum16*dpsibbbsdr(i);
	drbdpsib(j,:) = costh/dum*wb;
	dzbdpsib(j,:) = sinth/dum*wb;
	wp = reshape(wz0*wr0,1,16);
	drbdpsip(j,:) = -costh/dum*wp;
	dzbdpsip(j,:) = -sinth/dum*wp;
      end
      
      
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
  end
end

