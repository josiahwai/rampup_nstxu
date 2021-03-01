%  USAGE:   gs_connection_length
%
%  PURPOSE: Calculate connection lengths from all grid points inside limiter
%
%  INPUTS: psizr, flux at grid points rg, zg
%
%  OUTPUTS: lconnfzr, connection lengths from grid points to limiter along phi
%           lconnbzr, connection lengths going opposite to phi direction
%
%  METHOD:  countourc finds contours for fluxes psizr inside limiter
	
%  NOTES:  
	
%  WRITTEN BY:  Anders Welander  ON	2015-06-18
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

brzr = -[zeros(1,nr);psizr(3:nz,:)-psizr(1:nz-2,:);zeros(1,nr)]./rgg/(2*pi);
bzzr = +[zeros(nz,1) psizr(:,3:nr)-psizr(:,1:nr-2) zeros(nz,1)]./rgg/(2*pi);
%Babszr = sqrt((psizr_r./rgg/twopi).^2 + (psizr_z./rgg/twopi).^2 + (fpolzr./rgg).^2);
brzr = -psizr_z./rgg/twopi;
bzzr = +psizr_z./rgg/twopi;
bpzr = sqrt(brzr.^2+bzzr.^2);
btzr = fpolzr./rgg;
nzr = zeros(nz,nr);
mzr = zeros(nz,nr);
reszr = zeros(nz,nr);
dk = [-nz-1 -nz -nz+1 1 nz+1 nz nz-1 -1 -nz-1];
igg = isinpoly(rgg,zgg);
hzr = zeros(nz,nr); % Flags for tracing
vzr = zeros(nz,nr); % Flags for tracing
lconnbzr = zeros(nz,nr);
lconnfzr = zeros(nz,nr);
lconnbzr(iplasma) = inf;
lconnfzr(iplasma) = inf;
clf
pb
hold on
for k = 1:ngg
  if igg(k) & ~iplasma(k)
    % trace contour until it closes on itself or exits vessel
    R1 = rgg(k);
    Z1 = zgg(k);
    BR1 = brzr(k);
    BZ1 = bzzr(k);
    BP1 = bpzr(k);
    BT1 = btzr(k);
    PSI = psizr(k);
    hzr(k) = k; % Interval rgg(k) <= r < rgg(k+nz) has been visited
    vzr(k) = k; % Interval zgg(k) <= z < zgg(k+1) visited (for flux psizr(k))
    for l = 1:8
      k1 = k+dk(l);
      k2 = k+dk(l+1);
      dpsi1 = psizr(k1)-PSI;
      dpsi2 = psizr(k2)-PSI;
      if dpsi1*dpsi2 <= 0
        dum = dpsi1/(psizr(k1)-psizr(k2));
	R2 = (1-dum)*rgg(k1) + dum*rgg(k2);
	Z2 = (1-dum)*zgg(k1) + dum*zgg(k2);
	k3 = min(k1,k2);
	if dum == 0
	  k3 = k1;
	  hzr(k3) = k;
	  vzr(k3) = k;
	elseif dum == 1
	  k3 = k2;
	  hzr(k3) = k;
	  vzr(k3) = k;
	elseif rgg(k1) == rgg(k2)
	  hzr(k3) = k;
	elseif zgg(k1) == zgg(k2)
	  vzr(k3) = k;
	end
	plot([R1 R2],[Z1 Z2])
	BP2 = (1-dum)*bpzr(k1) + dum*bpzr(k2);
	BT2 = (1-dum)*btzr(k1) + dum*btzr(k2);
	if [R2-R1 Z2-Z1]*[BR1;BZ1]*BT1 > 0
	  % Tracing field line in positive B-direction
	  mzr(i,j) = mzr(i,j)+1;
	  % Major radius in middle of adjacent contour points
	  R = (R1+R2)/2;
	  % Poloidal B-field in middle of adjacent contour points
	  BP = (BP1+BP2)/2;
	  % Toroidal B-field in middle of adjacent contour points
	  BT = (BT1+BT2)/2;
	  % Total B-field in middle of adjacent contour points
	  BTOT = sqrt(BP^2+BT^2);
	  % Poloidal distance between contour points
	  DPOL = sqrt((R2-R1)^2+(Z2-Z1)^2);
	  % Toroidal distance between contour points
	  DTOR = abs(BT/BP*DPOL);
	  % Toroidal angle between contour points
	  DPHI = DTOR/R;
	  % Total distance between contour points
	  DLEN = sqrt(DTOR^2+DPOL^2);
        else
	  nzr(i,j) = nzr(i,j)+1;
	end
      end
    end
  end
end

return

% Magnitude of poloidal B-field in middle of adjacent contour points
bpolc = (bpcont(1:npola-1,:)+bpcont(2:npola,:))/2;

% Toroidal B-field in middle of adjacent contour points
btorc = (btcont(1:npola-1,:)+btcont(2:npola,:))/2;

% Total B-field in middle of adjacent contour points
btotc = sqrt(bpolc.^2+btorc.^2);

% Toroidal distance between contour points
dtorc = abs(btorc./bpolc.*dpolc);
dtorc(:,1) = 2*pi*rmaxis*qc(1)/size(dtorc,1);

% Toroidal angle between contour points
dphic = dtorc./(rcont(1:npola-1,:)+rcont(2:npola,:))*2;

% Total distance between contour points
dlenc = sqrt(dtorc.^2+dpolc.^2);

% Parallel current density along B*sign(bzero)
jparc = (ones(npola-1,1)*fprimc'/mu0.*btotc + ...
         ones(npola-1,1)*(fpolc.*pprimec)'./btotc)*sign(bzero);

% Resistive voltage along B*sign(bzero) between contour points
dvolt = dlenc.*jparc.*(etacont(1:npola-1,:)+etacont(2:npola,:))/2;

% Average resistive voltage along B*sign(bzero) through one toroidal turn
vresc = sum(dvolt)'./qc;
for i = 1:ncont
  if psibarc(i) == 0
    vresc(i) = 2*pi*rmaxis*mean(etacont(:,i))*mean(jparc(:,i));
  end
  if psibarc(i) == 1 && lim == 0
    jbdef = rbdef*pprime(nr) + ffprim(nr)/rbdef/mu0;
    vresc(i) = 2*pi*rbdef*etacont(1,i)*jbdef;    
  end
end

if isfield(index_in_y,'vres')
  lae.y(index_in_y.vres) = vresc;
end

if isfield(index_in_y,'jpar')
  lae.y(index_in_y.jpar) = sum(dlenc.*jparc)./sum(dlenc);
end

% Inductive voltage at stationary rcont, zcont as function of xdot
dvpdxdot = zeros(ncont,nx); % Voltage from poloidal flux
dvtdxdot = zeros(ncont,nx); % Voltage from toroidal flux
dhalffpolsquaredcdsf = c0(iknotc,:) + ...
  c1(iknotc,:).*(psibarc(:)*ones(1,nkn+2)) + ...
  c2(iknotc,:).*(psibarc(:).^2*ones(1,nkn+2)) +...
  c3(iknotc,:).*(psibarc(:).^3*ones(1,nkn+2));
dfpolcdsf = sign(rzero*bzero)./sqrt(2*halffpolsquaredc)*...
  ones(1,nkn+2).*dhalffpolsquaredcdsf;
for j = 1:ncont
  if psibarc(j) == 0
    dvpdxdot(j,:) = qc(j)*dpsimagdx;
  elseif psibarc(j) == 1 & lim == 0
    dvpdxdot(j,:) = qc(j)*dpsibrydx;
  else
    n = npola-1; % index to the dphic before point i
    for i = 1:npola-1 % index to points on a contour
      dum = (dphic(n,j)+dphic(i,j))/4/pi; % toroidal fraction for i
      n = i; % now updated to index the dphic before next i
      % Poloidal flux contributions to the voltage
      k = i+(j-1)*npola; % index into collapsed matrix
      dvpdxdot(j,:) = dvpdxdot(j,:) + ...
	dum*wcont(k,:)*dpsizrdx(icont(k,:),:);
    end
    % Toroidal flux contributions to the voltage
    %Tc = 0.50*[0; cumsum((fpolc(2:end)+fpolc(1:end-1)).*diff(Lc))];
    if j > 1
      dvtdxdot(j,indsf) = dvtdxdot(j-1,indsf) + ...
        (Lc(j)-Lc(j-1))/2*(dfpolcdsf(j,:)+dfpolcdsf(j-1,:));
    end
  end
end
% vindc = average of induced voltage per toroidal turn on contours
dvindcdxdot = (dvpdxdot + dvtdxdot)./(qc*ones(1,nx));

% Project onto nkn+2 equations
%Projcurprofeqs = pinv(dvindcdxdot(:,indsf));
Projcurprofeqs = zeros(nkn+2,nr);
Projcurprofeqs(1,1) = 1;
for i = 2:nkn+1
  j = 2;
  while psibarc(j) < psikn(i)
    Projcurprofeqs(i,j-1) = Projcurprofeqs(i,j-1) + Ac(j)-Ac(j-1);
    Projcurprofeqs(i,j) = Projcurprofeqs(i,j) + Ac(j)-Ac(j-1);
    j = j+1;
  end
  dum = (psikn(i)-psibarc(j-1))/(psibarc(j)-psibarc(j-1));
  Projcurprofeqs(i,j-1) = Projcurprofeqs(i,j-1) + dum*(Ac(j)-Ac(j-1));
  Projcurprofeqs(i,j) = Projcurprofeqs(i,j) + dum*(Ac(j)-Ac(j-1));
end
Projcurprofeqs(nkn+2,nr) = 1;

return

% Confirmation plot
gs_profiles
etaav = spline(psibarc,etacont(1,:),psibar);
plot(psibarc,vresc,psibar,etaav(:).*jtav(:)*2*pi*rmaxis)

Pnkn = zeros(nkn+2,nr);
for i = 1:nkn+1
  k = 1;
  for j = 2:nr-1
    if psibarc(j) < psikn(i)
      k = j;
    end
  end
  x = (psikn(i)-psibarc(k))/(psibarc(k+1)-psibarc(k));
  Pnkn(i,k+1) = x;
  Pnkn(i,k) = 1-x;
  if i ~= nkn+1
    %Pnkn(i,1:k-1) = 1;
  else
    Pnkn(i,:) = 0;
    Pnkn(i,nr) = 1;
  end
end
Pnkn(nkn+2,:) = 1;
if 0
  Pnkn = zeros(nkn+2,nr);
  Pnkn(1,1) = 1;
  Pnkn(nkn+2,nr) = 1;
  for i = 2:nkn+1
    k = 1;
    Pnkn(i,1) = Ac(2)/2;
    for j = 2:nr-2
      if psibarc(j) < psikn(i)
	Pnkn(i,j) = (Ac(j+1)-Ac(j-1))/2;
	k = j;
      end
    end
    x = (psikn(i)-psibarc(k))/(psibarc(k+1)-psibarc(k));
    Pnkn(i,k+1) = x*(Ac(k+2)-Ac(k))/2;
    Pnkn(i,k) = (1-x)*(Ac(k+1)-Ac(k-1))/2;
  end
end
xdot = zeros(nx,1);
xdot(indsf) = pinv(Pnkn*dvindcdxdot(:,indsf))*(Pnkn*vresc);
%xdot(1:nx-1) = pinv(Pnkn*dvindcdxdot(:,1:nx-1))*(Pnkn*vresc);
%xdot(1:nx-1) = pinv(dvindcdxdot(:,1:nx-1))*(vresc);
xdot(indsf) = pinv(Projcurprofeqs*dvindcdxdot(:,indsf))*...
  (Projcurprofeqs*vresc);
xdot(1:nx-1) = pinv(Projcurprofeqs*dvindcdxdot(:,1:nx-1))*...
  (Projcurprofeqs*vresc);
vp=dvpdxdot*xdot;
vt=dvtdxdot*xdot;
clf
hold on
plot(psibar,vresc,psibar,vp,psibar,vp+vt)
a = axis;
for i = 1:nkn+1
  plot([psikn(i) psikn(i)],a(3:4),'k--')
end
hold off
psilim =  sum(wl.*psizr(iil),2); % Flux at limiter
psidlim = sum(wld.*psizr(iil),2);
psiblim = sum(wlb.*psizr(iil),2);
psitlim = sum(wlt.*psizr(iil),2);
drzl = turnin*[diff(rl); diff(zl)];
% To confirm drzl
% clf,plot(rl,zl),hold
% for i=1:nl-1,plot(rl(i)+[0 drzl(1,i)],zl(i)+[0 drzl(2,i)]),end,axis image
ltflag = zeros(nl,1);
ltflag(touches.limiter(1:touches.count)) = 1;
% Field lines are traced in only one phi direction
clf
ncont = 0;
pb
hold on
plot(rl,zl,'rx','linew',2)
i1 = nl-1;
for i2 = 1:nl-1
  i3 = i2+1;
  if psidlim(i2) > 0
    trace_status = -1;
    xh = nan(ngg,3); % Solutions between horizontal neighbors
    xv = nan(ngg,3); % Solutions between vertical neighbors
    n = 1;
    j = max(2,min(nr-2,floor((rl(i2)-rg(1))/dr)+1));
    i = max(2,min(nz-2,floor((zl(i2)-zg(1))/dz)+1));
    k = i+(j-1)*nz;
    if isnan(xh(k,1))
      xh(k,:) = gs_solve_hermite_cubic(psizr(i,j-1:j+2), psilim(i2));
    end
    if isnan(xv(k,1))
      xv(k,:) = gs_solve_hermite_cubic(psizr(i-1:i+2,j), psilim(i2));
    end
    if isnan(xh(k+1,1))
      xh(k+1,:) = gs_solve_hermite_cubic(psizr(i+1,j-1:j+2), psilim(i2));
    end
    if isnan(xv(k+nz,1))
      xv(k+nz,:) = gs_solve_hermite_cubic(psizr(i-1:i+2,j+1), psilim(i2));
    end
    fh = logical(ones(ngg,3)); % Flags for okay to visit
    fv = logical(ones(ngg,3)); % Flags for okay to visit
    % Trace from touch point to edge of cell
    done_tracing = false;
    irect = 1;
    ir = (rl(i2)-rg(j))/dr;
    iz = (zl(i2)-zg(i))/dz;
    ir0 = ir;
    iz0 = iz;
    counter = 0;
    while ~done_tracing
      counter = counter+1;
      if counter > 4
	trace_status = 1;
	done_tracing = true;
      end
      if irect == 1
	[ir, iz, m] = gs_trace_step(...
	  ir, iz, 0, 1, 0, iz0, psizr(i-1:i+2,j-1:j+2), psilim(i2));
	if m == 1
          irect = 2;
	else
          done_tracing = true;
	end
      else % irect = 2
	[ir, iz, m] = gs_trace_step(...
	  ir, iz, 0, 1, iz0, 1, psizr(i-1:i+2,j-1:j+2), psilim(i2));
	if m == -1
          irect = 1;
	else
          done_tracing = true;
	end
      end
    end
    if m == 0 % Try chopping the cell with vertical cut
      done_tracing = false;
      counter = 0;
      ir = ir0;
      iz = iz0;
      while ~done_tracing
	counter = counter+1;
	if counter > 4
	  trace_status = 1;
	  done_tracing = true;
	end
	if irect == 1
	  [ir, iz, m] = gs_trace_step(...
	    ir, iz, 0, ir0, 0, 1, psizr(i-1:i+2,j-1:j+2), psilim(i2));
	  if m == 4
            irect = 2;
	  else
            done_tracing = true;
	  end
	else % irect = 2
	  [ir, iz, m] = gs_trace_step(...
	    ir, iz, ir0, 1, 0, 1, psizr(i-1:i+2,j-1:j+2), psilim(i2));
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
    if p < 4
      fh(k,p) = 0;
      i = i-1;
      r = rg(j) + xh(k,p)*dr;
      z = zg(i+1);
    elseif p < 7
      fv(k,p-3) = 0;
      j = j-1;
      r = rg(j+1);
      z = zg(i) + xv(k,p-3)*dz;
    elseif p < 10
      i = i+1;
      fh(k+1,p-6) = 0;
      r = rg(j) + xh(k+1,p-6)*dr;
      z = zg(i-1);
    else
      j = j+1;
      fv(k+nz,p-9) = 0;
      r = rg(j-1);
      z = zg(i) + xv(k+nz,p-9)*dz;
    end
    plot([rl(i2) r],[zl(i2) z])
  end
  i1 = i2;
end






return
