%  USAGE:   gs_helical_voltage
%
%  PURPOSE: Calculate voltage along B*sign(bzero) for contours
%
%  INPUTS: etacont, parallel resistivity at points rcont, zcont
%          Create rcont, zcont and other inputs with:
%            gs_eq_analysis
%            gs_response
%            gs_trace_contours
%            gs_contour_response
%            gs_contour_profiles
%          Then create etacont before calling this script
%
%  OUTPUTS: vresc, resistive voltage per toroidal turn
%           dvindcdxdot, inductive voltage as function of xdot
%	    Projcurprofeqs, projection onto nkn+2 equations
%
%  METHOD:  The contours are stationary. The outputs can be used
%           to calculate surface-integrals of dB/dt + curl E = 0
%           for surfaces bounded by field lines at fluxes psibarc,
%           using that E_para = eta_para*j_para along field lines
	
%  NOTES:  
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	2015-01-04
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
dlenc = sqrt(drcont.^2+dzcont.^2+dtorc.^2);
lae.dlenc = dlenc;

% Components of unit vectors in field direction between contour points
urc = drcont./dlenc;
uzc = dzcont./dlenc;
utc = dtorc./dlenc;

% fprim, ffprim and pprime at contour points
fprimcont = ones(npola,1)*fprimc';
ffprimcont = ones(npola,1)*ffprimc';
pprimecont = ones(npola,1)*pprimec';

% Components of current density at contour points
jrcont = fprimcont.*brcont/mu0;
jzcont = fprimcont.*bzcont/mu0;
jtcont = rcont.*pprimecont+ffprimcont./rcont/mu0;

% Components of current density between contour points
jrc = (jrcont(1:end-1,:)+jrcont(2:end,:))/2;
jzc = (jzcont(1:end-1,:)+jzcont(2:end,:))/2;
jtc = (jtcont(1:end-1,:)+jtcont(2:end,:))/2;

% Parallel current density along B*sign(bzero)
% jparc2 = (ones(npola-1,1)*fprimc'/mu0.*btotc + ...
%          ones(npola-1,1)*(fpolc.*pprimec)'./btotc)*sign(bzero);
jparc = urc.*jrc + uzc.*jzc + utc.*jtc; % jparc2 is same thing
lae.jparc = jparc;
lae.jbdef = rbdef*pprime(nr) + ffprim(nr)/rbdef/mu0;

% Parallel resistivity
if exist('u','var') & isfield(u,'etacont') % User sent eta for contours
  etacont = u.etacont;
else % Calculate eta versus psibar from eta vs rhot
  if exist('u','var') &  isfield(u,'eta_vs_rhot')
    eta_vs_rhot = u.eta_vs_rhot;
  end
  if ~exist('eta_vs_rhot') | isempty(eta_vs_rhot)
    eta_vs_rhot = zeros(nr,1);
  end
  rhoe = linspace(0,1,length(eta_vs_rhot))';
  eta = spline(rhoe, eta_vs_rhot, rhot);
  etacont = ones(npola,1)*eta';
end
% Resistivity between contour points
etac = (etacont(1:npola-1,:)+etacont(2:npola,:))/2;
lae.etac = etac;

% How jpar, vresc change AT FIXED rcont, zcont with x
djpardx = zeros(ncont,nx);
dvrescFdx = zeros(ncont,nx);
for j = 1:ncont
  dffprimcdxj = dffprimcdpsimag(j)*dpsimagdx+dffprimcdpsibry(j)*dpsibrydx;
  dpprimecdxj = dpprimecdpsimag(j)*dpsimagdx+dpprimecdpsibry(j)*dpsibrydx;
  k = 1+(j-1)*npola;
  dumnx = psibarc(j)*dpsibrydx + (1-psibarc(j))*dpsimagdx;
  dpsibarcdx1 = wcont(k,:)*dpsizrdx(icont(k,:),:) - dumnx;
  dfp1dx = dfprimcdpsibarc(j)*dpsibarcdx1;
  dff1dx = dffprimcdpsibarc(j)*dpsibarcdx1+dffprimcdxj;
  dpp1dx = dpprimecdpsibarc(j)*dpsibarcdx1+dpprimecdxj;
  dbr1dx = -wcont_z(k,:)*dpsizrdx(icont(k,:),:)/rcont(1,j)/twopi;
  dbz1dx = +wcont_r(k,:)*dpsizrdx(icont(k,:),:)/rcont(1,j)/twopi;
  djr1dx = fprimc(j)/mu0*dbr1dx + brcont(k)/mu0*dfp1dx;
  djz1dx = fprimc(j)/mu0*dbz1dx + bzcont(k)/mu0*dfp1dx;
  djt1dx = rcont(k)*dpp1dx + dff1dx/(mu0*rcont(k));
  ltot = sum(dlenc(:,j));
  for i = 2:npola
    k = i+(j-1)*npola;
    dpsibarcdx2 = wcont(k,:)*dpsizrdx(icont(k,:),:) - dumnx;
    dfp2dx = dfprimcdpsibarc(j)*dpsibarcdx2;
    dff2dx = dffprimcdpsibarc(j)*dpsibarcdx2+dffprimcdxj;
    dpp2dx = dpprimecdpsibarc(j)*dpsibarcdx2+dpprimecdxj;
    dbr2dx = -wcont_z(k,:)*dpsizrdx(icont(k,:),:)/rcont(1,j)/twopi;
    dbz2dx = +wcont_r(k,:)*dpsizrdx(icont(k,:),:)/rcont(1,j)/twopi;
    djr2dx = fprimc(j)/mu0*dbr2dx + brcont(k)/mu0*dfp2dx;
    djz2dx = fprimc(j)/mu0*dbz2dx + bzcont(k)/mu0*dfp2dx;
    djt2dx = rcont(k)*dpp2dx + dff2dx/(mu0*rcont(k));
    kc = i-1+(j-1)*(npola-1);
    djpar1dx = urc(kc)/2*(djr1dx+djr2dx) + uzc(kc)/2*(djz1dx+djz2dx) + ...
               utc(kc)/2*(djt1dx+djt2dx);
    dvrescFdx(j,:) = dvrescFdx(j,:) + etac(kc)*dlenc(kc)/qc(j)*djpar1dx;
    djpardx(j,:) = djpardx(j,:) + dlenc(kc)/ltot*djpar1dx;
  end
  djr1dx = djr2dx;
  djz1dx = djz2dx;
  djt1dx = djt2dx;
end
dvrescFdetac = dlenc.*jparc./(ones(npola-1,1)*qc');

if isfield(index_in_y,'jpar')
  lae.y(index_in_y.jpar) = sum(dlenc.*jparc)./sum(dlenc);
  dydx(index_in_y.jpar,:) = djpardx;
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

return

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
Projcurprofeqs(nkn+2,nr-1) = 1;

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
