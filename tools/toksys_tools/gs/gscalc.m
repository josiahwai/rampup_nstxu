function [psizr_plas, brs, bzs] = gscalc(jphis,rg,zg,mpp,nrb,nzb)
%
%  USAGE: [psizr_pla, br, bz] = gscalc(jphi,rg,zg,mpp,nrb,nzb)
%
%  PURPOSE: Calculate contribution to flux and fields from plasma
%
%  INPUTS: jphi, current density in plasma [MA/m2]
%            rg, vector of R coordinates for the grid [m]
%            zg, vector of Z coordinates for the grid [m]
%           mpp, (optional) supply (nr*nz,nr) mutuals [H] if available.
%       nrb,nzb, (optional) number of boundary calculations per edge.
%                Values from 2 to [nr,nz]. Default: sqrt([nr,nz]).
%
%  OUTPUTS: psizr_pla, plasma flux [Wb]
%                  br, plasma radial field [T]
%                  bz, plasma vertical field [T]
%
%  METHOD: The Grad-Shafranov equation is solved on finite-difference form
%          with boundary calculated by mutual inductances.
%          Tne FD equations are solved by sine transforms in Z
%          and tridiagonal matrix solutions in R.

%	
%  VERSION @(#)gscalc.m	1.5 04/11/13
%
%  WRITTEN BY:  Anders Welander  ON	6/11/11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mu0 = .4e-6*pi;
  
  nr = length(rg);
  nz = length(zg);
  ngg = nr*nz;
  rg = rg(:)';
  rgg = ones(nz,1)*rg;
  dr = (rg(nr)-rg(1))/(nr-1);
  dz = (zg(nz)-zg(1))/(nz-1);  
  alpha = dr/dz;
  beta = dz/dr;
  
  [n1,n2] = size(jphis);
  if n1*n2 == ngg % That means we only have 1 applied flux and it is nz*nr rather than ngg*1
    jphis = jphis(:); % Collapse it to ngg*1
    nfluxes = 1;
  else
    nfluxes = n2;
  end
  
  for iflux = 1:nfluxes
  
    jphi = reshape(jphis(:,iflux),nz,nr);

    % Calculate the boundary condition (i.e. flux at edge of grid)
    if exist('mpp','var') & size(mpp,1) == ngg % Simply use mutuals to figure out boundary flux
      psibot = 1e6*dr*dz*jphi(:)'*mpp(:,2:nr-1);
      psitop = 1e6*dr*dz*jphi((1:ngg)+mod(nz-(1:ngg),nz)-mod(0:ngg-1,nz))*mpp(:,2:nr-1);
      for j = 1:nz
	izshift = (1+abs(j-(1:nz))-(1:nz))'*ones(1,nr);
	psiinn(j,1) = 1e6*dr*dz*jphi(:)'*mpp((1:ngg)'+izshift(:),1);
	psiout(j,1) = 1e6*dr*dz*jphi(:)'*mpp((1:ngg)'+izshift(:),end);
      end
    else % Calculate boundary using fine1 
      % Decide interval between full calculations of boundary points
      if exist('nrb','var')
	nrb = max(min(nr-1,nrb),2);
      else
	nrb = round(sqrt(nr));
      end
      if exist('nzb','var')
	nzb = max(min(nz-1,nzb),2);
      else
	nzb = round(sqrt(nz));
      end
      % Various variables needed for fine1
      zgg = zg(:)*ones(1,nr);
      jphi([1 nz],:) = 0; jphi(:,[1 nr]) = 0; % To avoid risk of nasty error in fine1
      iplasma = find(jphi);
      ri = rgg(iplasma)-dr/2; ro = rgg(iplasma)+dr/2;
      zl = zgg(iplasma)-dz/2; zu = zgg(iplasma)+dz/2;
      ig = 1e6*dr*dz*mu0*jphi(iplasma); zs = iplasma*0;    
      % Define boundary coordinates
      rb = linspace(rg(1),rg(nr),nrb+1); drb = rb(2)-rb(1);
      zb = linspace(zg(1),zg(nz),nzb+1); dzb = zb(2)-zb(1);
      % Coordinates for a round trip beginning at rmin, zmin then up, out, down, in
      rbb = [rg(1)+zeros(1,nzb) rb rg(nr)+zeros(1,nzb) rb(nrb:-1:1)];
      zbb = [zb zg(nz)+zeros(1,nrb) zb(nzb:-1:1) zg(1)+zeros(1,nrb)];
      nbb = length(rbb);
      % Take a round trip to measure psi, br, bz along grid edge, around the plasma
      for j=1:nbb-1
	[hr,hz,fl] = fine1(ri,ro,zl,zu,ig,rbb(j)+zs,zbb(j)+zs);
	psi(j) = sum(fl);
	yr(j) = +2*pi*rbb(j)*sum(hz);
	yz(j) = -2*pi*rbb(j)*sum(hr);
      end
      psi(nbb) = psi(1); yr(nbb) = yr(1); yz(nbb) = yz(1);
      [hr,hz,fl] = fine1(ri,ro,zl,zu,ig,rbb(1)+zs,zbb(1)-1e-4+zs);
      yzz = (yz(1)+2*pi*rbb(1)*sum(hr))/1e-4;
      % Determine the first 4 poly coeffs assuming the intervals are from x = 0 to 1
      k0 = 1:nzb; k1 = nzb+(1:nrb); k2 = nzb+nrb+k0; k3 = nzb+nrb+k1;
      db = [dzb+zeros(1,nzb) drb+zeros(1,nrb) dzb+zeros(1,nzb) drb+zeros(1,nrb)];
      y1 = psi(2:end);
      yd1 = [yz(k0+1)*dzb yr(k1+1)*drb -yz(k2+1)*dzb -yr(k3+1)*drb];
      p0 = psi(1:end-1);
      p1 = [yz(k0)*dzb yr(k1)*drb -yz(k2)*dzb -yr(k3)*drb];
      p2 = p1-yd1+3*(y1-p1-p0);
      p3 = y1-p0-p1-p2;    
      % Make second derivative continuous at break points and delstar(psi)=0 at corners
      % Add p4*x^2*(1-x)^2. Second derivative of this term is 2*p4 at both x=0 and x=1
      % The equation for p4 will be: p4 = c0+c1*p4(end). We will find c0 and c1 first and p4(end) at the end
      % delstar(psi) for rmin, zmin is: (polynomials are for -r direction in bottom layer)
      % 2*p2(1)/dzb^2+2*p4(1)/dzb^2+2*p2(end)/drb^2+6*p3(end)/drb^2+2*p4(end)/drb^2-yr(end)/rb(1) = 0
      c0 = -(2*p2(1)/dzb^2+2*p2(end)/drb^2+6*p3(end)/drb^2-yr(end)/rb(1))*dzb^2/2;
      c1 = -dzb^2/drb^2;
      % Going up on the inside, p2(j-1)+3*p3(j-1)+p4(j-1) = p2(j)+p4(j)
      for j = 2:nzb
	c0(j) = c0(j-1) + p2(j-1)+3*p3(j-1)-p2(j);
	c1(j) = c1(j-1);
      end
      % delstar(psi) at rmin, zmax: j=nzb
      % 2*p2(j)/dzb^2+6*p3(j)/dzb^2+2*p4(j)/dzb^2 + 2*p2(j+1)/drb^2+2*p4(j+1)/drb^2-yr(j+1)/rb(1) = 0
      c0(j+1) = -(2*p2(j)/dzb^2+6*p3(j)/dzb^2+2*c0(j)/dzb^2+2*p2(j+1)/drb^2-yr(j+1)/rb(1))*drb^2/2;
      c1(j+1) = -c1(j)/dzb^2*drb^2;
      % Going out on the upside, p2(k-1)+3*p3(k-1)+p4(k-1) = p2(k)+p4(k)
      for k = nzb+(2:nrb)
	c0(k) = c0(k-1) + p2(k-1)+3*p3(k-1)-p2(k);
	c1(k) = c1(k-1);
      end
      % delstar(psi) at rmax, zmax: k=nzb+nrb, the new poly is for -z direction
      % 2*p2(k)/drb^2+6*p3(k)/drb^2+2*p4(k)/drb^2-yr(k+1)/rg(nr) + 2*p2(k+1)/dzb^2+2*p4(k+1)/dzb^2 = 0
      c0(k+1) = -(2*p2(k)/drb^2+6*p3(k)/drb^2+2*c0(k)/drb^2-yr(k+1)/rg(nr)+2*p2(k+1)/dzb^2)*dzb^2/2;
      c1(k+1) = -c1(k)/drb^2*dzb^2;
      % Going down on the outside
      for j = nzb+nrb+(2:nzb)
	c0(j) = c0(j-1) + p2(j-1)+3*p3(j-1)-p2(j);
	c1(j) = c1(j-1);
      end
      % delstar(psi) at rmax, zmin: j=nzb+nrb+nzb
      % 2*p2(j)/dzb^2+6*p3(j)/dzb^2+2*p4(j)/dzb^2 + 2*p2(j+1)/drb^2+2*p4(j+1)/drb^2-yr(j+1)/rg(nr) = 0
      c0(j+1) = -(2*p2(j)/dzb^2+6*p3(j)/dzb^2+2*c0(j)/dzb^2+2*p2(j+1)/drb^2-yr(j+1)/rg(nr))*drb^2/2;
      c1(j+1) = -c1(j)/dzb^2*drb^2;
      % Going in on the downside
      for k = nzb+nrb+nzb+(2:nrb)
	c0(k) = c0(k-1) + p2(k-1)+3*p3(k-1)-p2(k);
	c1(k) = c1(k-1);
      end
      p4end = ((yzz-2*p2(1)/dzb^2)*dzb^2/2-c0(1))/c1(1);
      p4 = c0+c1*p4end;
      pp = mkpp(zb,           [p4(k0)/dzb^4; (p3(k0)-2*p4(k0))/dzb^3; (p2(k0)+p4(k0))/dzb^2; p1(k0)/dzb; p0(k0);]');
      psiinn = ppval(pp,zg);
      pp = mkpp(rb,           [p4(k1)/drb^4; (p3(k1)-2*p4(k1))/drb^3; (p2(k1)+p4(k1))/drb^2; p1(k1)/drb; p0(k1);]');
      psitop = ppval(pp,rg(2:nr-1));
      pp = mkpp(-zb(end:-1:1),[p4(k2)/dzb^4; (p3(k2)-2*p4(k2))/dzb^3; (p2(k2)+p4(k2))/dzb^2; p1(k2)/dzb; p0(k2);]');
      psiout = ppval(pp,-zg);
      pp = mkpp(-rb(end:-1:1),[p4(k3)/drb^4; (p3(k3)-2*p4(k3))/drb^3; (p2(k3)+p4(k3))/drb^2; p1(k3)/drb; p0(k3);]');
      psibot = ppval(pp,-rg(2:nr-1));    
    end % Done calculating the boundary condition

    % Right hand side of finite-difference equation
    source = 1e6*dr*dz*mu0*rgg(2:nz-1,2:nr-1).*jphi(2:nz-1,2:nr-1);

    % Add missing pieces of FD equations at edges to source term
    source(1,:) = psibot*alpha/2/pi;
    source(end,:) = psitop*alpha/2/pi;
    source(:,1) = source(:,1)+psiinn(2:nz-1)*beta*(1+dr/2/rg(2))/2/pi;
    source(:,end) = source(:,end)+psiout(2:nz-1)*beta*(1-dr/2/rg(nr-1))/2/pi;

    % Inverse discrete sine transform in z direction of source
    isource = idst(source);

    % Use sparse matrices for the FD equations
    ammd = spparms('autommd');
    spparms('autommd',0);

    % Create system of equations
    isz = 2*(alpha+beta)-2*alpha*cos(pi*(1:nz-2)'/(nz-1)); % sine-transform solution for z
    % tridiagonal equations for r derivatives
    trid = sparse(1:nr-3,2:nr-2,-beta*(ones(1,nr-3)-dr/2./rg(2:nr-2)),nr-2,nr-2)+...
           sparse(2:nr-2,1:nr-3,-beta*(ones(1,nr-3)+dr/2./rg(3:nr-1)),nr-2,nr-2);

    % Solve all rows for inverse discrete sine transform
    for j = 1:nz-2
      trid(1:nr-1:(nr-2)^2) = isz(j)*ones(1,nr-2);
      ipsi(j,:) = (trid\isource(j,:)')';
    end
    spparms('autommd',ammd);

    % Sine-transform solution (ipsi) and surround by boundary values
    psizr_pla = [psiinn [psibot; 2*pi*dst(ipsi); psitop] psiout];

    psizr_plas(:,iflux) = psizr_pla(:);
    if nfluxes == 1
      psizr_plas = reshape(psizr_plas,n1,n2);
    end
    
    if nargout > 1 % return br
      br = [zeros(1,nr); (psizr_pla(1:end-2,:)-psizr_pla(3:end,:))/4/pi/dz./rgg(2:end-1,:); zeros(1,nr)];
      br(1,:) = 5*br(2,:)-7*br(3,:)+3*br(4,:); % Taylor expansion to order 2 around br(3,:)
      br(nz,:) = 5*br(nz-1,:)-7*br(nz-2,:)+3*br(nz-3,:);
      brs(:,iflux) = br(:);
      if nfluxes == 1
	brs = reshape(brs,n1,n2);
      end
    end

    if nargout > 2 % return bz
      bz = [zeros(nz,1) (psizr_pla(:,3:end)-psizr_pla(:,1:end-2))/4/pi/dr./rgg(:,2:end-1) zeros(nz,1)];
      bz(:,1) = 5*bz(:,2)-7*bz(:,3)+3*bz(:,4); % Taylor expansion to order 2 around bz(:,3)
      bz(:,nr) = 5*bz(:,nr-1)-7*bz(:,nr-2)+3*bz(:,nr-3);
      bzs(:,iflux) = bz(:);
      if nfluxes == 1
	bzs = reshape(bzs,n1,n2);
      end
    end
    
  end
    
  function b=dst(a)
  [n,m] = size(a);
  y=zeros(2*(n+1),m);
  y(2:n+1,:)=a;
  y(n+3:2*(n+1),:)=-flipud(a);
  yy=fft(y);
  b=yy(2:n+1,:)/(-2*sqrt(-1));
  if isreal(a), b = real(b); end

  function b=idst(a)
  n=size(a,1);
  nn=n+1;
  b=2/nn*dst(a);

