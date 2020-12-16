%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_nulls_xlim
%
%  PURPOSE: Find all nulls of interest for gsdesign of x point limits
%
%  INPUTS:  psizr, flux on grid
%
%  OUTPUTS:  nulls_xlim, a structure with information about nulls
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Anders Welander ON 2015-10-13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% One time assignments:
azr = zeros(nz,nr);
bzr = zeros(nz,nr);
czr = zeros(nz,nr);
dzr = zeros(nz,nr);
xzr0 = zeros(nz,nr);
xzr1 = zeros(nz,nr);
xzr2 = zeros(nz,nr);
xzr3 = zeros(nz,nr);
xzr4 = zeros(nz,nr);
xzr5 = zeros(nz,nr);
xzr6 = zeros(nz,nr);
xzr7 = zeros(nz,nr);
xzr8 = zeros(nz,nr);
xzr9 = zeros(nz,nr);
dazr = zeros(nz,nr);
dbzr = zeros(nz,nr);
dczr = zeros(nz,nr);
ddzr = zeros(nz,nr);
nbr0 = zeros(nz,nr); % intersections of Br=0 into rgg, zgg, rgg+dr, zgg+dz
nbz0 = zeros(nz,nr); % intersections of Bz=0 into rgg, zgg, rgg+dr, zgg+dz

% VERTICAL GRID LINES

azr(2:nz-2,:) = ( -psizr(1:nz-3,:) + 3*psizr(2:nz-2,:) - 3*psizr(3:nz-1,:)+psizr(4:nz,:))/2;
bzr(2:nz-2,:) = (2*psizr(1:nz-3,:) - 5*psizr(2:nz-2,:) + 4*psizr(3:nz-1,:)-psizr(4:nz,:))/2;
czr(2:nz-2,:) = ( -psizr(1:nz-3,:)                     +   psizr(3:nz-1,:)              )/2;
dzr(2:nz-2,:) =                        psizr(2:nz-2,:) -   psibry;

% Points with dpsizrdz = 0 along vertical grid lines:
xzr0 = (-bzr+sqrt(bzr.^2-3*azr.*czr))./(3*azr);
xzr1 = (-bzr-sqrt(bzr.^2-3*azr.*czr))./(3*azr);

% Find where dpsizrdr = 0 along vertical grid lines
dazr(:,2:nr-1) = azr(:,3:nr) - azr(:,1:nr-2);
dbzr(:,2:nr-1) = bzr(:,3:nr) - bzr(:,1:nr-2);
dczr(:,2:nr-1) = czr(:,3:nr) - czr(:,1:nr-2);
ddzr(:,2:nr-1) = dzr(:,3:nr) - dzr(:,1:nr-2);
[xzr2 xzr3 xzr4] = cubicroots(dazr, dbzr, dczr, ddzr);

% HORIZONTAL GRID LINES

azr(:,2:nr-2) = ( -psizr(:,1:nr-3) + 3*psizr(:,2:nr-2) - 3*psizr(:,3:nr-1)+psizr(:,4:nr))/2;
bzr(:,2:nr-2) = (2*psizr(:,1:nr-3) - 5*psizr(:,2:nr-2) + 4*psizr(:,3:nr-1)-psizr(:,4:nr))/2;
czr(:,2:nr-2) = ( -psizr(:,1:nr-3)                     +   psizr(:,3:nr-1)              )/2;
dzr(:,2:nr-2) =                        psizr(:,2:nr-2) -   psibry;

% Points with dpsizrdr = 0 along horizontal grid lines:
xzr5 = (-bzr+sqrt(bzr.^2-3*azr.*czr))./(3*azr);
xzr6 = (-bzr-sqrt(bzr.^2-3*azr.*czr))./(3*azr);

% Find where dpsizrdz = 0 along horizontal grid lines
dazr(2:nz-1,:) = azr(3:nz,:) - azr(1:nz-2,:);
dbzr(2:nz-1,:) = bzr(3:nz,:) - bzr(1:nz-2,:);
dczr(2:nz-1,:) = czr(3:nz,:) - czr(1:nz-2,:);
ddzr(2:nz-1,:) = dzr(3:nz,:) - dzr(1:nz-2,:);
[xzr7 xzr8 xzr9] = cubicroots(dazr, dbzr, dczr, ddzr);

xzr0(~isreal(xzr0)) = nan;
xzr1(~isreal(xzr1)) = nan;
xzr2(~isreal(xzr2)) = nan;
xzr3(~isreal(xzr3)) = nan;
xzr4(~isreal(xzr4)) = nan;
xzr5(~isreal(xzr5)) = nan;
xzr6(~isreal(xzr6)) = nan;
xzr7(~isreal(xzr7)) = nan;
xzr8(~isreal(xzr8)) = nan;
xzr9(~isreal(xzr9)) = nan;


nbr0(1:nz-1,1:nr-1) = ...
  (xzr0(1:nz-1,1:nr-1) >= 0 & xzr0(1:nz-1,1:nr-1) < 1) + ...
  (xzr1(1:nz-1,1:nr-1) >= 0 & xzr1(1:nz-1,1:nr-1) < 1) + ...
  (xzr7(1:nz-1,1:nr-1) >= 0 & xzr7(1:nz-1,1:nr-1) < 1) + ...
  (xzr8(1:nz-1,1:nr-1) >= 0 & xzr8(1:nz-1,1:nr-1) < 1) + ...
  (xzr9(1:nz-1,1:nr-1) >= 0 & xzr9(1:nz-1,1:nr-1) < 1) + ...
  (xzr0(1:nz-1,2:nr  ) >= 0 & xzr0(1:nz-1,2:nr  ) < 1) + ...
  (xzr1(1:nz-1,2:nr  ) >= 0 & xzr1(1:nz-1,2:nr  ) < 1) + ...
  (xzr7(2:nz  ,1:nr-1) >= 0 & xzr7(2:nz  ,1:nr-1) < 1) + ...
  (xzr8(2:nz  ,1:nr-1) >= 0 & xzr8(2:nz  ,1:nr-1) < 1) + ...
  (xzr9(2:nz  ,1:nr-1) >= 0 & xzr9(2:nz  ,1:nr-1) < 1);
nbz0(1:nz-1,1:nr-1) = ...
  (xzr5(1:nz-1,1:nr-1) >= 0 & xzr5(1:nz-1,1:nr-1) < 1) + ...
  (xzr6(1:nz-1,1:nr-1) >= 0 & xzr6(1:nz-1,1:nr-1) < 1) + ...
  (xzr2(1:nz-1,1:nr-1) >= 0 & xzr2(1:nz-1,1:nr-1) < 1) + ...
  (xzr3(1:nz-1,1:nr-1) >= 0 & xzr3(1:nz-1,1:nr-1) < 1) + ...
  (xzr4(1:nz-1,1:nr-1) >= 0 & xzr4(1:nz-1,1:nr-1) < 1) + ...
  (xzr5(2:nz  ,1:nr-1) >= 0 & xzr5(2:nz  ,1:nr-1) < 1) + ...
  (xzr6(2:nz  ,1:nr-1) >= 0 & xzr6(2:nz  ,1:nr-1) < 1) + ...
  (xzr2(1:nz-1,2:nr  ) >= 0 & xzr2(1:nz-1,2:nr  ) < 1) + ...
  (xzr3(1:nz-1,2:nr  ) >= 0 & xzr3(1:nz-1,2:nr  ) < 1) + ...
  (xzr4(1:nz-1,2:nr  ) >= 0 & xzr4(1:nz-1,2:nr  ) < 1);


nulls_xlim.count = 0;
for i = 2:nz-2
  for j = 2:nr-2
    if ((rg(j)-rmaxis)/dr)^2+((zg(i)-zmaxis)/dz)^2 > 4 & ...
       nbr0(i,j) > 0 & nbz0(i,j) > 0 % null inside this cell?
     
      pp = reshape(psizr(i-1:i+2,j-1:j+2),16,1);
      m = 0;

      if isreal(xzr5(i,j)) & xzr5(i,j) >= 0 & xzr5(i,j) < 1
	m = m+1;
	bzx(m) = xzr5(i,j);
	bzy(m) = 0;
      end
      if isreal(xzr5(i+1,j)) & xzr5(i+1,j) >= 0 & xzr5(i+1,j) < 1
	m = m+1;
	bzx(m) = xzr5(i+1,j);
	bzy(m) = 1;
      end
      if isreal(xzr6(i,j)) & xzr6(i,j) >= 0 & xzr6(i,j) < 1
	m = m+1;
	bzx(m) = xzr6(i,j);
	bzy(m) = 0;
      end
      if isreal(xzr6(i+1,j)) & xzr6(i+1,j) >= 0 & xzr6(i+1,j) < 1
	m = m+1;
	bzx(m) = xzr6(i+1,j);
	bzy(m) = 1;
      end
      if isreal(xzr2(i,j)) & xzr2(i,j) >= 0 & xzr2(i,j) < 1
	m = m+1;
	bzx(m) = 0;
	bzy(m) = xzr2(i,j);
      end
      if isreal(xzr2(i,j+1)) & xzr2(i,j+1) >= 0 & xzr2(i,j+1) < 1
	m = m+1;
	bzx(m) = 0+1;
	bzy(m) = xzr2(i,j+1);
      end
      if isreal(xzr3(i,j)) & xzr3(i,j) >= 0 & xzr3(i,j) < 1
	m = m+1;
	bzx(m) = 0;
	bzy(m) = xzr3(i,j);
      end
      if isreal(xzr3(i,j+1)) & xzr3(i,j+1) >= 0 & xzr3(i,j+1) < 1
	m = m+1;
	bzx(m) = 0+1;
	bzy(m) = xzr3(i,j+1);
      end
      if isreal(xzr4(i,j)) & xzr4(i,j) >= 0 & xzr4(i,j) < 1
	m = m+1;
	bzx(m) = 0;
	bzy(m) = xzr4(i,j);
      end
      if isreal(xzr4(i,j+1)) & xzr4(i,j+1) >= 0 & xzr4(i,j+1) < 1
	m = m+1;
	bzx(m) = 0+1;
	bzy(m) = xzr4(i,j+1);
      end

      for k = 1:m % Trace the contour bz=0 inside this cell from all starting points on sides
	[tr, tz, p, vx, vz] = ...
	  gs_trace_bznull_in_cell(bzx(k),bzy(k),0,1,0,1,psizr(i-1:i+2,j-1:j+2));
	if p == 0 % A null was found
	  rnull = rg(j)+tr*dr;
	  znull = zg(i)+tz*dz;
	  % Check if this null was already found
	  new_null = true;
	  for l = 1:nulls_xlim.count
	    if ((nulls_xlim.r(l)-rnull)/dr)^2+((nulls_xlim.z(l)-znull)/dz)^2 < 1e-2
	      new_null = false;
	    end
	  end
	  if new_null
	    w =   reshape(([1 tz   tz^2   tz^3]*mx)'*[1 tr   tr^2   tr^3]*mx,1,16);
	    wr =  reshape(([1 tz   tz^2   tz^3]*mx)'*[0 1  2*tr   3*tr^2]*mx,1,16);
	    wz =  reshape(([0 1  2*tz   3*tz^2]*mx)'*[1 tr   tr^2   tr^3]*mx,1,16);
	    wrr = reshape(([1 tz   tz^2   tz^3]*mx)'*[0 0  2      6*tr  ]*mx,1,16);
	    wrz = reshape(([0 1  2*tz   3*tz^2]*mx)'*[0 1  2*tr   3*tr^2]*mx,1,16);
	    wzz = reshape(([0 0    2    6*tz  ]*mx)'*[1 tr   tr^2   tr^3]*mx,1,16);
	    yr = wr*pp;
	    yz = wz*pp;
	    yrr = wrr*pp;
	    yrz = wrz*pp;
	    yzz = wzz*pp;
	    drzdgp = -inv([yrr yrz; yrz yzz]);
	    cnull = drzdgp*[yr; yz];
	    nulls_xlim.count = nulls_xlim.count + 1;
	    nulls_xlim.r(nulls_xlim.count) = rnull;
	    nulls_xlim.z(nulls_xlim.count) = znull;
	    nulls_xlim.psi(nulls_xlim.count) = w*pp;
            nulls_xlim.psi_r(nulls_xlim.count) = yr;
            nulls_xlim.psi_z(nulls_xlim.count) = yz;
            nulls_xlim.psi_rr(nulls_xlim.count) = yrr;
            nulls_xlim.psi_rz(nulls_xlim.count) = yrz;
            nulls_xlim.psi_zz(nulls_xlim.count) = yzz;
	    nulls_xlim.drdpsi(nulls_xlim.count,:) = (drzdgp(1,1)*wr+drzdgp(1,2)*wz)*dr;
            nulls_xlim.dzdpsi(nulls_xlim.count,:) = (drzdgp(2,1)*wr+drzdgp(2,2)*wz)*dz;
	    nulls_xlim.i(nulls_xlim.count) = i;
	    nulls_xlim.j(nulls_xlim.count) = j;
	    nulls_xlim.k(nulls_xlim.count) = i+nz*(j-1);
	    nulls_xlim.w(nulls_xlim.count,:) = w;
	    nulls_xlim.ii(nulls_xlim.count,:) = i+nz*(j-1)+neighbors;
            nulls_xlim.r_uncertainty(nulls_xlim.count) = cnull(1);
            nulls_xlim.z_uncertainty(nulls_xlim.count) = cnull(2);
	  end
	end
      end
    end
  end
end

% Sort the nulls so that the most relevant comes first
n = nulls_xlim.count;
if n > 0 % needed since the field psi may not exist if n == 0
  [~,kk(1:n)] = sort((nulls_xlim.psi(1:n)-psimag)/(psibry-psimag));
  nulls_xlim.r = nulls_xlim.r(kk(1:n));
  nulls_xlim.z = nulls_xlim.z(kk(1:n));
  nulls_xlim.psi = nulls_xlim.psi(kk(1:n));
  nulls_xlim.psi_r = nulls_xlim.psi_r(kk(1:n));
  nulls_xlim.psi_z = nulls_xlim.psi_z(kk(1:n));
  nulls_xlim.psi_rr = nulls_xlim.psi_rr(kk(1:n));
  nulls_xlim.psi_rz = nulls_xlim.psi_rz(kk(1:n));
  nulls_xlim.psi_zz = nulls_xlim.psi_zz(kk(1:n));
  nulls_xlim.drdpsi = nulls_xlim.drdpsi(kk(1:n),:);
  nulls_xlim.dzdpsi = nulls_xlim.dzdpsi(kk(1:n),:);
  nulls_xlim.i = nulls_xlim.i(kk(1:n));
  nulls_xlim.j = nulls_xlim.j(kk(1:n));
  nulls_xlim.k = nulls_xlim.k(kk(1:n));
  nulls_xlim.w = nulls_xlim.w(kk(1:n),:);
  nulls_xlim.ii = nulls_xlim.ii(kk(1:n),:);
  nulls_xlim.r_uncertainty = nulls_xlim.r_uncertainty(kk(1:n));
  nulls_xlim.z_uncertainty = nulls_xlim.z_uncertainty(kk(1:n));
elseif ~isfield(nulls,'r')
  nulls_xlim.r = [];
  nulls_xlim.z = [];
  nulls_xlim.psi = [];
  nulls_xlim.psi_r = [];
  nulls_xlim.psi_z = [];
  nulls_xlim.psi_rr = [];
  nulls_xlim.psi_rz = [];
  nulls_xlim.psi_zz = [];
  nulls_xlim.drdpsi = zeros(0,16);
  nulls_xlim.dzdpsi = zeros(0,16);
  nulls_xlim.i = [];
  nulls_xlim.j = [];
  nulls_xlim.k = [];
  nulls_xlim.w = zeros(0,16);
  nulls_xlim.ii = zeros(0,16);
  nulls_xlim.r_uncertainty = [];
  nulls_xlim.z_uncertainty = [];
end

return


clf
contour(rg,zg,psizr,99)
hold on
plot(rgg,zgg,'x')
for i = 1:nz
  for j = 1:nr
    if isreal(xzr0(i,j)) & xzr0(i,j) >= 0 & xzr0(i,j) < 1
      plot(rg(j), zg(i)+real(xzr0(i,j))*dz,'ro','linew',4,'markers',12)
    end
    if isreal(xzr1(i,j)) & xzr1(i,j) >= 0 & xzr1(i,j) < 1
      plot(rg(j), zg(i)+real(xzr1(i,j))*dz,'ro','linew',4,'markers',12)
    end
    if isreal(xzr7(i,j)) & xzr7(i,j) >= 0 & xzr7(i,j) < 1
      plot(rg(j)+real(xzr7(i,j))*dr, zg(i),'ro','linew',4,'markers',12)
    end
    if isreal(xzr8(i,j)) & xzr8(i,j) >= 0 & xzr8(i,j) < 1
      plot(rg(j)+real(xzr8(i,j))*dr, zg(i),'ro','linew',4,'markers',12)
    end
    if isreal(xzr9(i,j)) & xzr9(i,j) >= 0 & xzr9(i,j) < 1
      plot(rg(j)+real(xzr9(i,j))*dr, zg(i),'ro','linew',4,'markers',12)
    end
    
    if isreal(xzr5(i,j)) & xzr5(i,j) >= 0 & xzr5(i,j) < 1
      plot(rg(j)+real(xzr5(i,j))*dr, zg(i),'bo','linew',4,'markers',12)
    end
    if isreal(xzr6(i,j)) & xzr6(i,j) >= 0 & xzr6(i,j) < 1
      plot(rg(j)+real(xzr6(i,j))*dr, zg(i),'bo','linew',4,'markers',12)
    end
    if isreal(xzr2(i,j)) & xzr2(i,j) >= 0 & xzr2(i,j) < 1
      plot(rg(j), zg(i)+real(xzr2(i,j))*dz,'bo','linew',4,'markers',12)
    end
    if isreal(xzr3(i,j)) & xzr3(i,j) >= 0 & xzr3(i,j) < 1
      plot(rg(j), zg(i)+real(xzr3(i,j))*dz,'bo','linew',4,'markers',12)
    end
    if isreal(xzr4(i,j)) & xzr4(i,j) >= 0 & xzr4(i,j) < 1
      plot(rg(j), zg(i)+real(xzr4(i,j))*dz,'bo','linew',4,'markers',12)
    end
  end
end

plot(rl,zl,'k')
hold off
