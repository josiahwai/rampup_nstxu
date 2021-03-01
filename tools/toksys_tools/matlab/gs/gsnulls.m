function nulls = gsnulls(y,nmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   nulls = gsnulls(y,nmax)
%
%  PURPOSE: Find nulls in matrix y (points where gradient = 0)
%
%  INPUTS:  y, a matrix
%           nmax, max number of returned nulls (default 9)
%
%  OUTPUTS:  nulls, structure with information about nulls
%            nulls.info describes all the fields
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Anders Welander ON 2016-11-03
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nz, nr] = size(y);

if nargin < 2
  nmax = 9; % Maximum number of archived nulls
end

% For cubic Hermite interpolation
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;
ir4 = (-1:2)'*nz; % For horizontal grid lines
iz4 = (-1:2)';    % For vertical grid lines
i16 = reshape(iz4*ones(1,4)+ones(4,1)*ir4',1,16);

% Allocate memory for the output
nulls.count = 0;
nulls.type = char(zeros(1,nmax)+' '); % OX*
nulls.r = nan(1,nmax);
nulls.z = nan(1,nmax);
nulls.y = nan(1,nmax);
nulls.i = nan(1,nmax);
nulls.j = nan(1,nmax);
nulls.k = nan(1,nmax);
nulls.yrr = nan(1,nmax);
nulls.yrz = nan(1,nmax);
nulls.yzz = nan(1,nmax);
nulls.yrrr = nan(1,nmax);
nulls.yrrz = nan(1,nmax);
nulls.yrzz = nan(1,nmax);
nulls.yzzz = nan(1,nmax);
nulls.ii = ones(nmax,16);
nulls.w = nan(nmax,16);
nulls.wrr = nan(nmax,16);
nulls.wrz = nan(nmax,16);
nulls.wzz = nan(nmax,16);
nulls.wrrr = nan(nmax,16);
nulls.wrrz = nan(nmax,16);
nulls.wrzz = nan(nmax,16);
nulls.wzzz = nan(nmax,16);
nulls.ur = nan(1,nmax);
nulls.uz = nan(1,nmax);
nulls.yuu = nan(1,nmax);
nulls.yvv = nan(1,nmax);
nulls.drdy = nan(nmax,16);
nulls.dzdy = nan(nmax,16);
nulls.dudy = nan(nmax,16);
nulls.durdy = nan(nmax,16);
nulls.duzdy = nan(nmax,16);
nulls.f = false(nz,nr);
ds.count = 'number of nulls';
ds.type = 'axis, x, or snowflake: OX*';
ds.r = 'floating index position of nulls';
ds.z = 'floating index position of nulls';
ds.y = 'value of y at r,z';
ds.i = 'floor(z)';
ds.j = 'floor(r)';
ds.k = 'y(k) = y(i,j)';
ds.yrr = 'value of d2y/dr2';
ds.yrz = 'value of d2y/drdz';
ds.yzz = 'value of d2y/dz2';
ds.yrrr = 'value of d3y/dr3';
ds.yrrz = 'value of d3y/dr2dz';
ds.yrzz = 'value of d3y/dr1dz2';
ds.yzzz = 'value of d3y/dz3';
ds.ii = 'Indices: nulls.y = sum(nulls.w.*y(nulls.ii),2)';
ds.w = 'weights to calculate y(r,z)';
ds.wrr = 'weights to calculate yrr(r,z)';
ds.wrz = 'weights to calculate yrz(r,z)';
ds.wzz = 'weights to calculate yzz(r,z)';
ds.wrrr = 'weights to calculate yrrr(r,z)';
ds.wrrz = 'weights to calculate yrrz(r,z)';
ds.wrzz = 'weights to calculate yrzz(r,z)';
ds.wzzz = 'weights to calculate yzzz(r,z)';
ds.ur = 'r of unit vector that maximizes d2y/du2';
ds.uz = 'z of unit vector that maximizes d2y/du2';
ds.yuu = 'max second derivative, d2y/du2';
ds.yvv = 'min second derivative, u*v''=0';
ds.drdy = 'r response to 4x4 values around it';
ds.dzdy = 'z response to 4x4 values around it';
ds.dudy = 'response of angle u (ur=cos(u),uz=sin(u))';
ds.durdy = 'response of ur to 4x4 values around it';
ds.duzdy = 'response of uz to 4x4 values around it';
ds.f = '1 means null may exist in range i:i+1,j:j+1';
nulls.info = ds;

% Variables to be used for finding grid cells that may contain nulls
tr = nan(nz,nr);
tz = nan(nz,nr);
ma = nan(nz,nr);
A = nan(nz,nr);
B = nan(nz,nr);
C = nan(nz,nr);
a = nan(nz,nr);
b = nan(nz,nr);
c = nan(nz,nr);
d0 = nan(nz,nr);
d1 = nan(nz,nr);
d2 = nan(nz,nr);
d3 = nan(nz,nr);
ra = zeros(nz,nr); % points where dydz = ha
rb = ones(nz,nr);  % points where dydz = hb
ha = nan(nz,nr);
hb = nan(nz,nr);
hc = nan(nz,nr);
h0 = nan(nz,nr);
h1 = nan(nz,nr);
h2 = nan(nz,nr);
r1 = nan(nz,nr);
r2 = nan(nz,nr);
r3 = zeros(nz,nr);
r4 = zeros(nz,nr);
r5 = ones(nz,nr);
za = zeros(nz,nr); % points where dydr = va
zb = ones(nz,nr);  % points where dydr = vb
va = nan(nz,nr);
vb = nan(nz,nr);
vc = nan(nz,nr);
v0 = nan(nz,nr);
v1 = nan(nz,nr);
v2 = nan(nz,nr);
z1 = nan(nz,nr);
z2 = nan(nz,nr);
z3 = zeros(nz,nr);
z4 = zeros(nz,nr);
z5 = ones(nz,nr);
da = nan(nz,nr);
db = nan(nz,nr);
dc = nan(nz,nr);
dd = nan(nz,nr);
nbr0 = nan(nz,nr); % intersections of Br=0 into rgg, zgg, rgg+dr, zgg+dz
nbz0 = nan(nz,nr); % intersections of Bz=0 into rgg, zgg, rgg+dr, zgg+dz


% HORIZONTAL GRID LINES
% y(r+dr) = A.*dr.^3+B.*dr.^2+C.*dr+y(r), dr={0,1}
% Cubic Hermite spline to get A, B, C
A(:,2:nr-2) = ( -y(:,1:nr-3) + 3*y(:,2:nr-2) - 3*y(:,3:nr-1)+y(:,4:nr))/2;
B(:,2:nr-2) = (2*y(:,1:nr-3) - 5*y(:,2:nr-2) + 4*y(:,3:nr-1)-y(:,4:nr))/2;
C(:,2:nr-2) = ( -y(:,1:nr-3)                 +   y(:,3:nr-1)          )/2;
% Solve for dy/dr = 0
d = -B./A/3;
s = d.^2-C./A/3;
o = s >= 0;
s(o) = sqrt(s(o));
r1(o) = d(o)-s(o); % Points with dydr = 0 along
r2(o) = d(o)+s(o); % horizontal grid lines
% 2*dydz(r+dr) = da.*dr.^3+db.*dr.^2+dc.*dr+dydz(r), dr={0,1}
da(2:nz-1,:) = A(3:nz,:) - A(1:nz-2,:);
db(2:nz-1,:) = B(3:nz,:) - B(1:nz-2,:);
dc(2:nz-1,:) = C(3:nz,:) - C(1:nz-2,:);
h0(2:nz-1,:) = y(3:nz,:) - y(1:nz-2,:);
% Solve for d2ydrdz(r) = 0 as first step in finding dydz(r) = 0
d = -db./da/3;
s = d.^2-dc./da/3;
o = s >= 0;
s(o) = sqrt(s(o));
ra(o) = d(o)-s(o); % Points with d2ydzdr = 0 along
rb(o) = d(o)+s(o); % horizontal grid lines
% dydz is monotonic between points 0, ra, rb, 1
ra(ra < 0 | ra > 1) = 0;
rb(rb < 0 | rb > 1) = 1;
ha = da.*ra.^3+db.*ra.^2+dc.*ra+h0;
hb = da.*rb.^3+db.*rb.^2+dc.*rb+h0;
hc(:,1:nr-1) = h0(:,2:nr);
o = h0.*ha < 0; % true when dydz==0 between 0 and ra
r3(o) = - h0(o).*ra(o)./(ha(o)-h0(o)); % Approximately where dydz=0
o = ha.*hb < 0; % true when dydz==0 between ra and rb
r4(o) = ra(o) - ha(o).*(rb(o)-ra(o))./(hb(o)-ha(o));
o = hb.*hc < 0; % true when dydz==0 between rb and 1
r5(o) = rb(o) - hb(o).*(1-rb(o))./(hc(o)-rb(o));

% VERTICAL GRID LINES
% y(z) = A.*z.^3+B.*z.^2+C.*z+y(0), z={0,1}
A(2:nz-2,:) = ( -y(1:nz-3,:) + 3*y(2:nz-2,:) - 3*y(3:nz-1,:)+y(4:nz,:))/2;
B(2:nz-2,:) = (2*y(1:nz-3,:) - 5*y(2:nz-2,:) + 4*y(3:nz-1,:)-y(4:nz,:))/2;
C(2:nz-2,:) = ( -y(1:nz-3,:)                 +   y(3:nz-1,:)          )/2;
d = -B./A/3;
s = d.^2-C./A/3;
o = s >= 0;
s(o) = sqrt(s(o));
z1(o) = d(o)-s(o); % Points with dydz = 0 along
z2(o) = d(o)+s(o); % vertical grid lines
% 2*dydr(z) = da.*z.^3+db.*z.^2+dc.*z+d0, z={0,1}
da(:,2:nr-1) = A(:,3:nr) - A(:,1:nr-2);
db(:,2:nr-1) = B(:,3:nr) - B(:,1:nr-2);
dc(:,2:nr-1) = C(:,3:nr) - C(:,1:nr-2);
v0(:,2:nr-1) = y(:,3:nr) - y(:,1:nr-2);
% Solve for d2ydrdz(z) = 0 as first step in finding dydz(z) = 0
d = -db./da/3;
s = d.^2-dc./da/3;
o = s >= 0;
s(o) = sqrt(s(o));
za(o) = d(o)-s(o); % Points with d2ydzdr = 0 along
zb(o) = d(o)+s(o); % vertical grid lines
za(za < 0 | za > 1) = 0;
zb(zb < 0 | zb > 1) = 1;
va = da.*za.^3+db.*za.^2+dc.*za+v0;
vb = da.*zb.^3+db.*zb.^2+dc.*zb+v0;
vc(1:nz-1,:) = v0(2:nz,:);
o = v0.*va < 0;
z3(o) = - v0(o).*za(o)./(va(o)-v0(o));
o = va.*vb < 0;
z4(o) = za(o) - va(o).*(zb(o)-za(o))./(vb(o)-va(o));
o = vb.*vc < 0;
z5(o) = zb(o) - vb(o).*(1-zb(o))./(vc(o)-vb(o));

nbr0(2:nz-2,2:nr-2) = ...
 (z1(2:nz-2,2:nr-2) > 0 & z1(2:nz-2,2:nr-2) < 1) + ...
 (z2(2:nz-2,2:nr-2) > 0 & z2(2:nz-2,2:nr-2) < 1) + ...
 (r3(2:nz-2,2:nr-2) > 0 & r3(2:nz-2,2:nr-2) < 1) + ...
 (r4(2:nz-2,2:nr-2) > 0 & r4(2:nz-2,2:nr-2) < 1) + ...
 (r5(2:nz-2,2:nr-2) > 0 & r5(2:nz-2,2:nr-2) < 1) + ...
 (z1(2:nz-2,3:nr-1) > 0 & z1(2:nz-2,3:nr-1) < 1) + ...
 (z2(2:nz-2,3:nr-1) > 0 & z2(2:nz-2,3:nr-1) < 1) + ...
 (r3(3:nz-1,2:nr-2) > 0 & r3(3:nz-1,2:nr-2) < 1) + ...
 (r4(3:nz-1,2:nr-2) > 0 & r4(3:nz-1,2:nr-2) < 1) + ...
 (r5(3:nz-1,2:nr-2) > 0 & r5(3:nz-1,2:nr-2) < 1);
nbz0(2:nz-2,2:nr-2) = ...
 (r1(2:nz-2,2:nr-2) > 0 & r1(2:nz-2,2:nr-2) < 1) + ...
 (r2(2:nz-2,2:nr-2) > 0 & r2(2:nz-2,2:nr-2) < 1) + ...
 (z3(2:nz-2,2:nr-2) > 0 & z3(2:nz-2,2:nr-2) < 1) + ...
 (z4(2:nz-2,2:nr-2) > 0 & z4(2:nz-2,2:nr-2) < 1) + ...
 (z5(2:nz-2,2:nr-2) > 0 & z5(2:nz-2,2:nr-2) < 1) + ...
 (r1(3:nz-1,2:nr-2) > 0 & r1(3:nz-1,2:nr-2) < 1) + ...
 (r2(3:nz-1,2:nr-2) > 0 & r2(3:nz-1,2:nr-2) < 1) + ...
 (z3(2:nz-2,3:nr-1) > 0 & z3(2:nz-2,3:nr-1) < 1) + ...
 (z4(2:nz-2,3:nr-1) > 0 & z4(2:nz-2,3:nr-1) < 1) + ...
 (z5(2:nz-2,3:nr-1) > 0 & z5(2:nz-2,3:nr-1) < 1);

nulls.f = nbz0 >= 2 & nbr0 >= 2;

% Each grid cell with nulls.f = true will be divided into nx*nx sub cells.
% This increases the chance of converging on a null in challenging cases
% that may have several nulls in the same cell. Failure can still occur.
nx = 11;
l = ones(nx*nx,1);
O = zeros(nx*nx,1);
UR = reshape(ones(nx,1)*linspace(0,1,nx),nx*nx,1);
UZ = reshape(linspace(0,1,nx)'*ones(1,nx),nx*nx,1);
Wr0 = [l UR UR.^2 UR.^3]*mx;
Wz0 = [l UZ UZ.^2 UZ.^3]*mx;
Wr1 = [O l 2*UR 3*UR.^2]*mx;
Wz1 = [O l 2*UZ 3*UZ.^2]*mx;
Wr = ...
[Wz0(:,1).*Wr1(:,1) Wz0(:,2).*Wr1(:,1) Wz0(:,3).*Wr1(:,1) Wz0(:,4).*Wr1(:,1) ...
 Wz0(:,1).*Wr1(:,2) Wz0(:,2).*Wr1(:,2) Wz0(:,3).*Wr1(:,2) Wz0(:,4).*Wr1(:,2) ...
 Wz0(:,1).*Wr1(:,3) Wz0(:,2).*Wr1(:,3) Wz0(:,3).*Wr1(:,3) Wz0(:,4).*Wr1(:,3) ...
 Wz0(:,1).*Wr1(:,4) Wz0(:,2).*Wr1(:,4) Wz0(:,3).*Wr1(:,4) Wz0(:,4).*Wr1(:,4)];
Wz = ...
[Wz1(:,1).*Wr0(:,1) Wz1(:,2).*Wr0(:,1) Wz1(:,3).*Wr0(:,1) Wz1(:,4).*Wr0(:,1) ...
 Wz1(:,1).*Wr0(:,2) Wz1(:,2).*Wr0(:,2) Wz1(:,3).*Wr0(:,2) Wz1(:,4).*Wr0(:,2) ...
 Wz1(:,1).*Wr0(:,3) Wz1(:,2).*Wr0(:,3) Wz1(:,3).*Wr0(:,3) Wz1(:,4).*Wr0(:,3) ...
 Wz1(:,1).*Wr0(:,4) Wz1(:,2).*Wr0(:,4) Wz1(:,3).*Wr0(:,4) Wz1(:,4).*Wr0(:,4)];

% Find all cells with nulls.f = true and pinpoint the nulls in them
for i = 1:nz
  for j = 1:nr
    if nulls.f(i,j) % Possibly nulls here
      Yr = reshape(Wr*y(i+(j-1)*nz+i16'),nx,nx);
      Yz = reshape(Wz*y(i+(j-1)*nz+i16'),nx,nx);
      mrx = nan(nx-1);
      mzx = nan(nx-1);
      for ii = 1:nx-1
        for jj = 1:nx-1
	  if (Yr(ii,jj)*Yr(ii+1,jj) <= 0 | Yr(ii+1,jj)*Yr(ii+1,jj+1) <= 0 | ...
	     Yr(ii+1,jj+1)*Yr(ii,jj+1) <= 0 | Yr(ii,jj+1)*Yr(ii,jj) <= 0) & ...
	     (Yz(ii,jj)*Yz(ii+1,jj) <= 0 | Yz(ii+1,jj)*Yz(ii+1,jj+1) <= 0 | ...
	     Yz(ii+1,jj+1)*Yz(ii,jj+1) <= 0 | Yz(ii,jj+1)*Yz(ii,jj) <= 0)
	    % This sub cell can contain a null
	    tr = (jj-0.5)/(nx-1);
	    tz = (ii-0.5)/(nx-1);
	    for k = 1:9
	      wz0 = [1 tz tz^2 tz^3]*mx;
	      wr0 = [1 tr tr^2 tr^3]*mx;
	      wz1 = [0 1 2*tz 3*tz^2]*mx;
	      wr1 = [0 1 2*tr 3*tr^2]*mx;
	      wz2 = [0 0 2 6*tz]*mx;
	      wr2 = [0 0 2 6*tr]*mx;
	      yr = wz0*y(i-1:i+2,j-1:j+2)*wr1';
	      yz = wz1*y(i-1:i+2,j-1:j+2)*wr0';
	      yrr = wz0*y(i-1:i+2,j-1:j+2)*wr2';
	      yrz = wz1*y(i-1:i+2,j-1:j+2)*wr1';
	      yzz = wz2*y(i-1:i+2,j-1:j+2)*wr0';
	      cnull = -[yrr yrz; yrz yzz]\[yr; yz];
	      ma = max(abs(cnull));
	      cnull = cnull/max(1,8*ma);
	      tr = tr+cnull(1);
	      tz = tz+cnull(2);
	      if ma < 1e-14
		break % converged on the null
	      end
	    end % End of k loop
	    if ma < 1e-4
	      iii = max(1,min(nx-1,1+floor((nx-1)*tz)));
	      jjj = max(1,min(nx-1,1+floor((nx-1)*tr)));
	      % Null might be slightly outside range because of rounding errors
	      if tz > -1e-4 & tz < 1+1e-4 & tr > -1e-4 & tr < 1+1e-4
		mrx(iii,jjj) = tr;
		mzx(iii,jjj) = tz;
	      end
	    end
	  end % End of if possible null in sub cell
        end % End of jj loop
      end % End of ii loop
      for ii = 1:nx-1
	for jj = 1:nx-1
	  if ~isnan(mrx(ii,jj)) & nulls.count < nmax
            nulls.count = nulls.count + 1;
	    nulls.i(nulls.count) = i;
            nulls.j(nulls.count) = j;
            nulls.z(nulls.count) = i+mzx(ii,jj);
            nulls.r(nulls.count) = j+mrx(ii,jj);
	  end
	end
      end
    end
  end
end

% Remove duplicates (points within a distance 1e-4 of each other)
count = nulls.count;
for n = 1:nulls.count-1
  for m = n+1:nulls.count
    if (nulls.r(n)-nulls.r(m))^2+(nulls.z(n)-nulls.z(m))^2 < 1e-8
      if nulls.r(n)-nulls.j(n) < 0 | nulls.r(n)-nulls.j(n) > 1 | ...
	 nulls.z(n)-nulls.i(n) < 0 | nulls.z(n)-nulls.i(n) > 1
	k = n; % point n is out of range
      else
	k = m; % maybe m is out of range, either way one has to go
      end
      % Mark point k for deletion
      count = count-1;
      nulls.i(k) = nan;
      nulls.j(k) = nan;
      nulls.z(k) = nan;
      nulls.r(k) = nan;
    end
  end
end
nulls.count = count;
[~,kn] = sort(nulls.i+nulls.j*nz); % moves deleted points to the end
nulls.i = nulls.i(kn);
nulls.j = nulls.j(kn);
nulls.z = nulls.z(kn);
nulls.r = nulls.r(kn);

v = [-3:3]/8*pi;
vr = cos(v);
vz = sin(v);

% Fill in information about the nulls
for n = 1:nulls.count

  % Indices
  i = nulls.i(n);
  j = nulls.j(n);
  
  % Relative r, z
  tr = nulls.r(n)-j;
  tz = nulls.z(n)-i;
  
  % z,r weights
  wz0 = [1 tz tz^2 tz^3]*mx;
  wr0 = [1 tr tr^2 tr^3]*mx;
  wz1 = [0 1 2*tz 3*tz^2]*mx;
  wr1 = [0 1 2*tr 3*tr^2]*mx;
  wz2 = [0 0 2 6*tz]*mx;
  wr2 = [0 0 2 6*tr]*mx;
  wz3 = [0 0 0 6]*mx;
  wr3 = [0 0 0 6]*mx;
  
  % Derivatives
  yr = wz0*y(i-1:i+2,j-1:j+2)*wr1';
  yz = wz1*y(i-1:i+2,j-1:j+2)*wr0';
  yrr = wz0*y(i-1:i+2,j-1:j+2)*wr2';
  yrz = wz1*y(i-1:i+2,j-1:j+2)*wr1';
  yzz = wz2*y(i-1:i+2,j-1:j+2)*wr0';
  yrrr = wz0*y(i-1:i+2,j-1:j+2)*wr3';
  yrrz = wz1*y(i-1:i+2,j-1:j+2)*wr2';
  yrzz = wz2*y(i-1:i+2,j-1:j+2)*wr1';
  yzzz = wz3*y(i-1:i+2,j-1:j+2)*wr0';
  
  % Weights
  w =  reshape(wz0'*wr0,1,16);
  wr =  reshape(wz0'*wr1,1,16);
  wz =  reshape(wz1'*wr0,1,16);
  wrr =  reshape(wz0'*wr2,1,16);
  wrz =  reshape(wz1'*wr1,1,16);
  wzz =  reshape(wz2'*wr0,1,16);
  wrrr =  reshape(wz0'*wr3,1,16);
  wrrz =  reshape(wz1'*wr2,1,16);
  wrzz =  reshape(wz2'*wr1,1,16);
  wzzz =  reshape(wz3'*wr0,1,16);
  
  % Principal axis
  Yvv = yrr*vz.^2-2*yrz*vr.*vz+yzz*vr.^2;
  for k = 2:6
    if (Yvv(k)-Yvv(k-1))*(Yvv(k)-Yvv(k+1)) >= 0
      break
    end
  end
  u = (k-4)/8*pi;
  for k = 1:3
    ur = cos(u);
    uz = sin(u);
    D3u = (yzz-yrr)*ur*uz+yrz*(ur^2-uz^2); % d(d2)/dv should be 0
    D4u = (yzz-yrr)*(ur^2-uz^2)-4*yrz*ur*uz;
    u = u-D3u/D4u;
  end
  ur = cos(u);
  uz = sin(u);
  yuu = yrr*ur^2+2*yrz*ur*uz+yzz*uz^2;
  yvv = yrr*uz^2-2*yrz*ur*uz+yzz*ur^2;
  
  % Response of x-point position
  drzdgp = -inv([yrr yrz; yrz yzz]);
  drdy = drzdgp(1,1)*wr+drzdgp(1,2)*wz;
  dzdy = drzdgp(2,1)*wr+drzdgp(2,2)*wz;
  
  % Response of principal axis (dD3udy = 0, solve for dudy)
  yrrw = wrr + yrrr*drdy + yrrz*dzdy;
  yrzw = wrz + yrrz*drdy + yrzz*dzdy;
  yzzw = wzz + yrzz*drdy + yzzz*dzdy;
  den = (yzz-yrr)*(ur^2-uz^2) - yrz*ur*uz*4;	  
  dudy = -((yzzw-yrrw)*ur*uz + yrzw*(ur^2-uz^2))/den;
  durdy = -uz*dudy;
  duzdy = +ur*dudy;

  % Archive the results
  if yuu*yvv > 0
    nulls.type(n) = 'O';
  else
    nulls.type(n) = 'X';
  end
  nulls.k(n) = i+nz*(j-1);
  nulls.ii(n,:) = i+nz*(j-1)+i16;
  nulls.y(n) = wz0*y(i-1:i+2,j-1:j+2)*wr0';
  nulls.yrr(n) = yrr;
  nulls.yrz(n) = yrz;
  nulls.yzz(n) = yzz;
  nulls.yrrr(n) = yrrr;
  nulls.yrrz(n) = yrrz;
  nulls.yrzz(n) = yrzz;
  nulls.yzzz(n) = yzzz;
  nulls.w(n,:) = w;
  nulls.wrr(n,:) = wrr;
  nulls.wrz(n,:) = wrz;
  nulls.wzz(n,:) = wzz;
  nulls.wrrr(n,:) = wrrr;
  nulls.wrrz(n,:) = wrrz;
  nulls.wrzz(n,:) = wrzz;
  nulls.wzzz(n,:) = wzzz;
  nulls.ur(n) = ur;
  nulls.uz(n) = uz;
  nulls.yuu(n) = yuu;
  nulls.yvv(n) = yvv;
  nulls.drdy(n,:) = drdy;
  nulls.dzdy(n,:) = dzdy;
  nulls.dudy(n,:) = dudy;
  nulls.durdy(n,:) = durdy;
  nulls.duzdy(n,:) = duzdy;
  
end

