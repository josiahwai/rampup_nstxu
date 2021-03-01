function [M, Mr, Mz] = mpolygon2point(rp, zp, r, z, D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  [M, Mr, Mz] = mpolygon2point(rpolygon, zpolygon, r, z, D)
%
%  PURPOSE:  Calculate mutual inductances between coaxial loops
%            First loop has cross section described by rpolygon, zpolygon
%              which contain coordinates for the corners
%            Cross section of second loop is a point (thin wire) at r, z
%
%  INPUTS: rpolygon, zpolygon = coordinates of first current loop [meters]
%          r, z = major radius, height of second loop [meters]
%          D, choice of derivatives in Mr, Mz, see below
%
%  OUTPUTS: M = mutual inductance between loops [H]
%           Mr, Mz, derivatives w.r.t. r, z
%           if D = 1, derivatives are w.r.t. rigid shifts of rp, zp
%           if D = 2, derivatives are w.r.t. r, z
%
%  RESTRICTIONS: The mutual M is for uniform current density across polygon

%  METHOD:  A grid is created around the polygon, fields inside
%   are calculated using the grid solution and outside using mutinds
%
%  WRITTEN BY: Anders Welander ON 2015-04-03
%
%  MODIFICATION: 2016-04-18, adding Mr, Mz as optional outputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('D','var')
  D = 1;
end

% Last polygon corner same as first
if rp(end) ~= rp(1) || zp(end) ~= zp(1)
  p.r = [rp(:)' rp(1)];
  p.z = [zp(:)' zp(1)];
else
  p.r = rp(:)';
  p.z = zp(:)';
end

% Remove any duplicates of polygon corners
ii = [1 diff(p.r)] | [1 diff(p.z)];
p.r = p.r(ii);
p.z = p.z(ii);

p.n = length(p.z);

% Unique z values
p.unique.z = unique(p.z);
p.unique.n = length(p.unique.z);

% Find ranges in r inside polygon at unique z
% nranges will hold maximum number of ranges
nranges = 1;
% rranges are valid just above unique z
p.rranges = nan(p.unique.n,2*nranges);
p.drdz = nan(p.unique.n,2*nranges);
for i = 1:p.unique.n
  zz = p.unique.z(i);
  for j = 1:p.n-1
    dz = p.z(j+1)-p.z(j);
    if p.z(j) <= zz && zz < p.z(j+1) || ...
       p.z(j) > zz && zz >= p.z(j+1)
      % interpolate to find rr at zz
      dr = p.r(j+1)-p.r(j);
      drdz = dr/dz;
      f = (zz-p.z(j));
      rr = p.r(j) + f*drdz;
      % Put rr & drdz in first free entry
      k = 1;
      while ~isnan(p.rranges(i,k))
        k = k+1;
	if k > 2*nranges
	  nranges = nranges + 1;
	  p.rranges(:,k:2*nranges) = nan;
	  p.drdz(:,k:2*nranges) = nan;
	end
      end
      p.rranges(i,k) = rr;
      p.drdz(i,k) = drdz;
    end % End of < = > tests
  end % End of j loop
  [p.rranges(i,:), kranges] = sort(p.rranges(i,:));
  p.drdz(i,:) = p.drdz(i,kranges);
  for j = 1:2:size(p.rranges,2)-1 % Sorting won't work for zero ranges
    if p.rranges(i,j) == p.rranges(i,j+1)
      p.drdz(i,j:j+1) = sort(p.drdz(i,j:j+1)); % Fix for zero range
    end
  end
end % End of i loop
p.nranges = nranges;

% Create a grid
nr = 65;
nz = 65;
ngg = nr*nz;
DR = 2*(max(rp)-min(rp));
DZ = 2*(max(zp)-min(zp));
dr = DR/(nr-1);
dz = DZ/(nz-1);
rgmin = max(dr,min(rp)-DR/4);
rgmax = max(rp)+DR/4;
DR = rgmax-rgmin;
dr = DR/(nr-1);
rg = linspace(rgmin,rgmax,nr);
zgmin = min(zp)-DZ/4;
zgmax = max(zp)+DZ/4;
zg = linspace(zgmin,zgmax,nz)';
rgg = ones(nz,1)*rg;
zgg = zg*ones(1,nr);
A = zeros(nz,nr);
RA = zeros(nz,nr);
ZA = zeros(nz,nr);

% Find z values where polygon crosses over to new grid rectangle
rpmin = min(rp);
rpmax = max(rp);
zpmin = min(zp);
zpmax = max(zp);
irp = (rp-rpmin)/(rpmax-rpmin)*(nr-1)/2;
izp = (zp-zpmin)/(zpmax-zpmin)*(nz-1)/2;
dirp = diff(irp);
dizp = diff(izp);
jr = round(irp);
iz = round(izp);
djr = diff(jr);
n = sum(abs(djr));
jzx = zeros(1,n);
izx = zeros(1,n);
j = 0;
for i = 1:length(rp)-1
  if djr(i) > 0
    jrx(j+(1:+djr(i))) = jr(i)+0.5:jr(i+1);
    izx(j+(1:+djr(i))) = izp(i)+((jr(i)+0.5:jr(i+1))-irp(i))/dirp(i)*dizp(i);
    j = j+djr(i);
  elseif djr(i) < 0
    jrx(j+(1:-djr(i))) = jr(i)-0.5:-1:jr(i+1);
    izx(j+(1:-djr(i))) = izp(i)+((jr(i)-0.5:-1:jr(i+1))-irp(i))/dirp(i)*dizp(i);
    j = j-djr(i);
  end
end
% Z values where polygon contour crosses vertical borders into new cells
zx = unique(izx*2/(nz-1)*(zpmax-zpmin)+zpmin);
zunique = unique([p.unique.z zx]);
n = 1; % Pointer to p.unique.z
for k = 1:length(zunique)-1
  while n < p.unique.n & p.unique.z(n+1) <= zunique(k)
    n = n+1;
  end
  for i = 1:nz
    if zg(i)+dz/2 > zunique(k) & zg(i)-dz/2 < zunique(k+1)
      zd = max(zunique(k),zg(i)-dz/2);
      zu = min(zunique(k+1),zg(i)+dz/2);
      d = zu-zd;
      for m = 1:p.nranges % Index to rranges
	rdi = p.rranges(n,m*2-1)+p.drdz(n,m*2-1)*(zd-p.unique.z(n));
	rui = p.rranges(n,m*2-1)+p.drdz(n,m*2-1)*(zu-p.unique.z(n));
	rdo = p.rranges(n,m*2-0)+p.drdz(n,m*2-0)*(zd-p.unique.z(n));
	ruo = p.rranges(n,m*2-0)+p.drdz(n,m*2-0)*(zu-p.unique.z(n));
	rmin = min(rdi,rui);
	rmax = max(rdo,ruo);
	for j = 1:nr
	  if rg(j)+dr/2 > rmin & rg(j)-dr/2 < rmax
            rd1 = max(rdi,rg(j)-dr/2);
            ru1 = max(rui,rg(j)-dr/2);
            rd2 = max(rd1,min(rdo,rg(j)+dr/2));
            ru2 = max(ru1,min(ruo,rg(j)+dr/2));
	    drdz1 = (ru1-rd1)/d;
	    drdz2 = (ru2-rd2)/d;
	    % integral(dr) = (rd2+drdz2*z)-(rd1+drdz1*z)
	    dA = (rd2+ru2-rd1-ru1)*d/2;
	    A(i,j) = A(i,j) + dA;
	    % integral(r*dr) = (rd2+drdz2*z)^2/2-(rd1+drdz1*z)^2/2
	    RA(i,j) = RA(i,j) + ...
	      (rd2^2 + rd2*drdz2*d + drdz2^2*d^2/3)*d/2 - ...
	      (rd1^2 + rd1*drdz1*d + drdz1^2*d^2/3)*d/2;
	    % integral[(rd2+drdz2*z)-(rd1+drdz1*z)]*(zd+z)
	    ZA(i,j) = ZA(i,j) + dA*zd + (rd2-rd1)*d^2/2 + (drdz2-drdz1)*d^3/3;
	  end
	end
      end
    end
  end
end

Atot = sum(A(:));
ia = A > 0;
rc = RA(ia)./A(ia);
zc = ZA(ia)./A(ia);
fc = A(ia)/Atot;

if numel(r) == 1 & numel(z) > 1
  r = r+0*z;
elseif numel(r) > 1 & numel(z) == 1
  z = z+0*r;
end
M = r; % Make M the right size
Mr = r;
Mz = r;
for i = 1:numel(M)
  if nargout > 1
    if D == 1
      [ms, mrs, mzs] = mutinds(r(i),z(i),rc,zc);
      M(i) = fc'*ms;
      Mr(i) = fc'*mrs;
      Mz(i) = fc'*mzs;
    else
      [ms, mrs, mzs] = mutinds(rc,zc,r(i),z(i));
      M(i) = fc'*ms;
      Mr(i) = fc'*mrs;
      Mz(i) = fc'*mzs;
    end
  else
    M(i) = fc'*mutinds(rc,zc,r(i),z(i));
  end
end

ongrid = r > rgmin & r < rgmax & z > zgmin & z < zgmax;
if any(ongrid(:))
  S = zeros(ngg,1);
  T = zeros(ngg);
  mu0 = 4e-7*pi;
  for i = 2:nz-1
    for j = 2:nr-1
      k = i+(j-1)*nz;
      T(k,k) = -2/dr^2-2/dz^2;
      T(k,k-nz) = 1/dr^2+1/2/rg(j)/dr;
      T(k,k+nz) = 1/dr^2-1/2/rg(j)/dr;
      T(k,k-1) = 1/dz^2;
      T(k,k+1) = 1/dz^2;
      if A(i,j) > 0
	S(k) = -2*pi*mu0*rg(j)*A(i,j)/dr/dz/Atot;
      end
    end
  end
  for i = 1:nz
    for j = [1 nr]
      k = i+(j-1)*nz;
      T(k,k) = 1;
      S(k) = fc'*mutinds(rc,zc,rg(j),zg(i));
    end
  end
  for i = [1 nz]
    for j = 2:nr-1
      k = i+(j-1)*nz;
      T(k,k) = 1;
      S(k) = fc'*mutinds(rc,zc,rg(j),zg(i));
    end
  end
  yzr = reshape(T\S,nz,nr);
  if nargout > 1
    [M(ongrid(:)), Mr(ongrid(:)), Mz(ongrid(:))] = ...
      gs_interp2(rg,zg,yzr,r(ongrid(:)),z(ongrid(:)));
  else
    M(ongrid(:)) = gs_interp2(rg,zg,yzr,r(ongrid(:)),z(ongrid(:)));
  end
end


return

% Tests

M = A(:)'*yzr(:)/Atot
M = fc'*gs_interp2(rg,zg,yzr,rc,zc)

disp(['Testing total area, relative discrepancy: ' ...
  num2str((sum(A(:))-polyarea(rp,zp))/Atot)])

yzr = reshape(T\S,nz,nr);
Me = zeros(nz,nr);
for i = [1 nz]
  for j = 1:nr
    Me(i,j) = fc'*mutinds(rc,zc,rg(j),zg(i));
  end
end
for i = 1:nz
  for j = [1 nr]
    Me(i,j) = fc'*mutinds(rc,zc,rg(j),zg(i));
  end
end
iedge = Me ~= 0;
disp(['Testing flux at grid edge, maximum relative discrepancy: ' ...
  num2str(max(abs((yzr(iedge)-Me(iedge))./Me(iedge))))])
