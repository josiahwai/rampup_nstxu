function [y, yr, yz, yrr, yrz, yzz, yrrr, yrrz, yrzz, yzzz, yrrrr] = ...
         gs_interp2(i1, i2, i3, i4, i5, i6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   y = gs_interp2(rg, zg, yg, r, z)
%           Derivatives also available with:
%           [y, yr, yz, yrr, yrz, yzz, yrrr, yrrz, yrzz, yzzz] = ...
%           gs_interp2(rg, zg, yg, r, z)
%           To return indices and weights:
%           [i, w, wr, wz, wrr, wrz, wzz, wrrr, wrrz, wrzz, wzzz] = ...
%           gs_interp2(rg, zg, yg, r, z, 'WEIGHTS')
%           where y = sum(w'.*yg(i)'), yr = sum(wr'.*yg(i)'), etc.
%
%  PURPOSE: Interpolate values of yg on grid rg, zg to get values at r, z
%           using cubic Hermite splines (which are used by gs codes)
%           http://en.wikipedia.org/wiki/Bicubic_interpolation
%
%  INPUTS: rg, zg,  grid point coordinates, (default 1:nr, 1:nz)
%          yg,      values at the nz x nr points on the grid
%          r, z,    coordinates for interpolation points 
%
%  OUTPUTS:  y,       interpolated values at points r, z
%            yr, etc, derivative of y w.r.t. to r, etc at points r, z
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Anders Welander  ON	1/29/14
%
%  MODIFICATION HISTORY: 9/10/15 return weights when flag set			
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
  error('Provide at least matrix and coordinates for interpolation')
end

if min(size(i1)) > 1
  % i1 must be yg and rg, zg were omitted in the call
  yg = i1;
  [nz, nr] = size(yg);
  rg = 1:nr;
  zg = 1:nz;
  r = i2;
  z = i3;
  if nargin > 3
    flag = i4;
  end
else
  % i3 must be yg
  rg = i1;
  zg = i2;
  yg = i3;
  r = i4;
  z = i5;
  if nargin > 5
    flag = i6;
  end
end

mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % Used for interpolation

nr = length(rg);
nz = length(zg);
if isempty(yg)
  yg = zeros(nz,nr);
end
yg = reshape(yg, nz ,nr);

dr = (rg(nr)-rg(1))/(nr-1);
dz = (zg(nz)-zg(1))/(nz-1);

% Allow r & z to be different sizes that are both not scalar
sr = size(r);
sz = size(z);
if sr(1) == 1 & sz(1) > 1
  r = ones(sz(1),1)*r;
end
if sz(1) == 1 & sr(1) > 1
  z = ones(sr(1),1)*z;
end
if sr(2) == 1 & sz(2) > 1
  r = r*ones(1,sz(2));
end
if sz(2) == 1 & sr(2) > 1
  z = z*ones(1,sr(2));
end

n = numel(r+z);
inan = isnan(r(:)+z(:));

% Expand the grid by one row or column in each direction if needed
% This ensures that gs_interp2(rg,zg,yg,rgg,zgg) equals yg also at the edges
rgmin = rg(1);
if min(r(:)) < rg(2)
  rgmin = rgmin-dr;
  nr = nr+1;
  yg = [3*yg(:,1)-3*yg(:,2)+yg(:,3) yg];
  extrainnercol = 1;
else
  extrainnercol = 0;
end
if max(r(:)) > rg(end-1)
  nr = nr+1;
  yg = [yg 3*yg(:,end)-3*yg(:,end-1)+yg(:,end-2)];
  extraoutercol = 1;
else
  extraoutercol = 0;
end

zgmin = zg(1);
if min(z(:)) < zg(2)
  zgmin = zgmin-dz;
  nz = nz+1;
  yg = [3*yg(1,:)-3*yg(2,:)+yg(3,:); yg];
  extrarowbelow = 1;
else
  extrarowbelow = 0;
end
if max(z(:)) > zg(end-1)
  nz = nz+1;
  yg = [yg; 3*yg(end,:)-3*yg(end-1,:)+yg(end-2,:)];
  extrarowabove = 1;
else
  extrarowabove = 0;
end

kr = floor((r(:)-rgmin)/dr);
kz = floor((z(:)-zgmin)/dz);

% We will have to extrapolate for any r, z outside rg, zg
kr(kr<1) = 1;
kz(kz<1) = 1;
kr(kr>nr-3) = nr-3;
kz(kz>nz-3) = nz-3;

k = nz*kr+kz+1;
p = ones(n,1);
o = zeros(n,1);

x = (r(:)-rgmin-kr*dr)/dr;
wr0 = [p x x.^2 x.^3  ]*mx;
if nargout > 1
wr1 = [o p 2*x  3*x.^2]*mx/dr;
end
if nargout > 3
wr2 = [o o 2*p  6*x   ]*mx/dr^2;
end
if nargout > 6
wr3 = [o o o    6*p   ]*mx/dr^3;
end

x = (z(:)-zgmin-kz*dz)/dz;
wz0 = [p x x.^2 x.^3  ]*mx;
if nargout > 2
wz1 = [o p 2*x  3*x.^2]*mx/dz;
end
if nargout > 4
wz2 = [o o 2*p  6*x   ]*mx/dz^2;
end
if nargout > 9
wz3 = [o o o    6*p   ]*mx/dz^3;
end

i = [k-(nz+1)   k-nz     k-(nz-1)   k-(nz-2) ...
     k-1        k        k+1        k+2      ...
     k+(nz-1)   k+nz     k+(nz+1)   k+(nz+2) ...
     k+(2*nz-1) k+(2*nz) k+(2*nz+1) k+(2*nz+2)];
i(inan,1:16) = 1;

wr = wr0;
wz = wz0;
W = ...
  [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
   wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
   wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
   wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
y = r+z;
y(:) = sum(W'.*yg(i)');

if nargout > 1
  wr = wr1;
  wz = wz0;
  Wr = ...
    [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
     wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
     wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
     wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
  yr = r+z;
  yr(:) = sum(Wr'.*yg(i)');
end

if nargout > 2
  wr = wr0;
  wz = wz1;
  Wz = ...
    [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
     wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
     wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
     wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
  yz = r+z;
  yz(:) = sum(Wz'.*yg(i)');
end

if nargout > 3
  wr = wr2;
  wz = wz0;
  Wrr = ...
    [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
     wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
     wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
     wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
  yrr = r+z;
  yrr(:) = sum(Wrr'.*yg(i)');
end

if nargout > 4
  wr = wr1;
  wz = wz1;
  Wrz = ...
    [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
     wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
     wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
     wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
  yrz = r+z;
  yrz(:) = sum(Wrz'.*yg(i)');
end

if nargout > 5
  wr = wr0;
  wz = wz2;
  Wzz = ...
    [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
     wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
     wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
     wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
  yzz = r+z;
  yzz(:) = sum(Wzz'.*yg(i)');
end

if nargout > 6
  wr = wr3;
  wz = wz0;
  Wrrr = ...
    [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
     wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
     wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
     wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
  yrrr = r+z;
  yrrr(:) = sum(Wrrr'.*yg(i)');
end

if nargout > 7
  wr = wr2;
  wz = wz1;
  Wrrz = ...
    [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
     wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
     wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
     wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
  yrrz = r+z;
  yrrz(:) = sum(Wrrz'.*yg(i)');
end

if nargout > 8
  wr = wr1;
  wz = wz2;
  Wrzz = ...
    [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
     wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
     wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
     wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
  yrzz = r+z;
  yrzz(:) = sum(Wrzz'.*yg(i)');
end

if nargout > 9
  wr = wr0;
  wz = wz3;
  Wzzz = ...
    [wz(:,1).*wr(:,1) wz(:,2).*wr(:,1) wz(:,3).*wr(:,1) wz(:,4).*wr(:,1) ...
     wz(:,1).*wr(:,2) wz(:,2).*wr(:,2) wz(:,3).*wr(:,2) wz(:,4).*wr(:,2) ...
     wz(:,1).*wr(:,3) wz(:,2).*wr(:,3) wz(:,3).*wr(:,3) wz(:,4).*wr(:,3) ...
     wz(:,1).*wr(:,4) wz(:,2).*wr(:,4) wz(:,3).*wr(:,4) wz(:,4).*wr(:,4)];
  yzzz = r+z;
  yzzz(:) = sum(Wzzz'.*yg(i)');
end

if exist('flag','var') & upper(flag(1)) == 'W'
  % Assign indices to first output and weights to the rest
  iz = 1+mod(i(:,6)-1,nz);
  ir = ceil(i(:,6)/nz);
  if extrainnercol
    nr = nr-1;
    ir = ir-1;
    f = ir == 1;
    ir(f) = 2;
    if nargout > 1
      w = W(f,1:4);
      W(f,1:4) = W(f,5:8)+3*w;
      W(f,5:8) = W(f,9:12)-3*w;
      W(f,9:12) = W(f,13:16)+w;
      W(f,13:16) = 0;
    end
    if nargout > 2
      w = Wr(f,1:4);
      Wr(f,1:4) = Wr(f,5:8)+3*w;
      Wr(f,5:8) = Wr(f,9:12)-3*w;
      Wr(f,9:12) = Wr(f,13:16)+w;
      Wr(f,13:16) = 0;
    end
    if nargout > 3
      w = Wz(f,1:4);
      Wz(f,1:4) = Wz(f,5:8)+3*w;
      Wz(f,5:8) = Wz(f,9:12)-3*w;
      Wz(f,9:12) = Wz(f,13:16)+w;
      Wz(f,13:16) = 0;
    end
    if nargout > 4
      w = Wrr(f,1:4);
      Wrr(f,1:4) = Wrr(f,5:8)+3*w;
      Wrr(f,5:8) = Wrr(f,9:12)-3*w;
      Wrr(f,9:12) = Wrr(f,13:16)+w;
      Wrr(f,13:16) = 0;
    end
    if nargout > 5
      w = Wrz(f,1:4);
      Wrz(f,1:4) = Wrz(f,5:8)+3*w;
      Wrz(f,5:8) = Wrz(f,9:12)-3*w;
      Wrz(f,9:12) = Wrz(f,13:16)+w;
      Wrz(f,13:16) = 0;
    end
    if nargout > 6
      w = Wzz(f,1:4);
      Wzz(f,1:4) = Wzz(f,5:8)+3*w;
      Wzz(f,5:8) = Wzz(f,9:12)-3*w;
      Wzz(f,9:12) = Wzz(f,13:16)+w;
      Wzz(f,13:16) = 0;
    end
    if nargout > 7
      w = Wrrr(f,1:4);
      Wrrr(f,1:4) = Wrrr(f,5:8)+3*w;
      Wrrr(f,5:8) = Wrrr(f,9:12)-3*w;
      Wrrr(f,9:12) = Wrrr(f,13:16)+w;
      Wrrr(f,13:16) = 0;
    end
    if nargout > 8
      w = Wrrz(f,1:4);
      Wrrz(f,1:4) = Wrrz(f,5:8)+3*w;
      Wrrz(f,5:8) = Wrrz(f,9:12)-3*w;
      Wrrz(f,9:12) = Wrrz(f,13:16)+w;
      Wrrz(f,13:16) = 0;
    end
    if nargout > 9
      w = Wrzz(f,1:4);
      Wrzz(f,1:4) = Wrzz(f,5:8)+3*w;
      Wrzz(f,5:8) = Wrzz(f,9:12)-3*w;
      Wrzz(f,9:12) = Wrzz(f,13:16)+w;
      Wrzz(f,13:16) = 0;
    end
    if nargout > 10
      w = Wzzz(f,1:4);
      Wzzz(f,1:4) = Wzzz(f,5:8)+3*w;
      Wzzz(f,5:8) = Wzzz(f,9:12)-3*w;
      Wzzz(f,9:12) = Wzzz(f,13:16)+w;
      Wzzz(f,13:16) = 0;
    end
  end
  if extraoutercol
    nr = nr-1;
    f = ir == nr-1;
    ir(f) = nr-2; 
    if nargout > 1
      w = W(f,13:16);
      W(f,13:16) = W(f,9:12)+3*w;
      W(f,9:12) = W(f,5:8)-3*w;
      W(f,5:8) = W(f,1:4)+w;
      W(f,1:4) = 0;
    end
    if nargout > 2
      w = Wr(f,13:16);
      Wr(f,13:16) = Wr(f,9:12)+3*w;
      Wr(f,9:12) = Wr(f,5:8)-3*w;
      Wr(f,5:8) = Wr(f,1:4)+w;
      Wr(f,1:4) = 0;
    end
    if nargout > 3
      w = Wz(f,13:16);
      Wz(f,13:16) = Wz(f,9:12)+3*w;
      Wz(f,9:12) = Wz(f,5:8)-3*w;
      Wz(f,5:8) = Wz(f,1:4)+w;
      Wz(f,1:4) = 0;
    end
    if nargout > 4
      w = Wrr(f,13:16);
      Wrr(f,13:16) = Wrr(f,9:12)+3*w;
      Wrr(f,9:12) = Wrr(f,5:8)-3*w;
      Wrr(f,5:8) = Wrr(f,1:4)+w;
      Wrr(f,1:4) = 0;
    end
    if nargout > 5
      w = Wrz(f,13:16);
      Wrz(f,13:16) = Wrz(f,9:12)+3*w;
      Wrz(f,9:12) = Wrz(f,5:8)-3*w;
      Wrz(f,5:8) = Wrz(f,1:4)+w;
      Wrz(f,1:4) = 0;
    end
    if nargout > 6
      w = Wzz(f,13:16);
      Wzz(f,13:16) = Wzz(f,9:12)+3*w;
      Wzz(f,9:12) = Wzz(f,5:8)-3*w;
      Wzz(f,5:8) = Wzz(f,1:4)+w;
      Wzz(f,1:4) = 0;
    end
    if nargout > 7
      w = Wrrr(f,13:16);
      Wrrr(f,13:16) = Wrrr(f,9:12)+3*w;
      Wrrr(f,9:12) = Wrrr(f,5:8)-3*w;
      Wrrr(f,5:8) = Wrrr(f,1:4)+w;
      Wrrr(f,1:4) = 0;
    end
    if nargout > 8
      w = Wrrz(f,13:16);
      Wrrz(f,13:16) = Wrrz(f,9:12)+3*w;
      Wrrz(f,9:12) = Wrrz(f,5:8)-3*w;
      Wrrz(f,5:8) = Wrrz(f,1:4)+w;
      Wrrz(f,1:4) = 0;
    end
    if nargout > 9
      w = Wrzz(f,13:16);
      Wrzz(f,13:16) = Wrzz(f,9:12)+3*w;
      Wrzz(f,9:12) = Wrzz(f,5:8)-3*w;
      Wrzz(f,5:8) = Wrzz(f,1:4)+w;
      Wrzz(f,1:4) = 0;
    end
    if nargout > 10
      w = Wzzz(f,13:16);
      Wzzz(f,13:16) = Wzzz(f,9:12)+3*w;
      Wzzz(f,9:12) = Wzzz(f,5:8)-3*w;
      Wzzz(f,5:8) = Wzzz(f,1:4)+w;
      Wzzz(f,1:4) = 0;
    end
  end
  if extrarowbelow
    nz = nz-1;
    iz = iz-1;
    f = iz == 1;
    iz(f) = 2;
    if nargout > 1
      w = W(f,[1 5 9 13]);
      W(f,[1 5 9 13]) = W(f,[2 6 10 14])+3*w;
      W(f,[2 6 10 14]) = W(f,[3 7 11 15])-3*w;
      W(f,[3 7 11 15]) = W(f,[4 8 12 16])+w;
      W(f,[4 8 12 16]) = 0;
    end
    if nargout > 2
      w = Wr(f,[1 5 9 13]);
      Wr(f,[1 5 9 13]) = Wr(f,[2 6 10 14])+3*w;
      Wr(f,[2 6 10 14]) = Wr(f,[3 7 11 15])-3*w;
      Wr(f,[3 7 11 15]) = Wr(f,[4 8 12 16])+w;
      Wr(f,[4 8 12 16]) = 0;
    end
    if nargout > 3
      w = Wz(f,[1 5 9 13]);
      Wz(f,[1 5 9 13]) = Wz(f,[2 6 10 14])+3*w;
      Wz(f,[2 6 10 14]) = Wz(f,[3 7 11 15])-3*w;
      Wz(f,[3 7 11 15]) = Wz(f,[4 8 12 16])+w;
      Wz(f,[4 8 12 16]) = 0;
    end
    if nargout > 4
      w = Wrr(f,[1 5 9 13]);
      Wrr(f,[1 5 9 13]) = Wrr(f,[2 6 10 14])+3*w;
      Wrr(f,[2 6 10 14]) = Wrr(f,[3 7 11 15])-3*w;
      Wrr(f,[3 7 11 15]) = Wrr(f,[4 8 12 16])+w;
      Wrr(f,[4 8 12 16]) = 0;
    end
    if nargout > 5
      w = Wrz(f,[1 5 9 13]);
      Wrz(f,[1 5 9 13]) = Wrz(f,[2 6 10 14])+3*w;
      Wrz(f,[2 6 10 14]) = Wrz(f,[3 7 11 15])-3*w;
      Wrz(f,[3 7 11 15]) = Wrz(f,[4 8 12 16])+w;
      Wrz(f,[4 8 12 16]) = 0;
    end
    if nargout > 6
      w = Wzz(f,[1 5 9 13]);
      Wzz(f,[1 5 9 13]) = Wzz(f,[2 6 10 14])+3*w;
      Wzz(f,[2 6 10 14]) = Wzz(f,[3 7 11 15])-3*w;
      Wzz(f,[3 7 11 15]) = Wzz(f,[4 8 12 16])+w;
      Wzz(f,[4 8 12 16]) = 0;
    end
    if nargout > 7
      w = Wrrr(f,[1 5 9 13]);
      Wrrr(f,[1 5 9 13]) = Wrrr(f,[2 6 10 14])+3*w;
      Wrrr(f,[2 6 10 14]) = Wrrr(f,[3 7 11 15])-3*w;
      Wrrr(f,[3 7 11 15]) = Wrrr(f,[4 8 12 16])+w;
      Wrrr(f,[4 8 12 16]) = 0;
    end
    if nargout > 8
      w = Wrrz(f,[1 5 9 13]);
      Wrrz(f,[1 5 9 13]) = Wrrz(f,[2 6 10 14])+3*w;
      Wrrz(f,[2 6 10 14]) = Wrrz(f,[3 7 11 15])-3*w;
      Wrrz(f,[3 7 11 15]) = Wrrz(f,[4 8 12 16])+w;
      Wrrz(f,[4 8 12 16]) = 0;
    end
    if nargout > 9
      w = Wrzz(f,[1 5 9 13]);
      Wrzz(f,[1 5 9 13]) = Wrzz(f,[2 6 10 14])+3*w;
      Wrzz(f,[2 6 10 14]) = Wrzz(f,[3 7 11 15])-3*w;
      Wrzz(f,[3 7 11 15]) = Wrzz(f,[4 8 12 16])+w;
      Wrzz(f,[4 8 12 16]) = 0;
    end
    if nargout > 10
      w = Wzzz(f,[1 5 9 13]);
      Wzzz(f,[1 5 9 13]) = Wzzz(f,[2 6 10 14])+3*w;
      Wzzz(f,[2 6 10 14]) = Wzzz(f,[3 7 11 15])-3*w;
      Wzzz(f,[3 7 11 15]) = Wzzz(f,[4 8 12 16])+w;
      Wzzz(f,[4 8 12 16]) = 0;
    end
  end
  if extrarowabove
    nz = nz-1;
    f = iz == nz-1;
    iz(f) = nz-2;
    if nargout > 1
      w = W(f,[4 8 12 16]);
      W(f,[4 8 12 16]) = W(f,[3 7 11 15])+3*w;
      W(f,[3 7 11 15]) = W(f,[2 6 10 14])-3*w;
      W(f,[2 6 10 14]) = W(f,[1 5  9 13])+w;
      W(f,[1 5 9 13]) = 0;
    end
    if nargout > 2
      w = Wr(f,[4 8 12 16]);
      Wr(f,[4 8 12 16]) = Wr(f,[3 7 11 15])+3*w;
      Wr(f,[3 7 11 15]) = Wr(f,[2 6 10 14])-3*w;
      Wr(f,[2 6 10 14]) = Wr(f,[1 5  9 13])+w;
      Wr(f,[1 5 9 13]) = 0;
    end
    if nargout > 3
      w = Wz(f,[4 8 12 16]);
      Wz(f,[4 8 12 16]) = Wz(f,[3 7 11 15])+3*w;
      Wz(f,[3 7 11 15]) = Wz(f,[2 6 10 14])-3*w;
      Wz(f,[2 6 10 14]) = Wz(f,[1 5  9 13])+w;
      Wz(f,[1 5 9 13]) = 0;
    end
    if nargout > 4
      w = Wrr(f,[4 8 12 16]);
      Wrr(f,[4 8 12 16]) = Wrr(f,[3 7 11 15])+3*w;
      Wrr(f,[3 7 11 15]) = Wrr(f,[2 6 10 14])-3*w;
      Wrr(f,[2 6 10 14]) = Wrr(f,[1 5  9 13])+w;
      Wrr(f,[1 5 9 13]) = 0;
    end
    if nargout > 5
      w = Wrz(f,[4 8 12 16]);
      Wrz(f,[4 8 12 16]) = Wrz(f,[3 7 11 15])+3*w;
      Wrz(f,[3 7 11 15]) = Wrz(f,[2 6 10 14])-3*w;
      Wrz(f,[2 6 10 14]) = Wrz(f,[1 5  9 13])+w;
      Wrz(f,[1 5 9 13]) = 0;
    end
    if nargout > 6
      w = Wzz(f,[4 8 12 16]);
      Wzz(f,[4 8 12 16]) = Wzz(f,[3 7 11 15])+3*w;
      Wzz(f,[3 7 11 15]) = Wzz(f,[2 6 10 14])-3*w;
      Wzz(f,[2 6 10 14]) = Wzz(f,[1 5  9 13])+w;
      Wzz(f,[1 5 9 13]) = 0;
    end
    if nargout > 7
      w = Wrrr(f,[4 8 12 16]);
      Wrrr(f,[4 8 12 16]) = Wrrr(f,[3 7 11 15])+3*w;
      Wrrr(f,[3 7 11 15]) = Wrrr(f,[2 6 10 14])-3*w;
      Wrrr(f,[2 6 10 14]) = Wrrr(f,[1 5  9 13])+w;
      Wrrr(f,[1 5 9 13]) = 0;
    end
    if nargout > 8
      w = Wrrz(f,[4 8 12 16]);
      Wrrz(f,[4 8 12 16]) = Wrrz(f,[3 7 11 15])+3*w;
      Wrrz(f,[3 7 11 15]) = Wrrz(f,[2 6 10 14])-3*w;
      Wrrz(f,[2 6 10 14]) = Wrrz(f,[1 5  9 13])+w;
      Wrrz(f,[1 5 9 13]) = 0;
    end
    if nargout > 9
      w = Wrzz(f,[4 8 12 16]);
      Wrzz(f,[4 8 12 16]) = Wrzz(f,[3 7 11 15])+3*w;
      Wrzz(f,[3 7 11 15]) = Wrzz(f,[2 6 10 14])-3*w;
      Wrzz(f,[2 6 10 14]) = Wrzz(f,[1 5  9 13])+w;
      Wrzz(f,[1 5 9 13]) = 0;
    end
    if nargout > 10
      w = Wzzz(f,[4 8 12 16]);
      Wzzz(f,[4 8 12 16]) = Wzzz(f,[3 7 11 15])+3*w;
      Wzzz(f,[3 7 11 15]) = Wzzz(f,[2 6 10 14])-3*w;
      Wzzz(f,[2 6 10 14]) = Wzzz(f,[1 5  9 13])+w;
      Wzzz(f,[1 5 9 13]) = 0;
    end
  end
  k = iz+nz*(ir-1);
  i = [k-(nz+1)   k-nz     k-(nz-1)   k-(nz-2) ...
     k-1        k        k+1        k+2      ...
     k+(nz-1)   k+nz     k+(nz+1)   k+(nz+2) ...
     k+(2*nz-1) k+(2*nz) k+(2*nz+1) k+(2*nz+2)];
  i(inan,1:16) = 1;
  y = i;
  if nargout > 1
    yr = W;
  end
  if nargout > 2
    yz = Wr;
  end
  if nargout > 3
    yrr = Wz;
  end
  if nargout > 4
    yrz = Wrr;
  end
  if nargout > 5
    yzz = Wrz;
  end
  if nargout > 6
    yrrr = Wzz;
  end
  if nargout > 7
    yrrz = Wrrr;
  end
  if nargout > 8
    yrzz = Wrrz;
  end
  if nargout > 9
    yzzz = Wrzz;
  end
  if nargout > 10
    yrrrr = Wzzz;
  end
end
