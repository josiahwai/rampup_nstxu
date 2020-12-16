function config = regrid(rg, zg, config0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:  config = regrid(rg, zg, config0)
%          config = regrid(nr, nz, config0)
%
%  PURPOSE: Replace the grid in config0 with the grid rg, zg
%           If first input is scalar (= nr) then:
%             rg = linspace(config0.rg(1),config0.rg(end),nr)'
%           If second input is scalar (= nz) then:
%             zg = linspace(config0.zg(1),config0.zg(end),nz)'
%
%  INPUTS: rg, radii of grid points [m]
%          zg, height of grid points [m]
%          config0, a.k.a. tok_data_struct or vac_objs
%            a structure containing (TokSys) tokamak information
%
%  OUTPUTS: config, like config0 but all grid quantities are for rg, zg
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  METHOD: To calculate mpp, a temporary grid is created with extra points
%          The flux is calculated at the boundary of this expanded grid
%          using mutind and then the flux at all grid points are found by
%          the FD method. The rest of the grid objects are obtained by 
%          interpolation (for now) from old grid to new
	
%
%  WRITTEN BY:  Anders Welander  ON	5/26/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rg','var')
  rg = [];
end
if ~exist('zg','var')
  zg = [];
end
if ~exist('config0','var')
  config0 = [];
end

% Check for situations that mean no-go
if length(rg) < 2 & ~isfield(config0,'rg')
  error('radial grid coordinates must be supplied in rg or config0')
end
if length(zg) < 2 & ~isfield(config0,'zg')
  error('vertical grid coordinates must be supplied in zg or config0')
end

if length(rg(:)) == 1
  nr = rg;
  rg = linspace(config0.rg(1),config0.rg(end),nr)';
end
if length(zg(:)) == 1
  nz = zg;
  zg = linspace(config0.zg(1),config0.zg(end),nz)';
end
if length(rg(:)) == 0
  nr = length(config0.rg);
  rg = linspace(config0.rg(1),config0.rg(end),nr)';
end
if length(zg(:)) == 0
  nz = length(config0.zg);
  zg = linspace(config0.zg(1),config0.zg(end),nz)';
end
nr = length(rg(:));
nz = length(zg(:));

% Check if requested grid is identical to the one in config0
if isfield(config0,'rg') & isfield(config0,'zg')
  if length(rg(:)) == length(config0.rg(:)) & ...
     length(zg(:)) == length(config0.zg(:))
    DRmax = max( rg(nr)-rg(1), config0.rg(nr)-config0.rg(1) );
    DZmax = max( zg(nz)-zg(1), config0.zg(nz)-config0.zg(1) );
    if max(abs((rg(:)-config0.rg(:))/DRmax)) < 1e-6 & ...
       max(abs((zg(:)-config0.zg(:))/DZmax)) < 1e-6
      config = config0;
      return
    end
  end
end


mu0 = 4e-7*pi;

% Number of subdivisions of grid cells at edge
nx = 1;
ny = 1;

ngg = nr*nz;

dr = (rg(nr)-rg(1))/(nr-1);
dz = (zg(nz)-zg(1))/(nz-1);

ddr = dr/nx;
ddz = dz/ny;
mdr = reshape(  ones(ny,1)*((1:nx)-0.5)*ddr-dr/2,   nx*ny,1);
mdz = reshape(  ((1:ny)-0.5)'*ddz*ones(1,nx)-dz/2,  nx*ny,1);

rg = linspace(rg(1),rg(nr),nr);
zg = linspace(zg(1),zg(nz),nz);

% Number of extra points on sides of the grid
nxi = min(2,floor(rg(1)/dr-1)); % Extra points inside
nxo = 2; % Extra points outside
nxl = 2; % Extra points below
nxu = 0; % Extra points above

rgx = rg(1)+dr*(-nxi:nr+nxo-1); 
zgx = zg(1)+dz*(-nxl:nz+nxu-1);

nrx = length(rgx);
nzx = length(zgx);
nggx = nrx*nzx;

T = zeros(nggx);
%T = sparse([],[],[],nggx,nggx,nggx*6);
ii = reshape((nxl+1:nz+nxl)'*ones(1,nr)+ones(nz,1)*(nxi:nr+nxi-1)*nzx,ngg,1);

for i = 1:nzx
  for j = 1:nrx
    k = i+(j-1)*nzx;
    T(k,k) = -2/dr^2-2/dz^2;
    if j > 1
      T(k,k-nzx) = 1/dr^2+1/2/rgx(j)/dr;
    end
    if j < nrx
      T(k,k+nzx) = 1/dr^2-1/2/rgx(j)/dr;
    end
    if i > 1
      T(k,k-1) = 1/dz^2;
    end
    if i < nzx
      T(k,k+1) = 1/dz^2;
    end
  end
end

Ti = inv(T);
a = 1e-6/dz^2/nx/ny;
b = 1e-6/dz^2/nx/ny;
c = 1e-6*(1/dr^2+1/2/dr/rgx(1))/nx/ny;
d = 1e-6*(1/dr^2-1/2/dr/rgx(nrx))/nx/ny;
ia = 1:nzx:nggx;
ib = nzx:nzx:nggx;
ic = 1:nzx;
id = (nrx-1)*nzx+ic;
rgxa = ones(nx*ny,1)*rgx + mdr*ones(1,nrx);
rgxb = ones(nx*ny,1)*rgx + mdr*ones(1,nrx);
rgxc = rgx(1)-dr         + mdr*ones(1,nzx);
rgxd = rgx(nrx)+dr       + mdr*ones(1,nzx);
zgxa = zgx(1)-dz         + mdz*ones(1,nrx);
zgxb = zgx(nzx)+dz       + mdz*ones(1,nrx);
zgxc = ones(nx*ny,1)*zgx + mdz*ones(1,nzx);
zgxd = ones(nx*ny,1)*zgx + mdz*ones(1,nzx);
if nx*ny > 1
  for j = 1:nr
    S = zeros(nggx,1);
    S(ia) =       - sum(mutind(rgxa, zgxa, rg(j), zg(1)))'*a;
    S(ib) =       - sum(mutind(rgxa, zgxb, rg(j), zg(1)))'*b;
    S(ic) = S(ic) - sum(mutind(rgxc, zgxc, rg(j), zg(1)))'*c;
    S(id) = S(id) - sum(mutind(rgxd, zgxd, rg(j), zg(1)))'*d;
    S(nxl+1+(j+nxi-1)*nzx) = S(nxl+1+(j+nxi-1)*nzx)-2*pi*mu0*rg(j)/dr/dz;
    dum = Ti*S;
    mpp(:,j) = dum(ii);
  end
else
  for j = 1:nr
    S = zeros(nggx,1);
    S(ia) =       - (mutind(rgxa, zgxa, rg(j), zg(1)))'*a;
    S(ib) =       - (mutind(rgxa, zgxb, rg(j), zg(1)))'*b;
    S(ic) = S(ic) - (mutind(rgxc, zgxc, rg(j), zg(1)))'*c;
    S(id) = S(id) - (mutind(rgxd, zgxd, rg(j), zg(1)))'*d;
    S(nxl+1+(j+nxi-1)*nzx) = S(nxl+1+(j+nxi-1)*nzx)-2*pi*mu0*rg(j)/dr/dz;
    dum = Ti*S;
    mpp(:,j) = dum(ii);
  end
end

config = config0;
config.mpp = mpp;
config.rg = rg(:);
config.zg = zg(:);
config.nr = nr;
config.nz = nz;
config.rgg = ones(nz,1)*rg;
config.zgg = zg'*ones(1,nr);
if isfield(config,'mpc')
  nc = size(config0.mpc,2);
  config.mpc = zeros(ngg,nc);
  for i = 1:nc
    config.mpc(:,i) = gs_interp2(config0.rg,config0.zg,config0.mpc(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'mpv')
  nv = size(config0.mpv,2);
  config.mpv = zeros(ngg,nv);
  for i = 1:nv
    config.mpv(:,i) = gs_interp2(config0.rg,config0.zg,config0.mpv(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'mpl')
  nl = size(config0.mpl,2);
  config.mpl = zeros(ngg,nl);
  for i = 1:nl
    config.mpl(:,i) = gs_interp2(config0.rg,config0.zg,config0.mpl(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'mph')
  nh = size(config0.mph,2);
  config.mph = zeros(ngg,nh);
  for i = 1:nh
    config.mph(:,i) = gs_interp2(config0.rg,config0.zg,config0.mph(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'mpt')
  nt = size(config0.mpt,2);
  config.mpt = zeros(ngg,nt);
  for i = 1:nt
    config.mpt(:,i) = gs_interp2(config0.rg,config0.zg,config0.mpt(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'gpb')
  nb = size(config0.gpb,2);
  config.gpb = zeros(ngg,nb);
  for i = 1:nb
    config.gpb(:,i) = gs_interp2(config0.rg,config0.zg,config0.gpb(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'gbr2c')
  n = size(config0.gbr2c,2);
  config.gbr2c = zeros(ngg,n);
  for i = 1:n
    config.gbr2c(:,i) = gs_interp2(config0.rg,config0.zg,config0.gbr2c(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'gbz2c')
  n = size(config0.gbz2c,2);
  config.gbz2c = zeros(ngg,n);
  for i = 1:n
    config.gbz2c(:,i) = gs_interp2(config0.rg,config0.zg,config0.gbz2c(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'gbr2v')
  n = size(config0.gbr2v,2);
  config.gbr2v = zeros(ngg,n);
  for i = 1:n
    config.gbr2v(:,i) = gs_interp2(config0.rg,config0.zg,config0.gbr2v(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'gbz2v')
  n = size(config0.gbz2v,2);
  config.gbz2v = zeros(ngg,n);
  for i = 1:n
    config.gbz2v(:,i) = gs_interp2(config0.rg,config0.zg,config0.gbz2v(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'gbr2p')
  n = size(config0.gbr2p,2);
  config.gbr2p = zeros(ngg,n);
  for i = 1:n
    config.gbr2p(:,i) = gs_interp2(config0.rg,config0.zg,config0.gbr2p(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'gbz2p')
  n = size(config0.gbz2p,2);
  config.gbz2p = zeros(ngg,n);
  for i = 1:n
    config.gbz2p(:,i) = gs_interp2(config0.rg,config0.zg,config0.gbz2p(:,i), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'gmsebrp')
  n = size(config0.gmsebrp,1);
  config.gmsebrp = zeros(n,ngg);
  for i = 1:n
    config.gmsebrp(i,:) = gs_interp2(config0.rg,config0.zg,config0.gmsebrp(i,:), ...
      config.rgg(:),config.zgg(:));
  end
end
if isfield(config,'gmsebzp')
  n = size(config0.gmsebzp,1);
  config.gmsebzp = zeros(n,ngg);
  for i = 1:n
    config.gmsebzp(i,:) = gs_interp2(config0.rg,config0.zg,config0.gmsebzp(i,:), ...
      config.rgg(:),config.zgg(:));
  end
end
