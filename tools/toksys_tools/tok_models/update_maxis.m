function [rmaxis, zmaxis, maxis, iia, wa, drmaxisdpsi, dzmaxisdpsi] = ...
           update_maxis(rmaxis, zmaxis, psizr, config)
%%
%  USAGE: Initial call must include config:
%         [rmaxis, zmaxis, maxis, iia, wa, drmaxisdpsi, dzmaxisdpsi] = ...
%           update_maxis(rmaxis, zmaxis, psizr, config)
%
%         Initial call that only configures:
%         update_maxis([], [], [], config)
%
%         Repeat calls execute faster if config is omitted:
%         [rmaxis, zmaxis, maxis, iia, wa, drmaxisdpsi, dzmaxisdpsi] = ...
%           update_maxis(rmaxis, zmaxis, psizr)
%
%  PURPOSE: Zoom in on magnetic axis
%
%  INPUTS: rmaxis, zmaxis, approximate axis position
%          psizr, flux on the grid
%          config, a structure with fields:
%            rg, zg, grid coordinates
%            limdata, limiter coordinates [Zlim; Rlim]
%
%  OUTPUTS: rmaxis, zmaxis, updated position of magnetic axis
%           maxis, boolean true if axis still inside limiter after update
%           iia, wa, indices and weights such that psimag = wa*psizr(iia)
%           drmaxisdpsi, weights so that drmaxis = drmaxisdpsi*dpsizr(iia)
%           dzmaxisdpsi, weights so that dzmaxis = dzmaxisdpsi*dpsizr(iia)
%	
%  METHOD: interpolation with bicubic Hermite splines, Newton-Rhapson to zoom,
%          private version of isinpoly to check that axis is inside limiter
	
%  VERSION @(#)update_maxis.m	1.2 07/13/14
%
%  WRITTEN BY:  Anders Welander  ON	7/13/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent rg zg nr nz dr dz mx neighbors limdata

% Configure for new grid and/or limiter if config is supplied
if exist('config','var')
  missing_fields = '';
  if isfield(config,'rg')
    rg = config.rg;
  else
    if isempty(rg)
      missing_fields = 'rg, ';
    end
  end
  if isfield(config,'zg')
    zg = config.zg;
  else
    if isempty(zg)
      missing_fields = [missing_fields 'zg, '];
    end
  end
  if isfield(config,'limdata')
    limdata = config.limdata;
  else
    if isempty(limdata)
      missing_fields = [missing_fields 'limdata, '];
    end
  end
  if ~isempty(missing_fields)
    error(['Fields ' missing_fields(1:end-2) ' missing from config'])
  end
  
  % Parse the limiter in limdata. Tolerate many variations in the format
  if size(limdata,2) > size(limdata,1)
    limdata = limdata';
  end
  if min(limdata(:,1)) < min(limdata(:,2))
    limdata = limdata(:,[2 1]);
  end
  if limdata(1,1) == limdata(end,1) & limdata(1,2) == limdata(end,2)
    Rlim = limdata(:,1);
    Zlim = limdata(:,2);
  else
    Rlim = [limdata(:,1); limdata(1,1)];
    Zlim = [limdata(:,2); limdata(1,2)];
  end
  nlim = length(Rlim);
  ilim = logical(ones(1,nlim)); % Flag for whether to include point
  % Remove points that lie on a straight line
  for i = 1:nlim-2
    v1 = [Rlim(i+1)-Rlim(i  ), Zlim(i+1)-Zlim(i  )];
    v2 = [Rlim(i+2)-Rlim(i+1), Zlim(i+2)-Zlim(i+1)];
    dum = v1*v2'/norm(v1)/norm(v2);
    if abs(dum-1) < 1e-9
      ilim(i+1) = false;
    end
    if abs(dum+1) < 1e-9 % The limiter data for DIII-D has this error (as of 2014)
      ilim(i+2) = false;
    end
  end
  nlim2 = sum(ilim);

  % Prepare isinpoly
  isinpoly_update_maxis([],[],Rlim(ilim),Zlim(ilim))

  % Variables for cubic interpolation on the grid
  nr = length(rg);
  nz = length(zg);
  dr = (rg(nr)-rg(1))/(nr-1);
  dz = (zg(nz)-zg(1))/(nz-1);
  mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % Used for interpolation
  ir4 = (-1:2)'*nz; % For finding flux on horizontal grid lines
  iz4 = (-1:2)';    % For finding flux on vertical grid lines
  ir8 = [iz4-nz;iz4+nz]; % For finding derivatives w.r.t. r on vertical grid lines
  iz8 = [ir4- 1;ir4+ 1]; % For finding derivatives w.r.t. z on horizontal grid lines
  neighbors = reshape(iz4*ones(1,4)+ones(4,1)*ir4',1,16);
  
  if nargout == 0
    clear rmaxis
    return
  end
  
else
  if nargout == 0 && nargin < 3
    clear rmaxis
    disp(' ')
    disp('********  update_maxis is a function to update magnetic axis ********')
    help update_maxis
    return
  end
end

% Done with configuration


% Search for magnetic axis with Newton-Rhapson
twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
iterations = 0;
while iterations < 19 && twopirbrzmax > 1e-10
  iz = (zmaxis-zg(1))/dz + 1;
  i = min(nz-2,max(2,floor(iz)));
  jr = (rmaxis-rg(1))/dr + 1;
  j = min(nr-2,max(2,floor(jr)));
  k = i + (j-1)*nz;
  iia = k + neighbors'; % iia indexes 16 points around magnetic axis
  pp = psizr(iia);
  tr = jr - j;
  tz = iz - i;
  wa = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16); % weights to get flux at axis
  war = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16); % weights to get dpsi/dr
  waz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16); % weights to get dpsi/dz
  warr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16); % weights to get d2psi/dr2
  wazz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16); % weights to get d2psi/dr2
  warz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16); % weights to get d2psi/drdz
  ashift_Ba = -inv([warr*pp warz*pp; warz*pp wazz*pp]);
  cmaxis = ashift_Ba*[war*pp; waz*pp]; % Correction of position
  dum = 4*((cmaxis(1)/dr)^2+(cmaxis(2)/dz)^2);
  if dum > 1
    cmaxis = cmaxis/sqrt(dum);
  end
  rmaxis = rmaxis+cmaxis(1);
  zmaxis = zmaxis+cmaxis(2);
  iterations = iterations+1;
  twopirbrzmax = max(abs([war*pp waz*pp]));
end

drmaxisdpsi = ashift_Ba(1,1)*war+ashift_Ba(1,2)*waz;
dzmaxisdpsi = ashift_Ba(2,1)*war+ashift_Ba(2,2)*waz;

maxis = norm(cmaxis) < dr && isinpoly_update_maxis(rmaxis,zmaxis);

% Done updating maxis


function ft = isinpoly_update_maxis(rt,zt,rp,zp)
%
%  USAGE:   ft = isinpoly_update_maxis(rt,zt,rp,zp)
%           ft = isinpoly_update_maxis(rt,zt)
%           Second call uses latest supplied polygon and executes faster
%
%  PURPOSE: Test if point(s) rt, zt are inside polygon rp, zp
%
%  INPUTS: rt, zt, coordinates of test points
%          rp, zp, coordinates of polygon corners
%
%  OUTPUTS: ft, flag(s): 0 = outside, 1 = inside
%
	
%  VERSION @(#)isinpoly.m	1.1 07/06/14
%
%  WRITTEN BY:  Anders Welander  ON	7/6/14
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent polygon

if exist('rp','var') && exist('zp','var')

  % Last polygon corner same as first
  if rp(end) ~= rp(1) || zp(end) ~= zp(1)
    polygon.r = [rp(:)' rp(1)];
    polygon.z = [zp(:)' zp(1)];
  else
    polygon.r = rp(:)';
    polygon.z = zp(:)';
  end
  
  % Remove any duplicates of polygon corners
  ii = [1 diff(polygon.r)] | [1 diff(polygon.z)];
  polygon.r = polygon.r(ii);
  polygon.z = polygon.z(ii);
  
  polygon.n = length(polygon.z);
  
  % Unique z values
  polygon.unique.z = unique(polygon.z);
  polygon.unique.n = length(polygon.unique.z);
    
  % Find ranges in r inside polygon at unique z
  % nranges will hold maximum number of ranges
  nranges = 1;
  % rranges are valid just above unique z
  polygon.rranges = nan(polygon.unique.n,2*nranges);
  polygon.drdz = nan(polygon.unique.n,2*nranges);
  for i = 1:polygon.unique.n
    z = polygon.unique.z(i);
    for j = 1:polygon.n-1
      dz = polygon.z(j+1)-polygon.z(j);
      if polygon.z(j) <= z && z < polygon.z(j+1) || ...
         polygon.z(j) > z && z >= polygon.z(j+1)
	% interpolate to find r at z
	dr = polygon.r(j+1)-polygon.r(j);
	drdz = dr/dz;
	f = (z-polygon.z(j));
	r = polygon.r(j) + f*drdz;
	% Put r & drdz in first free entry
        k = 1;
        while ~isnan(polygon.rranges(i,k))
          k = k+1;
	  if k > 2*nranges
	    nranges = nranges + 1;
	    polygon.rranges(:,k:2*nranges) = nan;
	    polygon.drdz(:,k:2*nranges) = nan;
	  end
        end
	polygon.rranges(i,k) = r;
	polygon.drdz(i,k) = drdz;
      end % End of < = > tests
    end % End of j loop
    [polygon.rranges(i,:), kranges] = sort(polygon.rranges(i,:));
    polygon.drdz(i,:) = polygon.drdz(i,kranges);
  end % End of i loop
  polygon.nranges = nranges;
end % End of processing polygon

if 0
  % Confirmation plot
  clf
  hold on
  for i = 1:polygon.unique.n
    z2 = polygon.unique.z(i)+[0 0];
    for k = 1:polygon.nranges
      kk = 2*k-[1 0];
      r2 = polygon.rranges(i,kk);
      plot(r2,z2,'LineWidth',3)
    end
  end
  plot(polygon.r,polygon.z,'r')
end

% Allow extraction of variables, i.e. 'polygon'
if ischar(rt)
  ft = eval(rt);
  return
end

% Allow only preparing polygon
if isempty(rt) && nargout == 0
  clear ft
  return
end

ft = rt ~= rt;
[m,n,p] = size(rt);
for i = 1:polygon.unique.n-1
  for j = 1:m*n*p
    if polygon.unique.z(i) < zt(j) && zt(j) <= polygon.unique.z(i+1)
      for k = 1:polygon.nranges
        dz = zt(j) - polygon.unique.z(i);
        rmin = polygon.rranges(i,2*k-1)+polygon.drdz(i,2*k-1)*dz;
        rmax = polygon.rranges(i,2*k  )+polygon.drdz(i,2*k  )*dz;
	ft(j) = ft(j) || rmin < rt(j) && rt(j) < rmax;
      end
    end
  end
end
