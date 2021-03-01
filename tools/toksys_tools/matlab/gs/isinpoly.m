function ft = isinpoly(rt,zt,rp,zp)
%
%  USAGE:   ft = isinpoly(rt,zt,rp,zp)
%           ft = isinpoly(rt,zt)
%           Second call uses latest supplied polygon and executes faster
%
%  PURPOSE: Test if point(s) rt, zt are inside polygon rp, zp
%
%  INPUTS: rt, zt, coordinates of test points
%          rp, zp, coordinates of polygon corners
%
%  OUTPUTS: ft, flag(s): 0 = outside, 1 = inside
%
	
%  VERSION @(#)isinpoly.m	1.2 02/22/15
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
    for j = 1:2:size(polygon.rranges,2)-1 % Sorting won't work for zero ranges
      if polygon.rranges(i,j) == polygon.rranges(i,j+1)
        polygon.drdz(i,j:j+1) = sort(polygon.drdz(i,j:j+1)); % Fix for zero range
      end
    end
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

% Allow extraction of variable polygon
if ischar(rt)
  ft = eval(rt);
  return
end

% Allow only preparing polygon
if isempty(rt) && nargout == 0
  clear ft
  return
end

ft = rt ~= rt; % This way rt, zt can have any kind of dimension
for i = 1:polygon.unique.n-1
  for j = 1:length(ft(:))
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
