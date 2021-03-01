%  USAGE:   gs_cell_coverage
%
%  PURPOSE: Calculate how much area of each grid cell is covered by plasma
%           A "cell" is a rectangle with a coordinate rg,zg at its center
%
%  INPUTS: data for points where the boundary intersects cell edges:
%            redge, zedge, nedge, iedge, fedge, xedge (see gs_trace_edge)
%          data for boundary points:
%            rbbbs, zbbbs, (see gs_trace_boundary)
%          Other precomputed variables:	
%            rgg, zgg, dr, dz, nr, nz (grid variables)
%
%  OUTPUTS: Acell, plasma-covered area of each cell on the grid
%	    RAcell, surface integral of R over the plasma in cell
%           ARcell, surface integral of 1/R over plasma in cell
%
%  METHOD:  Areas of polygons are computed by integrating R*dZ from corner to corner
%           going ccw around the cell. When integrating along boundary: R(x)*Z'(x)dx
	
%  VERSION @(#)gs_cell_coverage.m	1.2 03/13/14
%
%  WRITTEN BY:  Anders Welander  ON	3/12/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncell = zeros(nz,nr); % Number of entry and exit points
fcell = zeros(ngg,6); % How grid index (ig) changes when exiting cell, going ccw
rcell = zeros(ngg,6); % The r of entry and exit points
zcell = zeros(ngg,6); % The z of entry and exit points

% Area of plasma within grid cell
Acell = zeros(nz,nr);

% Surface integral of R over plasma within grid cell
RAcell = zeros(nz,nr);

% Surface integral of 1/R over plasma within grid cell
ARcell = zeros(nz,nr);

for j = 1:nedge
  ig = gedge(j); % Entering grid cell ig:
  if ncell(ig) > 0 % Calculate areas from previous exit (a) to new entry (b)
    ri = rgg(ig)-dr/2;
    ro = rgg(ig)+dr/2;
    zd = zgg(ig)-dz/2;
    zu = zgg(ig)+dz/2;
    za = zcell(ig,ncell(ig));
    zb = zedge(j);
    fa = fcell(ig,ncell(ig));
    fb = fedge(j);
    if fa == -nz % boundary went out through inner edge of the cell
      if fb == nz % boundary came in through inner edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ri;
        RAcell(ig) = RAcell(ig) + (zb-za)*ri^2/2;
        ARcell(ig) = ARcell(ig) + (zb-za)*log(ri);
      else % Lower inner corner must be inside the plasma
        if ncell(ig-1) == 0
	  Acell(ig-1) = Ag;
	end
        if ncell(ig-1-nz) == 0
	  Acell(ig-1-nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zd-za)*ri;
        RAcell(ig) = RAcell(ig) + (zd-za)*ri^2/2;
        ARcell(ig) = ARcell(ig) + (zd-za)*log(ri);
	if fb == -nz % boundary came in from outer edge
	   Acell(ig) =  Acell(ig) + (zb-zd)*ro;
	  RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
	  ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	elseif fb == -1 % boundary came in from top
	   Acell(ig) =  Acell(ig) + dz*ro;
	  RAcell(ig) = RAcell(ig) + dz*ro^2/2;
	  ARcell(ig) = ARcell(ig) + dz*log(ro);
	end
      end
    elseif fa == nz % boundary went out through outer edge of the cell
      if fb == -nz % boundary came in through outer edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ro;
        RAcell(ig) = RAcell(ig) + (zb-za)*ro^2/2;
        ARcell(ig) = ARcell(ig) + (zb-za)*log(ro);
      else % Upper outer corner must be inside the plasma
        if ncell(ig+1) == 0
	  Acell(ig+1) = Ag;
	end
        if ncell(ig+1+nz) == 0
	  Acell(ig+1+nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zu-za)*ro;
        RAcell(ig) = RAcell(ig) + (zu-za)*ro^2/2;
        ARcell(ig) = ARcell(ig) + (zu-za)*log(ro);
	if fb == nz % boundary came in from inner edge
	   Acell(ig) =  Acell(ig) + (zb-zu)*ri;
	  RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
	  ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	elseif fb ==  1
	   Acell(ig) =  Acell(ig) - dz*ri;
	  RAcell(ig) = RAcell(ig) - dz*ri^2/2;
	  ARcell(ig) = ARcell(ig) - dz*log(ri);
	end
      end
    elseif fa == -1 % boundary went out through bottom of the cell
      if fb ~= 1 % boundary didn't come in through bottom of the cell
	% Lower outer corner must be inside the plasma
        if ncell(ig-1) == 0
	  Acell(ig-1) = Ag;
	end
        if ncell(ig-1+nz) == 0
	  Acell(ig-1+nz) = Ag;
	end
	if fb == -nz % boundary came in through outer edge of the cell
           Acell(ig) =  Acell(ig) + (zb-zd)*ro;
          RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
          ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	else
           Acell(ig) =  Acell(ig) + dz*ro;
          RAcell(ig) = RAcell(ig) + dz*ro^2/2;
          ARcell(ig) = ARcell(ig) + dz*log(ro);
	  if fb == nz % boundary came in from inner edge
             Acell(ig) =  Acell(ig) + (zb-zu)*ri;
            RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
            ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	  end
	end
      end
    elseif fa == 1 % boundary went out through top of the cell
      if fb ~= -1 % boundary didn't come in through top of the cell
	% Upper inner corner must be inside the plasma
        if ncell(ig+1) == 0
	  Acell(ig+1) = Ag;
	end
        if ncell(ig+1-nz) == 0
	  Acell(ig+1-nz) = Ag;
	end
	if fb == nz % boundary came in from inner edge
           Acell(ig) =  Acell(ig) + (zb-zu)*ri;
          RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
          ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	else
           Acell(ig) =  Acell(ig) - dz*ri;
          RAcell(ig) = RAcell(ig) - dz*ri^2/2;
          ARcell(ig) = ARcell(ig) - dz*log(ri);
	  if fb == -nz % boundary came in through outer edge of the cell
             Acell(ig) =  Acell(ig) + (zb-zd)*ro;
            RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
            ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	  end
	end
      end
    end    
  end
  ncell(ig) = ncell(ig)+1;
  rcell(ig,ncell(ig)) = redge(j);
  zcell(ig,ncell(ig)) = zedge(j);
  fcell(ig,ncell(ig)) = fedge(j);
  xa = xedge(j);
  k = iedge(j); % k is index to bbbs point
  keep_going = true;
  while keep_going % scan iedge(j):iedge(j+1)
    if k == iedge(j+1)
      xb = xedge(j+1);
    else
      xb = 1;
    end
    % The coefficients of r(x)*dz(x):
    xx(1) = rbbbs(k)*DY(k);
    xx(2) = DX(k)*DY(k);
    ya = 1;
    yb = 1;
    dA = 0;
    for i = 1:2
      ya = ya*xa;
      yb = yb*xb;
      dA = dA + xx(i)*(yb-ya)/i;
    end
    Acell(ig) = Acell(ig) + dA;
    % The coefficients of r(x)*r(x)*dz(x):
    x1(3) = DX(k)*xx(2);
    x1(2) = DX(k)*xx(1)+rbbbs(k)*xx(2);
    x1(1) = rbbbs(k)*xx(1);
    ya = 1;
    yb = 1;
    dRA = 0;
    for i = 1:3
      ya = ya*xa;
      yb = yb*xb;
      dRA = dRA + x1(i)*(yb-ya)/i/2;
    end
    RAcell(ig) = RAcell(ig) + dRA;
    if DX(k) ~= 0
      ra = rbbbs(k)+DX(k)*xa;
      rb = rbbbs(k)+DX(k)*xb;
      dAR = DY(k)*(rb*log(rb)-ra*log(ra)-DX(k)*(xb-xa))/DX(k);
    else
      dAR = DY(k)*log(rbbbs(k))*(xb-xa);
    end
    ARcell(ig) = ARcell(ig) + dAR;
    xa = 0; % If the k loop continues with next bbbs point then xa will be 0
    if k == iedge(j+1)
      keep_going = false;
    else
      k = k+1;
      if k == nbbbs
	k = 1;
      end
    end
  end
  % Exiting cell ig:
  ncell(ig) = ncell(ig)+1;
  rcell(ig,ncell(ig)) = redge(j+1);
  zcell(ig,ncell(ig)) = zedge(j+1);
  fcell(ig,ncell(ig)) = fedge(j+1);
end


% Walk ccw along the edges of the cells from the exit point (a) to the entry (b)
for ig = 1:ngg
  if ncell(ig) > 0
    ri = rgg(ig)-dr/2;
    ro = rgg(ig)+dr/2;
    zd = zgg(ig)-dz/2;
    zu = zgg(ig)+dz/2;
    za = zcell(ig,ncell(ig));
    zb = zcell(ig,1);
    fa = fcell(ig,ncell(ig));
    fb = fcell(ig,1);
    if fa == -nz % boundary went out through inner edge of the cell
      if fb == nz % boundary came in through inner edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ri;
        RAcell(ig) = RAcell(ig) + (zb-za)*ri^2/2;
        ARcell(ig) = ARcell(ig) + (zb-za)*log(ri);
	if zb > za % In this case we just calculated what was nicked out
           Acell(ig) =  Acell(ig) + Ag;
          RAcell(ig) = RAcell(ig) + RA(ig);
          ARcell(ig) = ARcell(ig) + AR(ig);	  
	end
      else % Lower inner corner must be inside the plasma
        if ncell(ig-1) == 0
	  Acell(ig-1) = Ag;
	end
        if ncell(ig-1-nz) == 0
	  Acell(ig-1-nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zd-za)*ri;
        RAcell(ig) = RAcell(ig) + (zd-za)*ri^2/2;
        ARcell(ig) = ARcell(ig) + (zd-za)*log(ri);
	if fb == -nz % boundary came in from outer edge
	   Acell(ig) =  Acell(ig) + (zb-zd)*ro;
	  RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
	  ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	elseif fb == -1 % boundary came in from top
	   Acell(ig) =  Acell(ig) + dz*ro;
	  RAcell(ig) = RAcell(ig) + dz*ro^2/2;
	  ARcell(ig) = ARcell(ig) + dz*log(ro);
	end
      end
    elseif fa == nz % boundary went out through outer edge of the cell
      if fb == -nz % boundary came in through outer edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ro;
        RAcell(ig) = RAcell(ig) + (zb-za)*ro^2/2;
        ARcell(ig) = ARcell(ig) + (zb-za)*log(ro);
      else % Upper outer corner must be inside the plasma
        if ncell(ig+1) == 0
	  Acell(ig+1) = Ag;
	end
        if ncell(ig+1+nz) == 0
	  Acell(ig+1+nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zu-za)*ro;
        RAcell(ig) = RAcell(ig) + (zu-za)*ro^2/2;
        ARcell(ig) = ARcell(ig) + (zu-za)*log(ro);
	if fb == nz % boundary came in from inner edge
	   Acell(ig) =  Acell(ig) + (zb-zu)*ri;
	  RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
	  ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	elseif fb ==  1
	   Acell(ig) =  Acell(ig) - dz*ri;
	  RAcell(ig) = RAcell(ig) - dz*ri^2/2;
	  ARcell(ig) = ARcell(ig) - dz*log(ri);
	end
      end
    elseif fa == -1 % boundary went out through bottom of the cell
      if fb ~= 1 % boundary didn't come in through bottom of the cell
	% Lower outer corner must be inside the plasma
        if ncell(ig-1) == 0
	  Acell(ig-1) = Ag;
	end
        if ncell(ig-1+nz) == 0
	  Acell(ig-1+nz) = Ag;
	end
	if fb == -nz % boundary came in through outer edge of the cell
           Acell(ig) =  Acell(ig) + (zb-zd)*ro;
          RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
          ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	else
           Acell(ig) =  Acell(ig) + dz*ro;
          RAcell(ig) = RAcell(ig) + dz*ro^2/2;
          ARcell(ig) = ARcell(ig) + dz*log(ro);
	  if fb == nz % boundary came in from inner edge
             Acell(ig) =  Acell(ig) + (zb-zu)*ri;
            RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
            ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	  end
	end
      end
    elseif fa == 1 % boundary went out through top of the cell
      if fb ~= -1 % boundary didn't come in through top of the cell
	% Upper inner corner must be inside the plasma
        if ncell(ig+1) == 0
	  Acell(ig+1) = Ag;
	end
        if ncell(ig+1-nz) == 0
	  Acell(ig+1-nz) = Ag;
	end
	if fb == nz % boundary came in from inner edge
           Acell(ig) =  Acell(ig) + (zb-zu)*ri;
          RAcell(ig) = RAcell(ig) + (zb-zu)*ri^2/2;
          ARcell(ig) = ARcell(ig) + (zb-zu)*log(ri);
	else
           Acell(ig) =  Acell(ig) - dz*ri;
          RAcell(ig) = RAcell(ig) - dz*ri^2/2;
          ARcell(ig) = ARcell(ig) - dz*log(ri);
	  if fb == -nz % boundary came in through outer edge of the cell
             Acell(ig) =  Acell(ig) + (zb-zd)*ro;
            RAcell(ig) = RAcell(ig) + (zb-zd)*ro^2/2;
            ARcell(ig) = ARcell(ig) + (zb-zd)*log(ro);
	  end
	end
      end
    end    
  end
end

% Find the rest of the cells that must be inside the plasma
% and also find approximate solutions for ARcell of partly covered cells
for j = 1:nz-1
  for k = 1:nr-1
    if ncell(j,k) == 0 & Acell(j,k) > 0 % if true then j,k is inside the plasma
      if ncell(j,k+1) == 0 % if true then j,k+1 is also inside the plasma
	Acell(j,k+1) = Ag;
      end
      if ncell(j+1,k) == 0 % if true then j+1,k is also inside the plasma
	Acell(j+1,k) = Ag;
      end
    end
  end
end
RAcell(Acell==Ag) = RA(Acell==Ag);
ARcell(Acell==Ag) = AR(Acell==Ag);

Acell(Acell < 0) = Acell(Acell < 0)+Ag;
ARcell(ARcell < 0) = ARcell(ARcell < 0)+AR(ARcell < 0);
RAcell(RAcell < 0) = RAcell(RAcell < 0)+RA(RAcell < 0);
