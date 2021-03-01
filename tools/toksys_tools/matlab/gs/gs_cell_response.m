%  USAGE:   gs_cell_response
%
%  PURPOSE: Calculate cell coverage and how it responds to changes in rbbbs, zbbbs
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
%            dAcelldrbbbs, dAcelldzbbbs, response of  Acell to rbbbs, zbbbs
%	    dRAcelldrbbbs,dRAcelldzbbbs, response of RAcell to rbbbs, zbbbs
%	    dARcelldrbbbs,dARcelldzbbbs, response of ARcell to rbbbs, zbbbs
%
%  METHOD:  Areas of polygons are computed by integrating R*dZ from corner to corner
%           going ccw around the cell. When integrating along boundary: R(x)*Z'(x)dx
%           Does all that gs_cell_coverage does *and* also response calculations
	
%  VERSION @(#)gs_cell_response.m	1.2 03/13/14
%
%  WRITTEN BY:  Anders Welander  ON	3/21/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Area response to change in rbbbs, zbbbs
dAcelldrbbbs = zeros(ngg,nbbbs_max);
dAcelldzbbbs = zeros(ngg,nbbbs_max);

% RAcell response to rbbbs, zbbbs
dRAcelldrbbbs = zeros(ngg,nbbbs_max);
dRAcelldzbbbs = zeros(ngg,nbbbs_max);

% ARcell response to rbbbs, zbbbs
dARcelldrbbbs = zeros(ngg,nbbbs_max);
dARcelldzbbbs = zeros(ngg,nbbbs_max);

ncell = zeros(nz,nr); % Number of entry and exit points
fcell = zeros(ngg,6); % How grid index (ig) changes when exiting cell, going ccw
rcell = zeros(ngg,6); % The r of entry and exit points
zcell = zeros(ngg,6); % The z of entry and exit points

ibc = zeros(ngg,6); % The lower of the two bbbs indices around points zcell
dzcelldrbbbs0 = zeros(ngg,6); % How zcell responds to rbbbs(ibc)
dzcelldzbbbs0 = zeros(ngg,6); % How zcell responds to zbbbs(ibc)
dzcelldrbbbs1 = zeros(ngg,6); % How zcell responds to rbbbs(ibc+1)
dzcelldzbbbs1 = zeros(ngg,6); % How zcell responds to zbbbs(ibc+1)

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
    ri22 = ri^2/2;
    ro22 = ro^2/2;
    logri = log(ri);
    logro = log(ro);
    zd = zgg(ig)-dz/2;
    zu = zgg(ig)+dz/2;
    za = zcell(ig,ncell(ig));
    zb = zedge(j);
    fa = fcell(ig,ncell(ig));
    fb = fedge(j);
    k = ibc(ig,ncell(ig));
    m = iedge(j);
    dzadrbbbs0 = dzcelldrbbbs0(ig,ncell(ig));
    dzadzbbbs0 = dzcelldzbbbs0(ig,ncell(ig));
    dzadrbbbs1 = dzcelldrbbbs1(ig,ncell(ig));
    dzadzbbbs1 = dzcelldzbbbs1(ig,ncell(ig));
    dzbdrbbbs0 = dxedgedrbbbs0(j)*DY(m);
    dzbdzbbbs0 = dxedgedzbbbs0(j)*DY(m)+1-xedge(j);
    dzbdrbbbs1 = dxedgedrbbbs1(j)*DY(m);
    dzbdzbbbs1 = dxedgedzbbbs1(j)*DY(m)+xedge(j);
    if fa == -nz % boundary went out through inner edge of the cell
      if fb == nz % boundary came in through inner edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ri;
        RAcell(ig) = RAcell(ig) + (zb-za)*ri22;
        ARcell(ig) = ARcell(ig) + (zb-za)*logri;
	dAcelldrbbbs(ig,k  ) = dAcelldrbbbs(ig,k  ) - ri*dzadrbbbs0;
	dAcelldzbbbs(ig,k  ) = dAcelldzbbbs(ig,k  ) - ri*dzadzbbbs0;
	dAcelldrbbbs(ig,k+1) = dAcelldrbbbs(ig,k+1) - ri*dzadrbbbs1;
	dAcelldzbbbs(ig,k+1) = dAcelldzbbbs(ig,k+1) - ri*dzadzbbbs1;
	dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ri*dzbdrbbbs0;
	dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ri*dzbdzbbbs0;
	dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ri*dzbdrbbbs1;
	dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ri*dzbdzbbbs1;
	dRAcelldrbbbs(ig,k  ) = dRAcelldrbbbs(ig,k  ) - ri22*dzadrbbbs0;
	dRAcelldzbbbs(ig,k  ) = dRAcelldzbbbs(ig,k  ) - ri22*dzadzbbbs0;
	dRAcelldrbbbs(ig,k+1) = dRAcelldrbbbs(ig,k+1) - ri22*dzadrbbbs1;
	dRAcelldzbbbs(ig,k+1) = dRAcelldzbbbs(ig,k+1) - ri22*dzadzbbbs1;
	dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ri22*dzbdrbbbs0;
	dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ri22*dzbdzbbbs0;
	dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ri22*dzbdrbbbs1;
	dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ri22*dzbdzbbbs1;
	dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) - logri*dzadrbbbs0;
	dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - logri*dzadzbbbs0;
	dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1) - logri*dzadrbbbs1;
	dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) - logri*dzadzbbbs1;
	dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logri*dzbdrbbbs0;
	dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logri*dzbdzbbbs0;
	dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logri*dzbdrbbbs1;
	dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logri*dzbdzbbbs1;
      else % Lower inner corner must be inside the plasma
        if ncell(ig-1) == 0
	  Acell(ig-1) = Ag;
	end
        if ncell(ig-1-nz) == 0
	  Acell(ig-1-nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zd-za)*ri;
        RAcell(ig) = RAcell(ig) + (zd-za)*ri22;
        ARcell(ig) = ARcell(ig) + (zd-za)*logri;
	dAcelldrbbbs(ig,k  ) = dAcelldrbbbs(ig,k  ) - ri*dzadrbbbs0;
	dAcelldzbbbs(ig,k  ) = dAcelldzbbbs(ig,k  ) - ri*dzadzbbbs0;
	dAcelldrbbbs(ig,k+1) = dAcelldrbbbs(ig,k+1) - ri*dzadrbbbs1;
	dAcelldzbbbs(ig,k+1) = dAcelldzbbbs(ig,k+1) - ri*dzadzbbbs1;
	dRAcelldrbbbs(ig,k  ) = dRAcelldrbbbs(ig,k  ) - ri22*dzadrbbbs0;
	dRAcelldzbbbs(ig,k  ) = dRAcelldzbbbs(ig,k  ) - ri22*dzadzbbbs0;
	dRAcelldrbbbs(ig,k+1) = dRAcelldrbbbs(ig,k+1) - ri22*dzadrbbbs1;
	dRAcelldzbbbs(ig,k+1) = dRAcelldzbbbs(ig,k+1) - ri22*dzadzbbbs1;
	dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) - logri*dzadrbbbs0;
	dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - logri*dzadzbbbs0;
	dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1) - logri*dzadrbbbs1;
	dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) - logri*dzadzbbbs1;
	if fb == -nz % boundary came in from outer edge
	   Acell(ig) =  Acell(ig) + (zb-zd)*ro;
	  RAcell(ig) = RAcell(ig) + (zb-zd)*ro22;
	  ARcell(ig) = ARcell(ig) + (zb-zd)*logro;
	  dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ro*dzbdrbbbs0;
	  dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ro*dzbdzbbbs0;
	  dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ro*dzbdrbbbs1;
	  dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ro*dzbdzbbbs1;
	  dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ro22*dzbdrbbbs0;
	  dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ro22*dzbdzbbbs0;
	  dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ro22*dzbdrbbbs1;
	  dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ro22*dzbdzbbbs1;
	  dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logro*dzbdrbbbs0;
	  dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logro*dzbdzbbbs0;
	  dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logro*dzbdrbbbs1;
	  dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logro*dzbdzbbbs1;
	elseif fb == -1 % boundary came in from top
	   Acell(ig) =  Acell(ig) + dz*ro;
	  RAcell(ig) = RAcell(ig) + dz*ro22;
	  ARcell(ig) = ARcell(ig) + dz*logro;
	end
      end
    elseif fa == nz % boundary went out through outer edge of the cell
      if fb == -nz % boundary came in through outer edge of the cell
         Acell(ig) =  Acell(ig) + (zb-za)*ro;
        RAcell(ig) = RAcell(ig) + (zb-za)*ro22;
        ARcell(ig) = ARcell(ig) + (zb-za)*logro;
	dAcelldrbbbs(ig,k  ) = dAcelldrbbbs(ig,k  ) - ro*dzadrbbbs0;
	dAcelldzbbbs(ig,k  ) = dAcelldzbbbs(ig,k  ) - ro*dzadzbbbs0;
	dAcelldrbbbs(ig,k+1) = dAcelldrbbbs(ig,k+1) - ro*dzadrbbbs1;
	dAcelldzbbbs(ig,k+1) = dAcelldzbbbs(ig,k+1) - ro*dzadzbbbs1;
	dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ro*dzbdrbbbs0;
	dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ro*dzbdzbbbs0;
	dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ro*dzbdrbbbs1;
	dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ro*dzbdzbbbs1;
	dRAcelldrbbbs(ig,k  ) = dRAcelldrbbbs(ig,k  ) - ro22*dzadrbbbs0;
	dRAcelldzbbbs(ig,k  ) = dRAcelldzbbbs(ig,k  ) - ro22*dzadzbbbs0;
	dRAcelldrbbbs(ig,k+1) = dRAcelldrbbbs(ig,k+1) - ro22*dzadrbbbs1;
	dRAcelldzbbbs(ig,k+1) = dRAcelldzbbbs(ig,k+1) - ro22*dzadzbbbs1;
	dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ro22*dzbdrbbbs0;
	dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ro22*dzbdzbbbs0;
	dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ro22*dzbdrbbbs1;
	dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ro22*dzbdzbbbs1;
	dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) - logro*dzadrbbbs0;
	dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - logro*dzadzbbbs0;
	dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1) - logro*dzadrbbbs1;
	dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) - logro*dzadzbbbs1;
	dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logro*dzbdrbbbs0;
	dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logro*dzbdzbbbs0;
	dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logro*dzbdrbbbs1;
	dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logro*dzbdzbbbs1;
      else % Upper outer corner must be inside the plasma
        if ncell(ig+1) == 0
	  Acell(ig+1) = Ag;
	end
        if ncell(ig+1+nz) == 0
	  Acell(ig+1+nz) = Ag;
	end
         Acell(ig) =  Acell(ig) + (zu-za)*ro;
        RAcell(ig) = RAcell(ig) + (zu-za)*ro22;
        ARcell(ig) = ARcell(ig) + (zu-za)*logro;
	dAcelldrbbbs(ig,k  ) = dAcelldrbbbs(ig,k  ) - ro*dzadrbbbs0;
	dAcelldzbbbs(ig,k  ) = dAcelldzbbbs(ig,k  ) - ro*dzadzbbbs0;
	dAcelldrbbbs(ig,k+1) = dAcelldrbbbs(ig,k+1) - ro*dzadrbbbs1;
	dAcelldzbbbs(ig,k+1) = dAcelldzbbbs(ig,k+1) - ro*dzadzbbbs1;
	dRAcelldrbbbs(ig,k  ) = dRAcelldrbbbs(ig,k  ) - ro22*dzadrbbbs0;
	dRAcelldzbbbs(ig,k  ) = dRAcelldzbbbs(ig,k  ) - ro22*dzadzbbbs0;
	dRAcelldrbbbs(ig,k+1) = dRAcelldrbbbs(ig,k+1) - ro22*dzadrbbbs1;
	dRAcelldzbbbs(ig,k+1) = dRAcelldzbbbs(ig,k+1) - ro22*dzadzbbbs1;
	dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) - logro*dzadrbbbs0;
	dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - logro*dzadzbbbs0;
	dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1) - logro*dzadrbbbs1;
	dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) - logro*dzadzbbbs1;
	if fb == nz % boundary came in from inner edge
	   Acell(ig) =  Acell(ig) + (zb-zu)*ri;
	  RAcell(ig) = RAcell(ig) + (zb-zu)*ri22;
	  ARcell(ig) = ARcell(ig) + (zb-zu)*logri;
	  dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ri*dzbdrbbbs0;
	  dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ri*dzbdzbbbs0;
	  dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ri*dzbdrbbbs1;
	  dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ri*dzbdzbbbs1;
	  dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ri22*dzbdrbbbs0;
	  dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ri22*dzbdzbbbs0;
	  dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ri22*dzbdrbbbs1;
	  dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ri22*dzbdzbbbs1;
	  dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logri*dzbdrbbbs0;
	  dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logri*dzbdzbbbs0;
	  dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logri*dzbdrbbbs1;
	  dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logri*dzbdzbbbs1;
	elseif fb ==  1
	   Acell(ig) =  Acell(ig) - dz*ri;
	  RAcell(ig) = RAcell(ig) - dz*ri22;
	  ARcell(ig) = ARcell(ig) - dz*logri;
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
          RAcell(ig) = RAcell(ig) + (zb-zd)*ro22;
          ARcell(ig) = ARcell(ig) + (zb-zd)*logro;
	  dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ro*dzbdrbbbs0;
	  dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ro*dzbdzbbbs0;
	  dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ro*dzbdrbbbs1;
	  dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ro*dzbdzbbbs1;
	  dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ro22*dzbdrbbbs0;
	  dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ro22*dzbdzbbbs0;
	  dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ro22*dzbdrbbbs1;
	  dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ro22*dzbdzbbbs1;
	  dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logro*dzbdrbbbs0;
	  dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logro*dzbdzbbbs0;
	  dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logro*dzbdrbbbs1;
	  dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logro*dzbdzbbbs1;
	else
           Acell(ig) =  Acell(ig) + dz*ro;
          RAcell(ig) = RAcell(ig) + dz*ro22;
          ARcell(ig) = ARcell(ig) + dz*logro;
	  if fb == nz % boundary came in from inner edge
             Acell(ig) =  Acell(ig) + (zb-zu)*ri;
            RAcell(ig) = RAcell(ig) + (zb-zu)*ri22;
            ARcell(ig) = ARcell(ig) + (zb-zu)*logri;
	    dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ri*dzbdrbbbs0;
	    dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ri*dzbdzbbbs0;
	    dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ri*dzbdrbbbs1;
	    dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ri*dzbdzbbbs1;
	    dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ri22*dzbdrbbbs0;
	    dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ri22*dzbdzbbbs0;
	    dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ri22*dzbdrbbbs1;
	    dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ri22*dzbdzbbbs1;
	    dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logri*dzbdrbbbs0;
	    dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logri*dzbdzbbbs0;
	    dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logri*dzbdrbbbs1;
	    dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logri*dzbdzbbbs1;
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
          RAcell(ig) = RAcell(ig) + (zb-zu)*ri22;
          ARcell(ig) = ARcell(ig) + (zb-zu)*logri;
	  dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ri*dzbdrbbbs0;
	  dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ri*dzbdzbbbs0;
	  dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ri*dzbdrbbbs1;
	  dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ri*dzbdzbbbs1;
	  dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ri22*dzbdrbbbs0;
	  dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ri22*dzbdzbbbs0;
	  dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ri22*dzbdrbbbs1;
	  dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ri22*dzbdzbbbs1;
	  dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logri*dzbdrbbbs0;
	  dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logri*dzbdzbbbs0;
	  dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logri*dzbdrbbbs1;
	  dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logri*dzbdzbbbs1;
	else
           Acell(ig) =  Acell(ig) - dz*ri;
          RAcell(ig) = RAcell(ig) - dz*ri22;
          ARcell(ig) = ARcell(ig) - dz*logri;
	  if fb == -nz % boundary came in through outer edge of the cell
             Acell(ig) =  Acell(ig) + (zb-zd)*ro;
            RAcell(ig) = RAcell(ig) + (zb-zd)*ro22;
            ARcell(ig) = ARcell(ig) + (zb-zd)*logro;
	    dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ro*dzbdrbbbs0;
	    dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ro*dzbdzbbbs0;
	    dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ro*dzbdrbbbs1;
	    dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ro*dzbdzbbbs1;
	    dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ro22*dzbdrbbbs0;
	    dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ro22*dzbdzbbbs0;
	    dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ro22*dzbdrbbbs1;
	    dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ro22*dzbdzbbbs1;
	    dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logro*dzbdrbbbs0;
	    dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logro*dzbdzbbbs0;
	    dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logro*dzbdrbbbs1;
	    dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logro*dzbdzbbbs1;
	  end
	end
      end
    end    
  end
  k = iedge(j); % k is index to bbbs point
  m = iedge(j+1); % m is index to bbbs point for next edge
  ncell(ig) = ncell(ig)+1;
  rcell(ig,ncell(ig)) = redge(j);
  zcell(ig,ncell(ig)) = zedge(j);
  fcell(ig,ncell(ig)) = fedge(j);
  xa = xedge(j);
  dxadrbbbs0 = dxedgedrbbbs0(j);
  dxadzbbbs0 = dxedgedzbbbs0(j);
  dxadrbbbs1 = dxedgedrbbbs1(j);
  dxadzbbbs1 = dxedgedzbbbs1(j);
  ibc(ig,ncell(ig)) = k;
  dzcelldrbbbs0(ig,ncell(ig)) = DY(k)*dxadrbbbs0;
  dzcelldzbbbs0(ig,ncell(ig)) = DY(k)*dxadzbbbs0+1-xedge(j);
  dzcelldrbbbs1(ig,ncell(ig)) = DY(k)*dxadrbbbs1;
  dzcelldzbbbs1(ig,ncell(ig)) = DY(k)*dxadzbbbs1+xedge(j);
  keep_going = true;
  while keep_going % scan iedge(j):iedge(j+1)
    if k == m
      xb = xedge(j+1);
      dxbdrbbbs0 = dxedgedrbbbs0(j+1);
      dxbdzbbbs0 = dxedgedzbbbs0(j+1);
      dxbdrbbbs1 = dxedgedrbbbs1(j+1);
      dxbdzbbbs1 = dxedgedzbbbs1(j+1);
    else
      xb = 1;
      dxbdrbbbs0 = 0;
      dxbdzbbbs0 = 0;
      dxbdrbbbs1 = 0;
      dxbdzbbbs1 = 0;
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
    
    duma = xx(1)+xx(2)*xa;
    dumb = xx(1)+xx(2)*xb;
    
    dAcelldrbbbs(ig,k  ) = dAcelldrbbbs(ig,k  ) + ...
      DY(k)*((xb-xa)-(yb-ya)/2) - ...
      duma*dxadrbbbs0;
      
    dAcelldrbbbs(ig,k+1) = dAcelldrbbbs(ig,k+1) + ...
      DY(k)*(yb-ya)/2 - ...
      duma*dxadrbbbs1;
        
    dAcelldzbbbs(ig,k  ) = dAcelldzbbbs(ig,k  ) - ...
      rbbbs(k)*(xb-xa)-DX(k)*(yb-ya)/2 - ...
      duma*dxadzbbbs0;      
      
    dAcelldzbbbs(ig,k+1) = dAcelldzbbbs(ig,k+1) + ...
      rbbbs(k)*(xb-xa)+DX(k)*(yb-ya)/2 - ...
      duma*dxadzbbbs1;
    
    dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + dumb*dxbdrbbbs0;
    
    dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + dumb*dxbdzbbbs0;
    
    dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + dumb*dxbdrbbbs1;
    
    dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + dumb*dxbdzbbbs1;
    
    % The coefficients of r(x)*r(x)*dz(x):
    x1(3) = DX(k)*xx(2); % DX(k)^2*DY(k)
    x1(2) = DX(k)*xx(1)+rbbbs(k)*xx(2); % 2*DX(k)*rbbbs(k)*DY(k)
    x1(1) = rbbbs(k)*xx(1); % rbbbs(k)^2*DY(k)
    ya = 1;
    yb = 1;
    dRA = 0;
    for i = 1:3
      ya = ya*xa;
      yb = yb*xb;
      dRA = dRA + x1(i)*(yb-ya)/i/2;
    end
    RAcell(ig) = RAcell(ig) + dRA;
    
    duma = (x1(1) + x1(2)*xa + x1(3)*xa^2)/2;
    dumb = (x1(1) + x1(2)*xb + x1(3)*xb^2)/2;
    
    dRAcelldrbbbs(ig,k  ) = dRAcelldrbbbs(ig,k  ) + ...
      rbbbs(k)*DY(k)*(xb-xa) + DY(k)*(DX(k)-rbbbs(k))*(xb^2-xa^2)/2 - DX(k)*DY(k)*(yb-ya)/3 - ...
      duma*dxadrbbbs0;
    
    dRAcelldrbbbs(ig,k+1) = dRAcelldrbbbs(ig,k+1) + ...
      DY(k)*rbbbs(k)*(xb^2-xa^2)/2 + DX(k)*DY(k)*(yb-ya)/3 - ...
      duma*dxadrbbbs1;
    
    dRAcelldzbbbs(ig,k  ) = dRAcelldzbbbs(ig,k  ) - ...
      rbbbs(k)^2/2*(xb-xa) - DX(k)*rbbbs(k)*(xb^2-xa^2)/2 - DX(k)^2/2*(yb-ya)/3 - ...
      duma*dxadzbbbs0;
    
    dRAcelldzbbbs(ig,k+1) = dRAcelldzbbbs(ig,k+1) + ...
      rbbbs(k)^2/2*(xb-xa) + DX(k)*rbbbs(k)*(xb^2-xa^2)/2 + DX(k)^2/2*(yb-ya)/3 - ...
      duma*dxadzbbbs1;
      
    dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + dumb*dxbdrbbbs0;

    dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + dumb*dxbdrbbbs1;

    dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + dumb*dxbdzbbbs0;

    dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + dumb*dxbdzbbbs1;
   
    
    ra = rbbbs(k)+DX(k)*xa;
    rb = rbbbs(k)+DX(k)*xb; % d(rb*log(rb))/dxb = DX(k)*(log(rb)+1)

    logra = log(ra);
    logrb = log(rb);
    duma = DY(k)*logra;
    dumb = DY(k)*logrb;

    % Changes due to changes of xb
    dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + dumb*dxbdrbbbs0;
    dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + dumb*dxbdrbbbs1;
    dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + dumb*dxbdzbbbs0;
    dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + dumb*dxbdzbbbs1;
    if DX(k) ~= 0
      dAR = DY(k)*(rb*logrb-ra*logra-DX(k)*(xb-xa))/DX(k);
      
      
      % Changes due to changes of rbbbs, zbbbs, DX, DY, xa
      
      dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) + ...
        DY(k)*((1-xb)*logrb-(1-xa)*logra)/DX(k)+dAR/DX(k) - duma*dxadrbbbs0;
           
      dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1) + ...
        DY(k)*(xb*logrb-xa*logra)/DX(k)-dAR/DX(k)         - duma*dxadrbbbs1;
      
      dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - ...
        (rb*logrb-ra*logra-DX(k)*(xb-xa))/DX(k)           - duma*dxadzbbbs0;
      
      dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) + ...
        (rb*logrb-ra*logra-DX(k)*(xb-xa))/DX(k)           - duma*dxadzbbbs1;
    
    else
      dAR = DY(k)*log(rbbbs(k))*(xb-xa);
      
      
      % Changes due to changes of rbbbs, zbbbs, DX, DY
      
      dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) + DY(k)/rbbbs(k)*(xb-xa) - duma*dxadrbbbs0;
    
      dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1)                          - duma*dxadrbbbs1;
      
      dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - log(rbbbs(k))*(xb-xa)  - duma*dxadzbbbs0;
      
      dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) + log(rbbbs(k))*(xb-xa)  - duma*dxadzbbbs1;
    end
    
    ARcell(ig) = ARcell(ig) + dAR;
    
    
    xa = 0; % If the k loop continues with next bbbs point then xa will be 0
    
    dxadrbbbs0 = 0;
    dxadzbbbs0 = 0;
    dxadrbbbs1 = 0;
    dxadzbbbs1 = 0;
    if k == m
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
  ibc(ig,ncell(ig)) = k;
  dzcelldrbbbs0(ig,ncell(ig)) = DY(k)*dxbdrbbbs0;
  dzcelldzbbbs0(ig,ncell(ig)) = DY(k)*dxbdzbbbs0+1-xedge(j+1);
  dzcelldrbbbs1(ig,ncell(ig)) = DY(k)*dxbdrbbbs1;
  dzcelldzbbbs1(ig,ncell(ig)) = DY(k)*dxbdzbbbs1+xedge(j+1);
end


% Walk ccw along the edges of the cells from the exit point (a) to the entry (b)
ig = 0;
for ir = 1:nr
  ri = rg(ir)-dr/2;
  ro = rg(ir)+dr/2;
  ri22 = ri^2/2;
  ro22 = ro^2/2;
  logri = log(ri);
  logro = log(ro);
  for iz = 1:nz
    ig = ig+1;
    if ncell(ig) > 0
      zd = zg(iz)-dz/2;
      zu = zg(iz)+dz/2;
      za = zcell(ig,ncell(ig));
      zb = zcell(ig,1);
      fa = fcell(ig,ncell(ig));
      fb = fcell(ig,1);
      k = ibc(ig,ncell(ig));
      m = ibc(ig,1);
      dzadrbbbs0 = dzcelldrbbbs0(ig,ncell(ig));
      dzadzbbbs0 = dzcelldzbbbs0(ig,ncell(ig));
      dzadrbbbs1 = dzcelldrbbbs1(ig,ncell(ig));
      dzadzbbbs1 = dzcelldzbbbs1(ig,ncell(ig));
      dzbdrbbbs0 = dzcelldrbbbs0(ig,1);
      dzbdzbbbs0 = dzcelldzbbbs0(ig,1);
      dzbdrbbbs1 = dzcelldrbbbs1(ig,1);
      dzbdzbbbs1 = dzcelldzbbbs1(ig,1);
      if fa == -nz % boundary went out through inner edge of the cell
	if fb == nz % boundary came in through inner edge of the cell
           Acell(ig) =  Acell(ig) + (zb-za)*ri;
          RAcell(ig) = RAcell(ig) + (zb-za)*ri22;
          ARcell(ig) = ARcell(ig) + (zb-za)*logri;
	  dAcelldrbbbs(ig,k  ) = dAcelldrbbbs(ig,k  ) - ri*dzadrbbbs0;
	  dAcelldzbbbs(ig,k  ) = dAcelldzbbbs(ig,k  ) - ri*dzadzbbbs0;
	  dAcelldrbbbs(ig,k+1) = dAcelldrbbbs(ig,k+1) - ri*dzadrbbbs1;
	  dAcelldzbbbs(ig,k+1) = dAcelldzbbbs(ig,k+1) - ri*dzadzbbbs1;
	  dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ri*dzbdrbbbs0;
	  dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ri*dzbdzbbbs0;
	  dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ri*dzbdrbbbs1;
	  dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ri*dzbdzbbbs1;
	  dRAcelldrbbbs(ig,k  ) = dRAcelldrbbbs(ig,k  ) - ri22*dzadrbbbs0;
	  dRAcelldzbbbs(ig,k  ) = dRAcelldzbbbs(ig,k  ) - ri22*dzadzbbbs0;
	  dRAcelldrbbbs(ig,k+1) = dRAcelldrbbbs(ig,k+1) - ri22*dzadrbbbs1;
	  dRAcelldzbbbs(ig,k+1) = dRAcelldzbbbs(ig,k+1) - ri22*dzadzbbbs1;
	  dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ri22*dzbdrbbbs0;
	  dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ri22*dzbdzbbbs0;
	  dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ri22*dzbdrbbbs1;
	  dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ri22*dzbdzbbbs1;
	  dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) - logri*dzadrbbbs0;
	  dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - logri*dzadzbbbs0;
	  dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1) - logri*dzadrbbbs1;
	  dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) - logri*dzadzbbbs1;
	  dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logri*dzbdrbbbs0;
	  dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logri*dzbdzbbbs0;
	  dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logri*dzbdrbbbs1;
	  dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logri*dzbdzbbbs1;
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
          RAcell(ig) = RAcell(ig) + (zd-za)*ri22;
          ARcell(ig) = ARcell(ig) + (zd-za)*logri;
	  dAcelldrbbbs(ig,k  ) = dAcelldrbbbs(ig,k  ) - ri*dzadrbbbs0;
	  dAcelldzbbbs(ig,k  ) = dAcelldzbbbs(ig,k  ) - ri*dzadzbbbs0;
	  dAcelldrbbbs(ig,k+1) = dAcelldrbbbs(ig,k+1) - ri*dzadrbbbs1;
	  dAcelldzbbbs(ig,k+1) = dAcelldzbbbs(ig,k+1) - ri*dzadzbbbs1;
	  dRAcelldrbbbs(ig,k  ) = dRAcelldrbbbs(ig,k  ) - ri22*dzadrbbbs0;
	  dRAcelldzbbbs(ig,k  ) = dRAcelldzbbbs(ig,k  ) - ri22*dzadzbbbs0;
	  dRAcelldrbbbs(ig,k+1) = dRAcelldrbbbs(ig,k+1) - ri22*dzadrbbbs1;
	  dRAcelldzbbbs(ig,k+1) = dRAcelldzbbbs(ig,k+1) - ri22*dzadzbbbs1;
	  dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) - logri*dzadrbbbs0;
	  dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - logri*dzadzbbbs0;
	  dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1) - logri*dzadrbbbs1;
	  dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) - logri*dzadzbbbs1;
	  if fb == -nz % boundary came in from outer edge
	     Acell(ig) =  Acell(ig) + (zb-zd)*ro;
	    RAcell(ig) = RAcell(ig) + (zb-zd)*ro22;
	    ARcell(ig) = ARcell(ig) + (zb-zd)*logro;
	    dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ro*dzbdrbbbs0;
	    dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ro*dzbdzbbbs0;
	    dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ro*dzbdrbbbs1;
	    dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ro*dzbdzbbbs1;
	    dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ro22*dzbdrbbbs0;
	    dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ro22*dzbdzbbbs0;
	    dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ro22*dzbdrbbbs1;
	    dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ro22*dzbdzbbbs1;
	    dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logro*dzbdrbbbs0;
	    dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logro*dzbdzbbbs0;
	    dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logro*dzbdrbbbs1;
	    dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logro*dzbdzbbbs1;
	  elseif fb == -1 % boundary came in from top
	     Acell(ig) =  Acell(ig) + dz*ro;
	    RAcell(ig) = RAcell(ig) + dz*ro22;
	    ARcell(ig) = ARcell(ig) + dz*logro;
	  end
	end
      elseif fa == nz % boundary went out through outer edge of the cell
	if fb == -nz % boundary came in through outer edge of the cell
           Acell(ig) =  Acell(ig) + (zb-za)*ro;
          RAcell(ig) = RAcell(ig) + (zb-za)*ro22;
          ARcell(ig) = ARcell(ig) + (zb-za)*logro;
	  dAcelldrbbbs(ig,k  ) = dAcelldrbbbs(ig,k  ) - ro*dzadrbbbs0;
	  dAcelldzbbbs(ig,k  ) = dAcelldzbbbs(ig,k  ) - ro*dzadzbbbs0;
	  dAcelldrbbbs(ig,k+1) = dAcelldrbbbs(ig,k+1) - ro*dzadrbbbs1;
	  dAcelldzbbbs(ig,k+1) = dAcelldzbbbs(ig,k+1) - ro*dzadzbbbs1;
	  dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ro*dzbdrbbbs0;
	  dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ro*dzbdzbbbs0;
	  dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ro*dzbdrbbbs1;
	  dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ro*dzbdzbbbs1;
	  dRAcelldrbbbs(ig,k  ) = dRAcelldrbbbs(ig,k  ) - ro22*dzadrbbbs0;
	  dRAcelldzbbbs(ig,k  ) = dRAcelldzbbbs(ig,k  ) - ro22*dzadzbbbs0;
	  dRAcelldrbbbs(ig,k+1) = dRAcelldrbbbs(ig,k+1) - ro22*dzadrbbbs1;
	  dRAcelldzbbbs(ig,k+1) = dRAcelldzbbbs(ig,k+1) - ro22*dzadzbbbs1;
	  dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ro22*dzbdrbbbs0;
	  dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ro22*dzbdzbbbs0;
	  dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ro22*dzbdrbbbs1;
	  dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ro22*dzbdzbbbs1;
	  dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) - logro*dzadrbbbs0;
	  dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - logro*dzadzbbbs0;
	  dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1) - logro*dzadrbbbs1;
	  dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) - logro*dzadzbbbs1;
	  dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logro*dzbdrbbbs0;
	  dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logro*dzbdzbbbs0;
	  dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logro*dzbdrbbbs1;
	  dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logro*dzbdzbbbs1;
	else % Upper outer corner must be inside the plasma
          if ncell(ig+1) == 0
	    Acell(ig+1) = Ag;
	  end
          if ncell(ig+1+nz) == 0
	    Acell(ig+1+nz) = Ag;
	  end
           Acell(ig) =  Acell(ig) + (zu-za)*ro;
          RAcell(ig) = RAcell(ig) + (zu-za)*ro22;
          ARcell(ig) = ARcell(ig) + (zu-za)*logro;
	  dAcelldrbbbs(ig,k  ) = dAcelldrbbbs(ig,k  ) - ro*dzadrbbbs0;
	  dAcelldzbbbs(ig,k  ) = dAcelldzbbbs(ig,k  ) - ro*dzadzbbbs0;
	  dAcelldrbbbs(ig,k+1) = dAcelldrbbbs(ig,k+1) - ro*dzadrbbbs1;
	  dAcelldzbbbs(ig,k+1) = dAcelldzbbbs(ig,k+1) - ro*dzadzbbbs1;
	  dRAcelldrbbbs(ig,k  ) = dRAcelldrbbbs(ig,k  ) - ro22*dzadrbbbs0;
	  dRAcelldzbbbs(ig,k  ) = dRAcelldzbbbs(ig,k  ) - ro22*dzadzbbbs0;
	  dRAcelldrbbbs(ig,k+1) = dRAcelldrbbbs(ig,k+1) - ro22*dzadrbbbs1;
	  dRAcelldzbbbs(ig,k+1) = dRAcelldzbbbs(ig,k+1) - ro22*dzadzbbbs1;
	  dARcelldrbbbs(ig,k  ) = dARcelldrbbbs(ig,k  ) - logro*dzadrbbbs0;
	  dARcelldzbbbs(ig,k  ) = dARcelldzbbbs(ig,k  ) - logro*dzadzbbbs0;
	  dARcelldrbbbs(ig,k+1) = dARcelldrbbbs(ig,k+1) - logro*dzadrbbbs1;
	  dARcelldzbbbs(ig,k+1) = dARcelldzbbbs(ig,k+1) - logro*dzadzbbbs1;
	  if fb == nz % boundary came in from inner edge
	     Acell(ig) =  Acell(ig) + (zb-zu)*ri;
	    RAcell(ig) = RAcell(ig) + (zb-zu)*ri22;
	    ARcell(ig) = ARcell(ig) + (zb-zu)*logri;
	    dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ri*dzbdrbbbs0;
	    dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ri*dzbdzbbbs0;
	    dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ri*dzbdrbbbs1;
	    dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ri*dzbdzbbbs1;
	    dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ri22*dzbdrbbbs0;
	    dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ri22*dzbdzbbbs0;
	    dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ri22*dzbdrbbbs1;
	    dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ri22*dzbdzbbbs1;
	    dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logri*dzbdrbbbs0;
	    dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logri*dzbdzbbbs0;
	    dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logri*dzbdrbbbs1;
	    dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logri*dzbdzbbbs1;
	  elseif fb ==  1
	     Acell(ig) =  Acell(ig) - dz*ri;
	    RAcell(ig) = RAcell(ig) - dz*ri22;
	    ARcell(ig) = ARcell(ig) - dz*logri;
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
            RAcell(ig) = RAcell(ig) + (zb-zd)*ro22;
            ARcell(ig) = ARcell(ig) + (zb-zd)*logro;
	    dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ro*dzbdrbbbs0;
	    dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ro*dzbdzbbbs0;
	    dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ro*dzbdrbbbs1;
	    dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ro*dzbdzbbbs1;
	    dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ro22*dzbdrbbbs0;
	    dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ro22*dzbdzbbbs0;
	    dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ro22*dzbdrbbbs1;
	    dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ro22*dzbdzbbbs1;
	    dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logro*dzbdrbbbs0;
	    dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logro*dzbdzbbbs0;
	    dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logro*dzbdrbbbs1;
	    dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logro*dzbdzbbbs1;
	  else
             Acell(ig) =  Acell(ig) + dz*ro;
            RAcell(ig) = RAcell(ig) + dz*ro22;
            ARcell(ig) = ARcell(ig) + dz*logro;
	    if fb == nz % boundary came in from inner edge
               Acell(ig) =  Acell(ig) + (zb-zu)*ri;
              RAcell(ig) = RAcell(ig) + (zb-zu)*ri22;
              ARcell(ig) = ARcell(ig) + (zb-zu)*logri;
	      dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ri*dzbdrbbbs0;
	      dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ri*dzbdzbbbs0;
	      dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ri*dzbdrbbbs1;
	      dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ri*dzbdzbbbs1;
	      dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ri22*dzbdrbbbs0;
	      dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ri22*dzbdzbbbs0;
	      dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ri22*dzbdrbbbs1;
	      dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ri22*dzbdzbbbs1;
	      dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logri*dzbdrbbbs0;
	      dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logri*dzbdzbbbs0;
	      dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logri*dzbdrbbbs1;
	      dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logri*dzbdzbbbs1;
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
            RAcell(ig) = RAcell(ig) + (zb-zu)*ri22;
            ARcell(ig) = ARcell(ig) + (zb-zu)*logri;
	    dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ri*dzbdrbbbs0;
	    dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ri*dzbdzbbbs0;
	    dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ri*dzbdrbbbs1;
	    dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ri*dzbdzbbbs1;
	    dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ri22*dzbdrbbbs0;
	    dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ri22*dzbdzbbbs0;
	    dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ri22*dzbdrbbbs1;
	    dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ri22*dzbdzbbbs1;
	    dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logri*dzbdrbbbs0;
	    dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logri*dzbdzbbbs0;
	    dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logri*dzbdrbbbs1;
	    dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logri*dzbdzbbbs1;
	  else
             Acell(ig) =  Acell(ig) - dz*ri;
            RAcell(ig) = RAcell(ig) - dz*ri22;
            ARcell(ig) = ARcell(ig) - dz*logri;
	    if fb == -nz % boundary came in through outer edge of the cell
               Acell(ig) =  Acell(ig) + (zb-zd)*ro;
              RAcell(ig) = RAcell(ig) + (zb-zd)*ro22;
              ARcell(ig) = ARcell(ig) + (zb-zd)*logro;
	      dAcelldrbbbs(ig,m  ) = dAcelldrbbbs(ig,m  ) + ro*dzbdrbbbs0;
	      dAcelldzbbbs(ig,m  ) = dAcelldzbbbs(ig,m  ) + ro*dzbdzbbbs0;
	      dAcelldrbbbs(ig,m+1) = dAcelldrbbbs(ig,m+1) + ro*dzbdrbbbs1;
	      dAcelldzbbbs(ig,m+1) = dAcelldzbbbs(ig,m+1) + ro*dzbdzbbbs1;
	      dRAcelldrbbbs(ig,m  ) = dRAcelldrbbbs(ig,m  ) + ro22*dzbdrbbbs0;
	      dRAcelldzbbbs(ig,m  ) = dRAcelldzbbbs(ig,m  ) + ro22*dzbdzbbbs0;
	      dRAcelldrbbbs(ig,m+1) = dRAcelldrbbbs(ig,m+1) + ro22*dzbdrbbbs1;
	      dRAcelldzbbbs(ig,m+1) = dRAcelldzbbbs(ig,m+1) + ro22*dzbdzbbbs1;
	      dARcelldrbbbs(ig,m  ) = dARcelldrbbbs(ig,m  ) + logro*dzbdrbbbs0;
	      dARcelldzbbbs(ig,m  ) = dARcelldzbbbs(ig,m  ) + logro*dzbdzbbbs0;
	      dARcelldrbbbs(ig,m+1) = dARcelldrbbbs(ig,m+1) + logro*dzbdrbbbs1;
	      dARcelldzbbbs(ig,m+1) = dARcelldzbbbs(ig,m+1) + logro*dzbdzbbbs1;
	    end
	  end
	end
      end    
    end
  end
end

% Find the rest of the cells that must be inside the plasma
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

