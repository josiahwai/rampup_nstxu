%  USAGE:   gs_trace_edge
%
%  PURPOSE: Find coordinates where the boundary intersects edges of grid cells
%           (used in calculation of how much of cells are covered by plasma)
%
%  INPUTS: rbbbs, zbbbs, nbbbs: R, Z of boundary and number of points
%          rmaxis, zmaxis, position of axis (used to order points w.r.t theta)
%          Other precomputed variables:	
%            rg, zg, dr, dz (grid variables)
%
%  OUTPUTS: redge, zedge, nedge: R, Z, and number of points that intersect grid cells
%           rhoedge, thedge, distance to magnetic axis and poloidal angle
%           fedge, direction when entering cells going ccw -1=down, 1=up, -nz=in, nz=out
%           xedge, values of x to use in interpolation formula to obtain redge, zedge
%           iedge, indices to the first of two bbbs points (j in the formula)
%           gedge, indices to cells entered at redge, zedge, going along bbbs
%	
%  METHOD:  The linear interpolation formula:
%           R = rbbbs(j) + (rbbbs(j+1)-rbbbs(j))*x
%           Z = zbbbs(j) + (zbbbs(j+1)-zbbbs(j))*x
%           is solved for x-values that make either R or Z
%           intersect a cell edge. All solutions with x between 0 and 1 are
%           returned in redge, zedge. The points are sorted according to thedge.
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	3/12/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DX = diff(rbbbs)';
DY = diff(zbbbs)';
irbbbs = (rbbbs-rg(1))/dr; % rbbbs in grid number units, starting at zero
izbbbs = (zbbbs-zg(1))/dz; % zbbbs in grid number units, starting at zero
ir1 = round(min([irbbbs(1:end-1)'; irbbbs(2:end)']));
ir2 = round(max([irbbbs(1:end-1)'; irbbbs(2:end)']))-1;
iz1 = round(min([izbbbs(1:end-1)'; izbbbs(2:end)']));
iz2 = round(max([izbbbs(1:end-1)'; izbbbs(2:end)']))-1;
nedge = 0;
for j = 1:nbbbs-1
  i1 = nedge+1;
  for k = ir1(j) : ir2(j)
    nedge = nedge+1;
    redge(nedge) = (k+0.5)*dr+rg(1);
    iedge(nedge) = j;
    xedge(nedge) = (redge(nedge)-rbbbs(j))/DX(j);
    zedge(nedge) = zbbbs(j)+xedge(nedge)*DY(j);
    if DX(j) > 0
      fedge(nedge) = nz;
      gedge(nedge) = (k+1)*nz+round((zedge(nedge)-zg(1))/dz)+1;
    else
      fedge(nedge) = -nz;
      gedge(nedge) = k*nz+round((zedge(nedge)-zg(1))/dz)+1;
    end
  end
  for k = iz1(j) : iz2(j)
    nedge = nedge+1;
    zedge(nedge) = (k+0.5)*dz+zg(1);
    iedge(nedge) = j;
    xedge(nedge) = (zedge(nedge)-zbbbs(j))/DY(j);
    redge(nedge) = rbbbs(j)+xedge(nedge)*DX(j);
    if DY(j) > 0
      fedge(nedge) = 1;
      gedge(nedge) = round((redge(nedge)-rg(1))/dr)*nz+k+2;
    else
      fedge(nedge) = -1;
      gedge(nedge) = round((redge(nedge)-rg(1))/dr)*nz+k+1;
    end
  end
  i2 = nedge;
  n = i2-i1+1;
  [xedge(i1:i2), kk(1:n)] = sort(xedge(i1:i2));
  fedge(i1:i2) = fedge(i1-1+kk(1:n));
  gedge(i1:i2) = gedge(i1-1+kk(1:n));
  iedge(i1:i2) = iedge(i1-1+kk(1:n));
  redge(i1:i2) = redge(i1-1+kk(1:n));
  zedge(i1:i2) = zedge(i1-1+kk(1:n));
end

thedge(1:nedge) = angle((redge(1:nedge)-rmaxis)+1i*(zedge(1:nedge)-zmaxis));
rhoedge(1:nedge) = sqrt((redge(1:nedge)-rmaxis).^2+(zedge(1:nedge)-zmaxis).^2);
  
 thedge(nedge+1) =  thedge(1)+2*pi;
  redge(nedge+1) =   redge(1);
  zedge(nedge+1) =   zedge(1);
  iedge(nedge+1) =   iedge(1);
  fedge(nedge+1) =   fedge(1);
  gedge(nedge+1) =   gedge(1);
  xedge(nedge+1) =   xedge(1);
rhoedge(nedge+1) = rhoedge(1);
