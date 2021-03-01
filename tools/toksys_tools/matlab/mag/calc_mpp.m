function [mpp,gbr2p,gbz2p] = calc_mpp(rgg,zgg,dr,dz,nz,nr)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  [mpp,gbr2p,gbz2p] = calc_mpp(rgg,zgg,dr,dz,nz,nr)
%
%  PURPOSE:  calculate mutuals between all plasma gridpoints
%
%  INPUT:
%	rgg = r coords of all grids
%	zgg = z coords of all grids
%	dr  = r width of all grid "current elements"
%	dz  = z width of all grid "current elements"
%	nz  = number of grids in z direction
%	nr  = number of grids in r direction
%
%  OUTPUT:
%	mpp = mutuals between all plasma gridpoints (uH)
%	gbr2p = Br greens functions from plasma grid to itself (T/MA)
%	gbz2p = Bz greens functions from plasma grid to itself (T/MA)
%
% !!NOTE!!: This data is stored in a compressed format using the fact that it
%	is only the difference in z coordinates which enters the mutual and
%	Greens' function calculations:  All output matrices are of
%	size (nz*nr) by nr, where columns represent different values of 
%	radial coordinate.  To extract a mutual inductance or Green's function
%	value for the kth grid point (measurement point) versus the jth grid 
%	point (current source), do the following:
%	(1) calculate the number of grids n that grid point k is above the 
%	    jth grid point in z direction (negative value n if kth is below jth)
%	(2) calculate the column number m (m=1,...,nr) that grid point k is in.
%	(3) mutual inductance = mpp((m-1)*nz+|n|+1,j)
%	(4) Br greens function = sign(n) * gbr2p((m-1)*nz+|n|+1,j)
%	    where sign(n) =  1, if n>=0
%		          = -1, if n<0.
%	(5) Bz greens function = gbz2p((m-1)*nz+|n|+1,j)
%	Note that it is only important to distinguish between the source and
%	measurement point for Br Green's function values.
%	This "uncompress" procedure is implemented in get_plasma_greens.m.
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	11/26/97
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rtmp = rgg(:);
ztmp = zgg(:);
s1 = size(rtmp,1);
onev = ones(s1,1);
rg = rgg(1,:)';
zg = zgg(:,1);
mpp = zeros(s1,length(rg));
gbr2p = mpp;
gbz2p = mpp;

for k=1:length(rg)
   [mut,gbr,gbz] = mindbf_gen(rg(k)*onev,rtmp,zg(1)*onev,ztmp,dr,dz);
   mpp(:,k) = mut;
   gbr2p(:,k) = gbr;
   gbz2p(:,k) = gbz;
end

mpp   = mpp * 1e6;		% H -> uH
gbr2p = gbr2p * 1e6;		% currents -> MA
gbz2p = gbz2p * 1e6;		% currents -> MA

%save mpp mpp gbr2p gbz2p
