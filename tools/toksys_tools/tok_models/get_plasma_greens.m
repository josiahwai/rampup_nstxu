function [mp,gbr,gbz]=get_plasma_greens(mpp,gbr2p,gbz2p,src_r_idx,src_z_idx,nz)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  [mp,gbr,gbz]=...
%	get_plasma_greens(mpp,gbr2p,gbz2p,src_r_idx,src_z_idx,nz)
%
%  PURPOSE:  Get mutual inductance, Br, and Bz green's function responses 
%	(plasma green's functions) to a given grid point current source.
%
%  INPUT:
%	mpp = plasma grid to plasma grid mutuals in compressed format
%	gbr2p = Br green's functions from plasma grid to plasma grid 
%		in compressed format
%	gbz2p = Bz green's functions from plasma grid to plasma grid 
%		in compressed format
%	src_r_idx = index of r coordinate of current source in overall grid
%			(0 < src_r_idx <= nr)
%	src_z_idx = index of z coordinate of current source in overall grid
%			(0 < src_z_idx <= nz)
%	nz    = number of grids in z dimension
%
%  OUTPUT:
%	mp  = mutuals from selected grid current source to all other grids
%	gbr = Br greens fns from selected grid current source to all other grids
%	gbz = Bz greens fns from selected grid current source to all other grids
%   (Output units are same as input objects mpp, gbr2p, gbz2p.)
%
%  SEE ALSO: mpp_x_vec
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	2/5/98
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(src_r_idx)~=length(src_z_idx)
   fprintf('get_plasma_greens: r and z indices must be same length\n')
   return
else
   len = length(src_r_idx);
end

ngrids = size(mpp,1);

mp = zeros(ngrids,len);
gbr = zeros(ngrids,len);
gbz = zeros(ngrids,len);
for j=1:len
   for k = 1:ngrids
      n = mod(k-1,nz)+1 - src_z_idx(j);
      m = floor((k-1)/nz)+1;
      mp(k,j) = mpp((m-1)*nz+abs(n)+1,src_r_idx(j));
      gbr(k,j) = sign(n+0.5)*gbr2p((m-1)*nz+abs(n)+1,src_r_idx(j));
      gbz(k,j) = gbz2p((m-1)*nz+abs(n)+1,src_r_idx(j));
   end
end
