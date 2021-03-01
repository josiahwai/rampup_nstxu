function product = mpp_x_vec(mpp,vec)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  product = mpp_x_vec(mpp,vec)
%
%  PURPOSE:  Compute product of matrix mpp and vector vec, where mpp is
%	     stored in compressed format. Can be used also for gbz2p in
%	     place of mpp. Use gbr2p_x_vec for multiplication with gbr2p.
%
%  INPUT:    mpp, grid to grid mutuals in compressed format
%            vec, vector or matrix with length(vec(:)) = size(mpp,1)
%
%  OUTPUT:   product, product of matrix multiplication, same size as vec
%
%  SEE ALSO: get_plasma_greens, gbr2p_x_vec
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  VERSION @(#)mpp_x_vec.m	1.3 11/06/12
%
%  WRITTEN BY:  Mike Walker 	ON 	6/11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ngrids,nr] = size(mpp);
nz = ngrids/nr;
if length(vec(:))~=ngrids
   disp('ERROR mpp_x_vec: length(vec(:)) must match mpp 1:st dimension')
   product = [];
   return;
end

% !!NOTE!!: This data is stored in a compressed format using the fact that it
%       is only the difference in z coordinates which enters the mutual and
%       Greens' function calculations:  All output matrices are of
%       size (nz*nr) by nr, where columns represent different values of
%       radial coordinate.  To extract a mutual inductance or Green's function
%       value for the kth grid point (measurement point) versus the jth grid
%       point (current source), do the following:
%       (1) calculate the number of grids n that grid point k is above the
%           jth grid point in z direction (negative value n if kth is below jth)
%       (2) calculate the column number m (m=1,...,nr) that grid point k is in.
%       (3) mutual inductance = mpp((m-1)*nz+|n|+1,j)
%       (4) Br greens function = sign(n) * gbr2p((m-1)*nz+|n|+1,j)
%           where sign(n) =  1, if n>=0
%                         = -1, if n<0.
%       (5) Bz greens function = gbz2p((m-1)*nz+|n|+1,j)
%       Note that it is only important to distinguish between the source and
%       measurement point for Br Green's function values.
%       This "uncompress" procedure is implemented in get_plasma_greens.m.

[nzvec,nrvec]=size(vec);

product = zeros(2,1);	% force it to be column vector

for j=1:ngrids		% index of output element of matrix multiply
   z_idx = mod(j-1,nz)+1;	% z index of that output element
   r_idx = ceil(j/nz);		% r index of that output element

   sum = 0;
   for m=1:nr			% r index of input element
      k1 = (r_idx-1)*nz;
      k2 = (m-1)*nz;
      for n1=1:nz		% z index of input element
         n = z_idx-n1;
         k = k2+n1;
         mp(k) = mpp(k1+abs(n)+1,m);
      end
   end

   sum = sum + mp*vec(:);
   product(j) = sum;
end
product = reshape(product,nzvec,nrvec);
