function det = det_poly(matrix)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  det = det_poly(matrix)
%
%  PURPOSE:  Compute determinant of polynomial matrix.  A poly matrix is
%	of the form matrix(P,M,M), which represents an MxM matrix of polynomials
%	of at most degree P-1. The P elements of matrix(:,k,j) represent a
%	poly in matlab representation.
%
%  INPUT:
%	matrix = matrix of polynomials
%
%  OUTPUT:
%	det = polynomial that is determinant of matrix
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	2/13/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P,M1,M] = size(matrix);
if(M1~=M)
   wait('ERROR det_poly: last 2 dimensions must be equal')
   return;
end

% Uses recursive calls to itself, expanding along the first column during
% each call.

if(M==1)
   det = matrix;
else
   det=0;
   for k=1:M
      det1 = det_poly(submatrix_poly(matrix,[1:k-1 k+1:M],2:M));
      det = det + (-1)^(k-1) * conv(matrix(:,k,1), det1);
   end
end
