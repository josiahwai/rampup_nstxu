function submatrix = submatrix_poly(matrix,rowindex,colindex)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  submatrix = submatrix_poly(matrix,rowindex,colindex)
%
%  PURPOSE:  Extract selected submatrix out of matrix of polynomials.
%
%  INPUT:
%	matrix	 = matrix of polys to extract from
%	rowindex = indices of rows to extract
%	colindex = indices of columns to extract
%
%  OUTPUT:
%	submatrix = extracted sub-matrix
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	2/13/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P,M,N] = size(matrix);
if(max(rowindex)>M)
   wait('ERROR submatrix_poly: selected row outside size of matrix')
   return;
end
if(max(colindex)>N)
   wait('ERROR submatrix_poly: selected column outside size of matrix')
   return;
end

for k=1:length(rowindex)
   for j=1:length(colindex)
      submatrix(:,k,j) = matrix(:,rowindex(k),colindex(j));
   end
end

