function [V,d] = eigsort(A)

%  SYNTAX:  eigsort(A), d = eigsort(A), OR [V,d] = eigsort(A) 
%
%  PURPOSE:  Calculate eigenvectors V and vector of eigenvalues d, 
%       sorted by real part.   Works the same as eig.m, except eigenvectors
%       are always returned as a vector (eig returns sometimes as a vector,
%       sometimes as a diagonal matrix) and eigenvalues are sorted.
%
%  INPUT:
%       A = matrix to calculate eigenvalues for
%
%  OUTPUT: (see eig.m)
%       V = eigenvectors of A
%       d = eigenvalues of A, i'th eigenvector corresponds to V(:,i)
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Dave Humphreys  ON      ??
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargout == 0)
  V = eig(A);    %V is now eigenvalue vector...
  [tmp,itmp]=sort(real(V));
  V = flipud(V(itmp)) %set V to eigenvalue vector so get only that
end
if(nargout == 1)
  V = eig(A);
  [tmp,itmp] = sort(real(V));
  V = flipud(V(itmp)); %set V to eigenvalue vector so get only that
end
if(nargout == 2)
  [V,D] = eig(A);
  d = diag(D);
  [tmp,itmp] = sort(real(d));
  d = flipud(d(itmp));
  V = V(:,flipud(itmp));
end

