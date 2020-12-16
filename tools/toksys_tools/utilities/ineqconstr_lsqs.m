function x = ineqconstr_lsqs(A,b,B,d)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  x = ineqconstr_lsqs(A,b,B,d)
%
%  PURPOSE:  Inequality constrained least squares solution.  Solves
%		min || Ax - b || subject to ||Bx - d|| <= a.
%
%	**************NOT WORKING YET**************
%
%  INPUT:
%	A = mxn matrix in LHS of Ax=b least squares problem
%	b = mx1 vector in RHS of Ax=b least squares problem
%	B = pxn matrix in constraint equation Bx=d
%	d = px1 vector (all elts >=0) in constraint equation Bx=d
%
%  OUTPUT:
%	x = solution vector for constrained least squares problem
%
%  RESTRICTIONS:  A and B must be full rank.
 
%  METHOD:  From Golub & Van Loan, section 12.1.1.
%
%  WRITTEN BY:  Mike Walker 	ON 	11/25/97
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[D1,D2,R,U,V,Q]=gsvd(A,B);

% Need to figure out how to convert matrices produced by LAPACK into form
% needed for Golub and Van Loan algorithm.

