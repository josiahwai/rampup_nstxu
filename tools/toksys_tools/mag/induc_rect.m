  function ind= induc_rect(a,b,r)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAXX:
%         ind= induc_rect(a,b,r)
%
%  PURPOSE:  Calculate one-turn, self inductance of plane rectangle
%            made from a circular cross section wire of radius r. 
%            The rectangle dimensions are a by b
%
%  INPUT:
%	a,b = dimensions of rectangle, (vector) [m]
%       r   = radius or wire cross (vector) [m]
%
%  OUTPUT:
%	ind = one turn, self inductance of coil [Henries]
%             Multiply by N^2 if coil has N turns
%
%  RESTRICTIONS:
%        All input vectors must be the same length (or a scalar)

%  METHOD:  
%
%  WRITTEN BY: Jim Leuer 8/25/99  
%  Based on approximate algorithm in AIP handbook.
%  (Grover has similar but slightly different form using sinh^-1 ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    d= sqrt(a.^2+b.^2); % diagonal
    ind= 4e-7*(a.*log(2*a.*b./(r.*(a+d))) + b.*log(2*a.*b./(r.*(b+d))) ...
              +2*d - 1.75*(a+b));
    return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test stuff below

% All looks good jal 8-25-99

% ---------------------
%TEST 1 for SQUARE
  a=1; b=1; r=.1;

  ind= induc_rect(a,b,r)

% Grover for square:

  ind_square= 8e-7*a*(log(a./r) -0.77401 + 1/4)

% Grover for Rectangle using asihn terms:

    d= sqrt(a.^2+b.^2); % diagonal
    ind_grov= 4e-7*(a.*log(2*a./r) + b.*log(2*b./r) ...
              -a.*asinh(a./b) -b.*asinh(b./a)...
              +2*d - 1.75*(a+b))
% ---------------------
%TEST 2 for 2/1 rectangle 

  a=1; b=2; r=.1;

  ind= induc_rect(a,b,r)

% Grover table 8: for b/a = beta = 0.5 ==>
  alpha= 2.962;

  l= 2*(a+b); % perimeter

  ind_grov= 2e-7*l*(log(2*l/r) - alpha +1/4)

% ---------------------
% TEST 3  2-Vectors & 1 scalar

  a= ones(2,1); b= (1:2)'; r= .1;

  ind= induc_rect(a,b,r)
 
