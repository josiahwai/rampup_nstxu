function [xb, rootin01, x01min, x01max] = gs_solve_hermite_cubic(y, yb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   [xb, rootin01, x01min, x01max] = gs_solve_hermite_cubic(y, yb)
%
%  PURPOSE: Find where the function y = yb, with y given at x = -1, 0, 1, 2,
%           and interpolated between points using cubic Hermite splines
%
%  INPUTS: y,  values at x = -1, 0, 1, 2
%          yb, value for which the roots xb are sought
%
%  OUTPUTS:  xb, either 1 or 3 roots are real, first 1 always real
%            rootin01, true if at least one root is between 0 and 1
%            x01min, smallest root in the interval 0 to 1, if existent
%            x01max, biggest root in the interval 0 to 1, if existent
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	2/14/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = ( -y(1) + 3*y(2) - 3*y(3)+y(4))/2;
b = (2*y(1) - 5*y(2) + 4*y(3)-y(4))/2;
c = ( -y(1)        +     y(3)     )/2;
d =             y(2)              -yb;

%xb = roots(y(:)'*[-1 2 -1 0; 3 -5 0 2; -3 4 1 0; 1 -1 0 0]/2-[0 0 0 yb]);
xb = roots([a b c d]);

% if a==0 we only get two roots and if also b==0 we get one
if length(xb) == 0
  xb = [-1;-1;-1];
elseif length(xb) == 1
  xb = [xb; -1; -1];
elseif length(xb) == 2
  xb(3) = -1;
end

ii = imag(xb) == 0 & xb >= 0 & xb < 1;

rootin01 = any(ii);
if rootin01
  x01min = min(xb(ii));
  x01max = max(xb(ii));
else
  x01min = 1;
  x01max = 0;
end
