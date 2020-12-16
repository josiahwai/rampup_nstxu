function m = slope(x,y)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  m = slope(x,y)
%
%  PURPOSE: Calculate slope of line specified by two points. 
%
%  INPUT:
%	x = length 2 vector giving x values
%	y = length 2 vector giving y values
%
%  OUTPUT:
%	m = slope of line
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	10/24/96
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = (y(2)-y(1))/(x(2)-x(1));

end
