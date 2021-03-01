function corners = elt_to_corners(element)
 %
%  SYNTAX:  corners = elt_to_corners(element)
%
%  PURPOSE:  Convert standard "element form" ([Z;R;dZ;dR;AC;AC2]) to set
%	of corners of the element.  (For a description of this standard form,
%	see Toksys users guide.)
%
%  INPUT:
%	element = [Z;R;dZ;dR;AC;AC2]
%
%  OUTPUT:
%	corners = [x y], where x = [x1 x2 x3 x4]^T, y = [y1 y2 y3 y4]^T
%		are coordinates of the four corners of the parallelogram.
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	3/22/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)elt_to_corners.m	1.1 03/03/15

AC=element(5)*pi/180;
AC2=element(6)*pi/180;

if AC~=0
   Z = element(1);
   R = element(2);
   dZ = element(3);
   dR = element(4);

   l = dR/cos(AC);

   x(1) = R + l/2*cos(AC);
   y(1) = Z + l/2*sin(AC) - dZ/2;

   x(2) = R + l/2*cos(AC);
   y(2) = Z + l/2*sin(AC) + dZ/2;

   x(3) = R - l/2*cos(AC);
   y(3) = Z - l/2*sin(AC) - dZ/2;

   x(4) = R - l/2*cos(AC);
   y(4) = Z - l/2*sin(AC) + dZ/2;

elseif AC2~=0
   Z = element(1);
   R = element(2);
   dZ = element(3);
   dR = element(4);

   l = dZ/sin(AC2);
   
   x(1) = R + l/2*cos(AC2) - dR/2;
   y(1) = Z + l/2*sin(AC2);

   x(2) = R + l/2*cos(AC2) + dR/2;
   y(2) = Z + l/2*sin(AC2);
   
   x(3) = R - l/2*cos(AC2) - dR/2;
   y(3) = Z - l/2*sin(AC2);

   x(4) = R - l/2*cos(AC2) + dR/2;
   y(4) = Z - l/2*sin(AC2);

else
   Z = element(1);
   R = element(2);
   dZ = element(3);
   dR = element(4);

   x(1) = R - dR/2;
   y(1) = Z + dZ/2;

   x(2) = R + dR/2;
   y(2) = Z + dZ/2;

   x(3) = R - dR/2;
   y(3) = Z - dZ/2;

   x(4) = R + dR/2;
   y(4) = Z - dZ/2;

end

corners = [x' y'];
