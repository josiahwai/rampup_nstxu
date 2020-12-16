  function mut=mutind(ra,za,rb,zb)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  mut=mutind(ra,za,rb,zb)
%
%  PURPOSE:  Calculate mutual inductances between loops centered at r=0.   
%
%  INPUT:
%	ra,za = vectors of (r,z) coordinates (meters) of first current loops  
%	rb,zb = vectors of (r,z) coordinates (meters) of 2nd current loops    
% Note arguments have different ordering convention than IDL function.
%
%  OUTPUT:
%	mut = vector of mutual inductances between current loops (in 
%	      micro-Henries), mut(k) = mutual from loop at (ra(k),za(k)) to 
%	      loop at (rb(k),zb(k)). 
%
%  RESTRICTIONS:  Coordinates (r,z) must represent concentric current loops.
%   All input vectors must be the same length.

%  METHOD:  
%
%  WRITTEN BY:  Dave Humphreys 	ON 	2/12/95  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    m=((rb+ra).^2 - (rb-ra).^2)./((rb+ra).^2 + (za-zb).^2);
    [kk,ek]=ellipke(m);
    mut=0.8*pi*sqrt(ra.*rb).*((1.-.5*m).*kk-ek)./sqrt(m);
