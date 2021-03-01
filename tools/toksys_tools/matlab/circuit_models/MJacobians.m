function [J]=MJacobians(ra,za,rb,zb,delr,delz)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  J=MJacobians(ra,za,rb,zb)
%
%  PURPOSE:  Compute Jacobians dm(ra,za;rb,zb)/d(ra,za) and 
%	dm(ra,za;rb,zb)/d(rb,zb), where m(ra,za;rb,zb)
%	is the mutual from (ra,za) to (rb,zb).
%
%  INPUT:
%	ra,za =(r,z) coordinates of current source (meters)
%	rb,zb =(r,z) coordinates of point where psi is measured (meters)
%	delr  = 
%	delz  = 
%  Legal input combinations are:
%	- ra,za scalar, rb,zb column vectors of same length
%	- rb,zb scalar, ra,za column vectors of same length
%	- ra,za,rb,zb all column vectors of same length
%
%  OUTPUT:
%	J = Jacobian (units = Tesla/Amp) matrix, kth row of which holds
%		[dm/d(ra) dm/d(za) dm/d(rb) dm/d(zb)]
%	for the kth elements in vectors ra,za and/or rb,zb.
%
%  RESTRICTIONS:  Note that this function does not handle carefully the
%	case where source and measurement location are very close.  It DOES
%	handle the case where source and measurement location are identical.
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	11/8/96
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    f=((rb+ra).^2 - (rb-ra).^2)./((rb+ra).^2 + (za-zb).^2);
    k=sqrt(f);
    
    [kk,ek]=ellipke(f);

    dkkdk  = ek./(k.*(1.-f))-kk./k;
    dekdk  = (ek-kk)./k;
    dkdza = (zb-za).*sqrt(((rb+ra).^2-(rb-ra).^2)./((rb+ra).^2+(za-zb).^2).^3);
    dkdzb = -dkdza;

    d = (rb+ra).^2 + (za-zb).^2;
    n = (rb+ra).^2 - (rb-ra).^2;
    dddr = 2*(rb + ra);
    dndra = 4*rb;
    dndrb = 4*ra;
    dkdra = (.5./k).*(d.*dndra - n.*dddr)./d.^2;
    dkdrb = (.5./k).*(d.*dndrb - n.*dddr)./d.^2;
    
%   (mut=0.8*pi*sqrt(ra.*rb).*((1.-.5*f).*kk-ek)./k)    		
    dmdk = 0.8*pi*sqrt(ra.*rb).*((ek-(1-.5*f).*kk)./f  ...
        				+((1-.5*f).*dkkdk-k.*kk-dekdk)./k);   
        				 
    dmdra = dmdk .* dkdra + 0.8*pi*0.5*sqrt(rb./ra).*((1.-.5*f).*kk-ek)./k;
    dmdza = dmdk .* dkdza;
    dmdrb = dmdk .* dkdrb + 0.8*pi*0.5*sqrt(ra./rb).*((1.-.5*f).*kk-ek)./k;
    dmdzb = dmdk .* dkdzb;

ind = find(ra==rb & za==zb);
if ~isempty(ind)
   fprintf('ind = %d\n',ind(1));
   dmdra(ind)= Jacsind(ra,delr,delz);
   dmdza(ind) = 0;
   dmdrb(ind) = dmdra(ind);
   dmdzb(ind) = 0;
end

J = [dmdra dmdza dmdrb dmdzb];

if (0)  		% testing only
mut=0.8*pi*sqrt(ra.*rb).*((1.-.5*f).*kk-ek)./k;

dmdra
calc = slope(ra,mut)

dmdza
calc = slope(za,mut)

dmdrb
calc = slope(rb,mut)

dmdzb
calc = slope(zb,mut)

end



end

              
