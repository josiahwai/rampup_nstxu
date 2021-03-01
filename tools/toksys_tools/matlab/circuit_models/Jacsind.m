function dmdr = Jacsind(r,delr,delz)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  dmdr = Jacsind(r,delr,delz)
%
%  PURPOSE:  Calculate Jacobian dmdr for source and msmnt point in same 
%	location.  Assumes rectangular coil.
%
%  INPUT:
%	r    = major radial location (meters)
%	delr = width of rectangular coil
%	delz = height of rectangular coil
%
%  OUTPUT:
%       J = Jacobian (units = Tesla/Amp) matrix, kth row of which holds
%               [dm/d(ra) dm/d(za) dm/d(rb) dm/d(zb)]
%       for the kth element in vectors r.  (uH/m)
%
%  RESTRICTIONS: Not accurate for delz/r > 2.  (Huh?)
%
%  METHOD:  Self inductance does not depend on z position, so derivatives
%	with respect to z are simply set to 0.  does dr depend on whether
%	source or measurement?
%
%  WRITTEN BY:  Mike Walker 	ON 	11/11/96
%	(derived from Jim Leuer's RECTL function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  a=r; b=delz; c=delr;
  amu = 1.2566371;        %mu0=0.4pi to give uH output

% Use fourth order expansion for rectangular section
   d= sqrt(b.*b+c.*c);
   u= b.*b./(c.*c).*log(d.*d./(b.*b));
   v= c.*c./(b.*b).*log(d.*d./(c.*c));
   w= b./c.*atan2(c,b);
   wp= c./b.*atan2(b,c);

   ba2= b.*b./(a.*a);
   ca2= c.*c./(a.*a);
   ba4= ba2.*ba2;
   ca4= ca2.*ca2;
   baca2= ba2.*ca2;
   al= log(8*a./d);

   a3 = a.^3;
   dt1dr = 1./a;		% (= dal/da)
   dt5adr = (3*b^2+c^2)./a3 - ((6*b^2+2*c^2)./a3).*al;
   dt5bdr = -u*b^2./a3;
   dt5cdr = v/5*c^2./a3;
   dt5ddr = 32.*w/5*b^2./a3;
   dt5edr = -69/10*b^2./a3;
   dt5fdr = -442/60*c^2./a3;
   dt5dr = (dt5adr + dt5bdr + dt5cdr + dt5ddr + dt5edr + dt5fdr)/96.0;

   a5 = a.^5;
   dba4dr = -4*b^4./a5;
   dca4dr = -4*c^4./a5;
   dbaca2dr = -4*b^2*c^2./a5;

   dt6adr = (-30.0*dba4dr+35.0*dbaca2dr+22.0*dca4dr/3.0).*al + ...
   		    	(-30.0*ba4+35.0*baca2+22.0*ca4/3.0)./a;
   dt6bdr = -(115.0*dba4dr-480.0*dbaca2dr).*u/12.0;
   dt6cdr = -23.0*dca4dr.*v/28.0;
   dt6ddr = 256.0*(6.0*dba4dr-7.0*dbaca2dr).*w/21.0;
   dt6edr = -(36590.0*dba4dr-2035.0*dbaca2dr-11442.0*dca4dr)/840.0;
   dt6dr  =  (dt6adr + dt6bdr + dt6cdr + dt6ddr + dt6edr )/30720.0;

       t1= al + (1.0+u+v)/12.0 - 2.0*(w+wp)/3.0;

       t5a= (3*ba2+ca2).*al;
       t5b= ba2.*u/2.0;
       t5c= -ca2.*v/10.0;
       t5d= -16.0*ba2.*w/5.0;
       t5e= 69.0*ba2/20.0;
       t5f= 221.0*ca2/60.0;
       t5=  (t5a+t5b+t5c+t5d+t5e+t5f)/96.0;

       t6a= (-30.0*ba4+35.0*baca2+22.0*ca4/3.0).*al;
       t6b= -(115.0*ba4-480.0*baca2).*u/12.0;
       t6c= -23.0*ca4.*v/28.0;
       t6d= 256.0*(6.0*ba4-7.0*baca2).*w/21.0;
       t6e= -(36590.0*ba4-2035.0*baca2-11442.0*ca4)/840.0;
       t6=  (t6a+t6b+t6c+t6d+t6e)/30720.0;

       selfind = amu*a.*(t1+t5+t6);

   dmdr = amu*(t1+t5+t6) + amu*a.*(dt1dr + dt5dr + dt6dr);
