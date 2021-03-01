function selfind = rectl(rc,delz,delr)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  selfind = rectl(rc,delz,delr)
%
%  PURPOSE:  Calculate self inductance of rectangular coil. 
%
%  INPUT:
%      rc   =  mean radius of the coil (meters)
%      delz =  axial height of the coil (meters)
%      delr =  radial width of the coil (meters)
%
%  OUTPUT:
%      selfind = self inductance (in uH)
%
%  RESTRICTIONS:  Accuracy of calculated self inductance will be better than
%		  1% for delz/(2*delr) < 1.0.  This is not guaranteed for
%		  delz/(2*delr) > 1.0.

%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	5/9/95
%   (Taken from Leuer's FORTRAN code RECTL.FOR, via Dave Humphrey's rectl.pro.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   a=rc; b=delz; c=delr;
%     amu = 4.e-7*!pi
     amu = 1.2566371E-06;

%  Caution: 1% accuracy for b/(2*a) < 1.0 (L_compar_fld.com)
   if(any(b./(2.*a) > 1.0)) 
       [mm,mi]=max(b./(2*a))
       fprintf('  CAUTION rectl: Worse than 1%% L accuracy, b/2a= %f\n',max(b./(2*a)))
   end;    
% ***************************************************************************
%                   USE GROVERS FORMULA FOR SQUARE SECTION
% ***************************************************************************
%   if(b eq c) then begin
%       CO2A= C/(2.*A)
%       CO2A= CO2A*CO2A
%       selfind = AMU*A*(0.5*(1.+CO2A/6.)*ALOG(8./CO2A)-0.84834+.2041*CO2A)
%       goto,theend
%   endif

% ***************************************************************************
%                   use fourth order expansion for rectangular section
% ***************************************************************************

       d= sqrt(b.*b+c.*c);
       u= b.*b./c./c.*log(d.*d./b./b);
       v= c.*c./b./b.*log(d.*d./c./c);
       w= b./c.*atan2(c,b);
       wp= c./b.*atan2(b,c);

       ba2= b.*b./(a.*a);
       ca2= c.*c./(a.*a);
       ba4= ba2.*ba2;
       ca4= ca2.*ca2;
       baca2= ba2.*ca2;
       al= log(8.0.*a./d);

       t1= al + (1.0+u+v)./12.0 - 2.0.*(w+wp)./3.0;

       t5a= (3.*ba2+ca2).*al;
       t5b= ba2.*u./2.0;
       t5c= -ca2.*v./10.0;
       t5d= -16.0.*ba2.*w./5.0;
       t5e= 69.0.*ba2./20.0;
       t5f= 221.0.*ca2./60.0;
       t5=  (t5a+t5b+t5c+t5d+t5e+t5f)./96.0;

       t6a= (-30.0.*ba4+35.0.*baca2+22.0.*ca4./3.0).*al;
       t6b= -(115.0.*ba4-480.0.*baca2).*u./12.0;
       t6c= -23.0.*ca4.*v./28.0;
       t6d= 256.0.*(6.0.*ba4-7.0.*baca2).*w./21.0;
       t6e= -(36590.0.*ba4-2035.0.*baca2-11442.0.*ca4)./840.0;
       t6=  (t6a+t6b+t6c+t6d+t6e)./30720.0;

       selfind = amu.*a.*(t1+t5+t6);
       selfind = 1e6*selfind;			% convert from H to uH
%----------------------------------------------------
%The End:
%    theend:
       return
