function [J]=BJacobians(ra,za,rb,zb)
% NOT TESTED
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  J=BJacobians_source(ra,za,rb,zb)
%
%  PURPOSE:  Compute Jacobian d(Br,Bz)/d(ra,za) , where Br and Bz are the
% 	radial field components evaluated at (rb,zb) due to current source
%	at (ra,za)%.
% 	!!! NOTE: THIS DOES NOT COMPUTE d(Br,Bz)/d(rb,zb) !!!
%
%  INPUT:
%	ra,za = (r,z) coordinates of current source (vectors)
%	rb,zb = (r,z) coordinates of point where (Br,Bz) are measured (scalar)
%   (units = meters)
%
%  OUTPUT:
%	J = Jacobian = [dbrdr dbrdz dbzdr dbzdz], (units = Tesla/Amp/m)
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	10/25/96
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debug=1;

% from Dave's code (modified)
    m=((rb+ra).^2 - (rb-ra).^2)./((rb+ra).^2 + (za-zb).^2);
    k=sqrt(m);
    
    [kk,ek]=ellipke(m);

    dkdk=ek./(k.*(1.-m))-kk./k;
    dedk=(ek-kk)./k;
    dkdz=(zb-za).*sqrt(((rb+ra).^2-(rb-ra).^2)./((rb+ra).^2+(za-zb).^2).^3);
    d = (rb+ra).^2 + (za-zb).^2;
    n = (rb+ra).^2 - (rb-ra).^2;
    dddr = 2*(rb + ra);
    dndr = 4*rb;
    dkdr = (.5./k).*(d.*dndr - n.*dddr)./d.^2;

dmdk=8*pi*sqrt(ra.*rb)*1.0e-7.*((ek-(1.-.5*m).*kk)./m  ...
		+((1.-.5.*m).*dkdk-k.*kk-dedk)./k);
twopir = 2*pi*rb;
dmdr = 8*pi*0.5*sqrt(ra./rb)*1.0e-7.*((1.-.5*m).*kk-ek)./k;

br0 = -dmdk.*dkdz./twopir;
bz0 = (dmdk.*dkdr + dmdr)./twopir;

% from psical.f
z = zb-za;
d = (rb+ra).^2 + (za-zb).^2;

kb = (ra-rb).^2 + z.^2;
k1 = (ra.^2+rb.^2+z.^2) ./ kb;
k2 = (ra.^2-rb.^2-z.^2) ./ kb;

br = z./(rb.*sqrt(d)) .* (-kk + k1.*ek)*2e-7;
bz = (kk + k2.*ek)./sqrt(d)*2e-7;

%if debug		% test that derivatives above all OK
%   fprintf('norm(br-br0) = %f\n', norm(br-br0)/norm(br));
%   fprintf('norm(bz-bz0) = %f\n', norm(bz-bz0)/norm(bz));
%end

% Extend to Jacobian calculations

dkbdr = 2*(ra-rb);
dkbdz = 2*(za-zb);
dddz = dkbdz;

dk1dr = (2*kb .* ra  - (ra.^2+rb.^2+z.^2).*dkbdr) ./ kb.^2;
dk2dr = (2*kb .* ra  - (ra.^2-rb.^2-z.^2).*dkbdr) ./ kb.^2;

dk1dz = (2*kb .* (za-zb) - (ra.^2+rb.^2+z.^2).*dkbdz) ./ kb.^2;
dk2dz = (2*kb .* (zb-za) - (ra.^2-rb.^2-z.^2).*dkbdz) ./ kb.^2;

dbrdd = -0.5 * z./(rb.* (sqrt(d)).^3) .* (-kk + k1.*ek) * 2e-7;
dbrdkk = -z./(rb.*sqrt(d)) * 2e-7;
dbrdk1 = z./(rb.*sqrt(d)) .* ek * 2e-7;
dbrdek = z./(rb.*sqrt(d)) .* k1 * 2e-7;

dbzdd = -0.5./(sqrt(d).^3) .* (kk + k2.*ek) * 2e-7;
dbzdkk = 1./sqrt(d) * 2e-7; 
dbzdk2 = ek./sqrt(d) * 2e-7;
dbzdek = k2./sqrt(d) * 2e-7;


dbrdr = dbrdd.*dddr + dbrdkk.*dkdk.*dkdr + dbrdek.*dedk.*dkdr + dbrdk1.*dk1dr;
dbrdz = dbrdd.*dddz + dbrdkk.*dkdk.*dkdz + dbrdek.*dedk.*dkdz +dbrdk1.*dk1dz ...
	 			- 1./(rb.*sqrt(d)) .* (-kk + k1.*ek)*2e-7; 
dbzdr = dbzdd.*dddr + dbzdkk.*dkdk.*dkdr + dbzdek.*dedk.*dkdr + dbzdk2.*dk2dr;
dbzdz = dbzdd.*dddz + dbzdkk.*dkdk.*dkdz + dbzdek.*dedk.*dkdz + dbzdk2.*dk2dz;

J = [dbrdr dbrdz  dbzdr dbzdz];

if (0) then 		% testing only
dbrdr
calc = slope(ra,br)

dbrdz
calc = slope(za,br)

dbzdr
calc = slope(ra,bz)

dbzdz
calc = slope(za,bz)

end



end

              
