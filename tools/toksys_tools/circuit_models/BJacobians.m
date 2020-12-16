function [J]=BJacobians(ra,za,rb,zb)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  J=BJacobians(ra,za,rb,zb)
%
%  PURPOSE:  Compute Jacobian d(Br,Bz)/d(rb,zb) , where Br and Bz are the
% 	radial field components evaluated at (rb,zb) due to current source
%	of 1 Amp at (ra,za).
%
% 	!!! NOTE: THIS DOES NOT COMPUTE d(Br,Bz)/d(ra,za) !!!
%
%  INPUT:
%	ra,za = (r,z) coordinates of current source (vectors)
%	rb,zb = (r,z) coordinates of point where (Br,Bz) are measured (scalar)
%   (units = meters)
%
%  OUTPUT:
%	J = Jacobian = [dbr/drb dbr/dzb dbz/drb dbz/dzb], (units = Tesla/Amp/m)
 
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

    kk_k=ek./(k.*(1.-m))-kk./k;
    ek_k=(ek-kk)./k;
    k_zb=(za-zb).*sqrt(((rb+ra).^2-(rb-ra).^2)./((rb+ra).^2+(za-zb).^2).^3);
    d = (rb+ra).^2 + (za-zb).^2;
    n = (rb+ra).^2 - (rb-ra).^2;
    d_rb = 2*(rb + ra);
    n_rb = 4*ra;
    k_rb = (.5./k).*(d.*n_rb - n.*d_rb)./d.^2;

psi_k=8*pi*sqrt(ra.*rb)*1.0e-7.*((ek-(1.-.5*m).*kk)./m  ...
		+((1.-.5.*m).*kk_k-k.*kk-ek_k)./k);
twopir = 2*pi*rb;
psi_rb = 8*pi*0.5*sqrt(ra./rb)*1.0e-7.*((1.-.5*m).*kk-ek)./k;

br0 = -psi_k.*k_zb./twopir;
bz0 = (psi_k.*k_rb + psi_rb)./twopir;

% from psical.f
z = zb-za;
d = (rb+ra).^2 + (za-zb).^2;

kb = (ra-rb).^2 + z.^2;
k1 = (ra.^2+rb.^2+z.^2) ./ kb;
k2 = (ra.^2-rb.^2-z.^2) ./ kb;

br = z./(rb.*sqrt(d)) .* (-kk + k1.*ek)*2e-7;
bz = (kk + k2.*ek)./sqrt(d)*2e-7;

if debug  % test that derivatives above all OK
   fprintf('norm(br-br0) = %f\n', norm(br-br0)/norm(br));
   fprintf('norm(bz-bz0) = %f\n', norm(bz-bz0)/norm(bz));
end

% Extend to Jacobian calculations

kb_rb = 2*(rb-ra);	% close
kb_zb = 2*(zb-za);	% close
d_zb = kb_zb;		% close

k1_rb = ( 2*kb .* rb  - (ra.^2+rb.^2+z.^2).*kb_rb) ./ kb.^2;	% close
k2_rb = (-2*kb .* rb  - (ra.^2-rb.^2-z.^2).*kb_rb) ./ kb.^2;	% close

k1_zb = (2*kb .* (zb-za) - (ra.^2+rb.^2+z.^2).*kb_zb) ./ kb.^2;	% very close
k2_zb = (2*kb .* (za-zb) - (ra.^2-rb.^2-z.^2).*kb_zb) ./ kb.^2;	% very close

br_d = -0.5 * z./(rb.* (sqrt(d)).^3) .* (-kk + k1.*ek) * 2e-7;	% perfect
br_kk = -z./(rb.*sqrt(d)) * 2e-7;
br_k1 = z./(rb.*sqrt(d)) .* ek * 2e-7;
br_ek = z./(rb.*sqrt(d)) .* k1 * 2e-7;

bz_d = -0.5./(sqrt(d).^3) .* (kk + k2.*ek) * 2e-7;		% perfect
bz_kk = 1./sqrt(d) * 2e-7; 
bz_k2 = ek./sqrt(d) * 2e-7;
bz_ek = k2./sqrt(d) * 2e-7;


br_rb = br_d.*d_rb + br_kk.*kk_k.*k_rb + br_ek.*ek_k.*k_rb + br_k1.*k1_rb ...
				- z./(rb.^2.*sqrt(d)).*(-kk + k1.*ek)*2e-7; %OK
br_zb = br_d.*d_zb + br_kk.*kk_k.*k_zb + br_ek.*ek_k.*k_zb + br_k1.*k1_zb ...
	 			+ 1./(rb.*sqrt(d)) .* (-kk + k1.*ek)*2e-7; %OK
bz_rb = bz_d.*d_rb + bz_kk.*kk_k.*k_rb + bz_ek.*ek_k.*k_rb + bz_k2.*k2_rb; %OK
bz_zb = bz_d.*d_zb + bz_kk.*kk_k.*k_zb + bz_ek.*ek_k.*k_zb + bz_k2.*k2_zb; %OK

J = [br_rb br_zb  bz_rb bz_zb];

if (0) then 		% testing only
br_rb
calc = slope(ra,br)

br_zb
calc = slope(za,br)

bz_rb
calc = slope(ra,bz)

bz_zb
calc = slope(za,bz)

end
