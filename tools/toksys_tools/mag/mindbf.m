function [mut,br,bz]=mindbf(ra,za,rb,zb)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  [mut,br,bz]=mindbf(ra,za,rb,zb) 
%
%  PURPOSE:  Calculate mutual inductances and B-field Green functions
%   in Matlab. Fields are calculated at point b due to current at point a.
%
%  INPUT:
%    ra,za = 	(meters)
%    rb,zb = 	(meters)
%  Allows input of vectors of ra,za,rb,zb.
%
%  OUTPUT:
%    mut = mutual inductance between 2 loops at locations (ra,za), (rb,zb)
%		(Henries)
%    br  = multiplier of current at pt(s) a to produce radial field at pt(s) b 
%		(Tesla/Amp)
%    bz  = multiplier of current at pt(s) a to produce vertical field at pt(s) b 
%		(Tesla/Amp)
%
%  RESTRICTIONS:
%  If vector input, output corresponds to pairings of elements in sets a and b.
%  For example, mut is calculated between (ra(k),za(k)) and (rb(k),zb(k)), but
%  not between (ra(k),za(k)) and (rb(j),zb(j)) for j different from k.

%  METHOD:  
%
%
%  WRITTEN BY:  Dave Humphreys 	ON 	2/12/95
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% jal20061031: mut & br blows up for rb=0 answer should be br=0 below is clug
%    id= find(rb<eps);
%    rb(id)= eps*ones(size(id));

    twopir = 2*pi*rb;
    m=((rb+ra).^2 - (rb-ra).^2)./((rb+ra).^2 + (za-zb).^2);
    k=sqrt(m);
    
    [kk,ek]=ellipke(m);
    mut=8*pi*sqrt(ra.*rb)*1.0e-7.*((1.-.5*m).*kk-ek)./k;

%B-field calculations:
 %Use elliptic integrals for all points:
      dkdk=ek./(k.*(1.-m))-kk./k;
      dedk=(ek-kk)./k;
      dmdk=8*pi*sqrt(ra.*rb)*1.0e-7.*((ek-(1.-.5*m).*kk)./m  ...
          +((1.-.5.*m).*dkdk-k.*kk-dedk)./k);
      dkdz=(za-zb).*sqrt(((rb+ra).^2-(rb-ra).^2)./((rb+ra).^2+(za-zb).^2).^3);
      br = -dmdk.*dkdz./twopir;
      d = (rb+ra).^2 + (za-zb).^2;
      n = (rb+ra).^2 - (rb-ra).^2;
      dddr = 2*(rb + ra);
      dndr = 4*ra;
      dkdr = (.5./k).*(d.*dndr - n.*dddr)./d.^2;
      dmdr = 8*pi*0.5*sqrt(ra./rb)*1.0e-7.*((1.-.5*m).*kk-ek)./k;
      bz = (dmdk.*dkdr + dmdr)./twopir;

     return
