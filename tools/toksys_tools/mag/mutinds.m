  function [M,Mr,Mz,Mrr,Mrz,Mzz,Mrrr,Mrrz,Mrzz,Mzzz] = mutinds(ra, za, rb, zb, D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  M = mutinds(ra,za,rb,zb)
%           [M,Mr,Mz] = mutinds(ra,za,rb,zb)
%           [M,Mr,Mz,Mrr,Mrz,Mzz] = mutinds(ra,za,rb,zb)
%           [M,Mr,Mz,Mrr,Mrz,Mzz,Mrrr,Mrrz,Mrzz,Mzzz] = mutinds(ra,za,rb,zb)
%           To obtain just a derivative, specify one of 'r','z', etc as fifth 
%             argument, example: Mzz = mutinds(ra,za,rb,zb,'zz')
%
%  PURPOSE:  Calculate mutual inductance between loops centered at r=0
%
%  INPUTS: ra,za = coordinates of first current loop [meters]
%	   rb,zb = coordinates of second current loop [meters]
%
%  OUTPUTS: M = mutual inductance between loops a and b [H]
%           Mr, Mz, etc = partial derivatives w.r.t. rb and zb

%  RESTRICTIONS: 

%  METHOD:  
%
%  VERSION @(#)mutinds.m	1.1 02/22/15
%
%  WRITTEN BY: Anders Welander ON 2015-02-10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = 8e-7*pi*sqrt(ra.*rb);
a = (rb+ra).^2 + (zb-za).^2;
b = (rb+ra).^2 - (rb-ra).^2;
m = b./a;
k = sqrt(m);

[K,E] = ellipke(m);
% dKdk = E./(k.*(1-m))-K./k;
% dEdk = (E-K)./k;

M = f.*((K-E)./k-k.*K/2);

% Calculating first derivatives
if nargout > 1 | exist('D','var') & length(D) > 0
  g = rb.*ra./a;
  h = sqrt(g);
  dMdk_fE = 1/2./(1-m) + 1./m;
  dMdk_fK = -1./m;
  dMdk = f.*(dMdk_fE.*E+dMdk_fK.*K);
  dkdr = (ra - 2*(rb+ra).*g)./h./a;
  dkdz = 2*h.*(za-zb)./a;
  Mz = dMdk.*dkdz;      
  Mr = M/2./rb + dMdk.*dkdr;
end

% Calculating second derivatives
if nargout > 3 | exist('D','var') & length(D) > 1
  d2Mdk2_fEk = 1./(1-m).^2 - 1/2./(1-m)./m - 2./m.^2;
  d2Mdk2_fKk = -1/2./(1-m) + 2./m.^2 - 1/2./m;
  d2Mdk2 = f.*k.*(d2Mdk2_fEk.*E+d2Mdk2_fKk.*K);
  d2kdrdz = (zb-za)./a.*(4*(rb+ra).*g./a./h - dkdr);
  d2kdz2 = (3*(za-zb).*dkdz - 2*h)./a;
  Mzz = d2Mdk2.*dkdz.^2 + dMdk.*d2kdz2;
  Mrz = Mz/2./rb + d2Mdk2.*dkdr.*dkdz + dMdk.*d2kdrdz;      
  Mrr = Mr./rb-Mzz;
end

% Calculating third derivatives
if nargout > 6 | exist('D','var') & length(D) > 2
  d3Mdk3_fE = 1/2./(1-m).^2 + 4*m./(1-m).^3 + 3/2./(1-m)./m + 6./m.^2;
  d3Mdk3_fK = -(1/2+3/2*m)./(1-m).^2 - (6-3/2*m)./m.^2;
  d3Mdk3 = f.*(d3Mdk3_fE.*E+d3Mdk3_fK.*K);
  dhdr = (1/2*ra./a - g./a.*(rb+ra))./h;
  dhdz = -g./a.*(zb-za)./h;
  d3kdrdz2 = (3*(za-zb).*d2kdrdz - 2*dhdr - 2*(rb+ra).*d2kdz2)./a;
  d3kdz3 = (5*(za-zb).*d2kdz2-3*dkdz - 2*dhdz)./a;
  Mzzz = d3Mdk3.*dkdz.^3 + 3*d2Mdk2.*dkdz.*d2kdz2 + dMdk.*d3kdz3;
  Mrzz = Mzz/2./rb + dMdk.*d3kdrdz2 + d3Mdk3.*dkdr.*dkdz.^2 + ...
    d2Mdk2.*(dkdr.*d2kdz2+ 2*d2kdrdz.*dkdz);
  Mrrz = Mrz./rb-Mzzz;
  Mrrr = (Mrr-Mr./rb)./rb-Mrzz;
end

if exist('D','var') & ~isempty(D)
  if strcmp(lower(D),'r')
    M = Mr;
  elseif strcmp(lower(D),'z')
    M = Mz;
  elseif strcmp(lower(D),'rr')
    M = Mrr;
  elseif strcmp(lower(D),'rz')
    M = Mrz;
  elseif strcmp(lower(D),'zz')
    M = Mzz;
  elseif strcmp(lower(D),'rrr')
    M = Mrrr;
  elseif strcmp(lower(D),'rrz')
    M = Mrrz;
  elseif strcmp(lower(D),'rzz')
    M = Mrzz;
  elseif strcmp(lower(D),'zzz')
    M = Mzzz;
  else
    disp('Warning mutinds: Fifth input not recognized. Allowed values are:')
    disp('''r'',''z'',''rr'',''rz'',''zz'',''rrr'',''rrz'',''rzz'',''zzz''')
  end
end
