function inductance=selfind(rc,radc,kappa)
 %
%  SYNTAX:  inductance=selfind(rc,radc,kappa)
%
%  PURPOSE:  Self-inductance expression for an axisymmetric coil of major radius RC
%   and minor radius RADC (in uH).
%
%  INPUT:
%	rc   = major radius [m]
%	radc = minor radius [m]
%	kappa = elongation (optional, default=1)
%
%  OUTPUT:
%	inductance = self inductance of coil
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  UNKNOWN  ON 	??
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)selfind.m	1.1 01/24/13

if nargin<3, kappa=1; end
rad = radc*sqrt(kappa);

    inductance=0.4*pi*rc.*((1.+rad.^2./(8*rc.^2))  ...
          .*log(8*rc./rad)+(rad.^2)./(24*rc.^2) - 1.75);

