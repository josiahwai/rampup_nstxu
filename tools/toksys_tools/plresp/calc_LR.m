   function [tau,Lp,Res,eta,Lp0] = calc_LR(R0,a,kappa,li,Te,Zeff,lnlam,ne)
 %
%  SYNTAX:  [tau,Lp,Res,eta,Lp0] = calc_LR(R0,a,kappa,li,Te,Zeff,lnlam,ne)
%
%  PURPOSE:  Calculate decay time for plasma
%
%  INPUTS:
%	R0 = major radius [m]
%	a = minor radius [m]
%	kappa = elongation
%	li = (optional) internal inductance (default = 0.85)
%	Te = electron temperature [eV]
%	Zeff = (optional) effective charge state (default=1)
%	lnlam = (optional) Coulomb logarithm ln(Lambda) (default=15)
%	ne = (optional) electron density [m^-3] (if present, calculates lnlam
%		from lnlam = 31 - ln(sqrt(ne)/Te) and ignores input lnlam)
%
%  OUTPUTS:
%	tau = resistive decay L/R time constant [s]
%	Lp = plasma self inductance (calculated with sind.m) [H]
%	Res = plasma resistance [Ohm]
%	eta = resistivity [Ohm-m]
%	Lp0 = Shafranov expression with kappa correction
%
%  RESTRICTIONS:
%	All quantities scalars now.
%
%  METHOD:  
%	From Chen, p. 183:
%    	    eta = 5.2e-5*Zeff*lnlam/(Te^1.5); %Pl.restvty ohm-m Spitzer
%	From Hutchinson, p. 19:
%	    lnlam = 31 - log(sqrt(ne)/Te);  %valid for Te>10 eV

%  WRITTEN BY:  Dave Humphreys 	ON	7/9/98
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)calc_LR.m	1.4 01/24/13

% Constants:
  mu0 = 0.4*pi;   %

% Defaults:
  if ~exist('li')
     li = 0.85;
  end
  if ~exist('Zeff')
    Zeff=1;
  end
  if ~exist('lnlam')
    lnlam = 15;
  end
  if nargin>7
    lnlam = 31 - log(sqrt(ne)/Te);
  end

% Resistivity, Resistance:
    eta = 5.2e-5*Zeff*lnlam/(Te^1.5);      %Pl.restvty ohm-m Spitzer
    Res = 2*pi*R0*eta/(pi*kappa*a^2);       %Pl.res Ohms(Spitzer)

% Inductance (using sind.m):
    Lp = selfind(R0,sqrt(kappa*a^2))*1e-6;     %Pl Inductance, H
    Lp0 = 1e-6*mu0*R0*(log(8*R0/(a*sqrt(kappa))) + 0.5*li -2);  %H

% Decay time:
    tau = Lp/Res;

             
