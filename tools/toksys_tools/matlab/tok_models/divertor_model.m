   function [nec,Wth,Qtard,nec2,nrc2,nrd2,output_data] = ...
                  divertor_model(Gamec,Gamrc,Gamrd,nec1,nrc1,nrd1,dt);
 %
%  SYNTAX:  [nec,Wth,Qtard,nec2,nrc2,nrd2,output_data]=
%		  divertor_model(Gamec,Gamrc,Gamrd,nec1,nrc1,nrd1,dt);
%
%  PURPOSE:  Model evolution of divertor (and core) plasma states in response
%		to changing gas flow rates. Outputs quantities useful to 
%		measure for active feedback control of divertor (and core)
%		properties. Evolves single time step on call
%
%  INPUTS:
%	Gamec = core fueling gas flow rate (hydrogenic, so same as electron)
%			[atoms/m^3/sec = electrons/m^3/sec]
%	Gamrc = core radiating impurity (RI) flow rate [atoms/m^3/sec]
%	Gamrd = divertor radiating impurity flow rate [atoms/m^3/sec]
%	nec1 = initial core electron density [electrons/m^3]
%	nrc1 = initial core RI density [ions/m^3] 
%	nrd1 = initial divertor RI density [ions/m^3] 
%	dt = time step [sec]
%
%  OUTPUTS:    
%	nec = present core electron density output (same as nec2) [electrons/m^3]
%	Wth = core thermal stored energy output [J]
%	Qtard = heat flux to divertor target output [W/s]
%	nec2 = core electron density after time step dt [electrons/m^3] 
%	nrc2 = core RI density after time step dt [ions/m^3] 
%	nrd2 = divertor RI density after time step dt [ions/m^3] 
%	output_data = structure with useful data to provide to analyze output:
%		Pradc, Pradd
%
%  RESTRICTIONS:
%
%  METHOD:  
%   State equations are the three continuity equations for nec, nrc, nrd.
%   Output equations describe Wth Qtard as function of input heating power Pheat and
%	radiated core power Pradc. nec is the third output variable. 
%   Pradc represented by Bremsstrahlung-like expression including both core
%	electron-electron and electron-impurity collisions. Pradd (radiation
%	from divertor) represented by empirical function to capture observed
%	behavior of ionization/density front up and along divertor legs, with
%	corresponding maximum in radiation as function of upstream density (nec).
%	Scale factor of 150 makes function produce ~ 1 MW total rad pwr
%	for 2 keV, 2e19 plasma w/ no impurity... 
%   Pradd coeff Crd is chosen so 2e19 core density makes ~1.5 MW rad in div
%   State equations evolved by implicit scheme (exponential solution over time step)
%   

%  WRITTEN BY:  Dave Humphreys 	ON	1/16/13
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters:
  tauec = 1.00;    %core e conf time [s]
  taurc = 2.00;    %core RI conf time [s]
  taurd = 1.00;   %div RI conf time [s]
  taurcd = 2.00;   %characteristic diff. time from div to core [s]
  Zr =  18;    %RI effective charge state
  etarcd = 0.1;   %RI div to core transport efficiency
  Tec = 2000;    %core electron temp [eV]
  Pheatc = 5.0e6;    %core heating power [W]
  fdiv = 0.75;     %frac of SOL power reaching divertor
  nec00 = 2.0e19;    %core e density for peak div radiation
  signec =  1.0e19;   %characteristic width of Pradd function [particles/m^3]
  Crd = 0.5e6/2e18^2;     %coeff of Pradd function [W-m^6]
  tauE = 0.100;      %core energy conf time [s]
  Vcore = 20;	    %core plasma volume [m^3]
  scalerad = 4;    %scale factor for Pradc from simple Bremsstrahlung
  
% Derived values:


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolve state equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Core electron density:
  RHS = Gamec + Zr*Gamrc + etarcd*Zr*nrd1/taurcd;
  nec2 = nec1 + (RHS - nec1)*(1 - exp(-dt/tauec));
  nec = nec2;   %set output value to latest state value
  
% Core RI density:
  RHS = Gamrc + etarcd*nrd1/taurcd;
  nrc2 = nrc1 + (RHS - nrc1)*(1 - exp(-dt/taurc));
  nrc = nrc2;
  
% Divertor RI density:
  RHS = Gamrd - etarcd*nrd1/taurcd;
  nrd2 = nrd1 + (RHS - nrd1)*(1 - exp(-dt/taurd));
  nrd = nrd2;
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Total radiated power from core:
  Pradc = scalerad*1.69e-38*nec*sqrt(Tec)*(nec + Zr^2*nrc)*Vcore;   %[W]
  
%Power radiated from divertor:
  Pradd = Crd*nrd^2*exp( -(nec-nec00)^2/(2*signec^2) );           %[W]
 
%Core stored thermal energy:
  Wth = tauE*(Pheatc - Pradc);
  
%Heat flux to divertor target:
  Qtard = fdiv*(Pheatc - Pradc) - Pradd;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct output_data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  output_data = struct( 'Pradc', Pradc, ...
  			'Pradd', Pradd    ); 
  




