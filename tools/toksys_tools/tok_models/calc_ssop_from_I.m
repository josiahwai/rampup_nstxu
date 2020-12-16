function ssop = calc_ssop_from_I(tok_system,ssop_inputs)
 %
%  SYNTAX:  ssop = calc_ssop_from_I(tok_system,ssop_inputs)
%
%  PURPOSE: Compute steady-state operating point (ynom,Inom,Vnom), where Inom contains
%	closest approximation to Icnom consistent with being a valid steady- state operating 
%	point.  Values of ynom are always nonzero when default values of zero nominal 
%	currents are used.  
%
%  INPUT:
%   tok_system = tokamak system model, as constructed by build_tokamak_system
%   ssop_inputs = structure containing all optional inputs:
%     Icnom = vector of PF coil currents defining center of coordinates for model.  To use
%	plasma equilibrium currents, let Icnom = Pcc*Ieq(1:ncx)   (optional, default = all 0)
%     Ipnom = value of plasma current defining center of coordinates for model (optional, 
%	default = 0). Nonzero value should be used only if non-ohmic drive is available to 
%	hold plasma current = Ipnom in steady-state.
%
%  OUTPUT:
%   ssop = steady state operating point, containing:
%      ynom = nominal values for outputs y
%      Inom = vector of nominal coil, vessel, and plasma currents (coil currents in 
%		Inom are best match to Icnom that is consistent with steady-state gain)
%      Vnom = vector of voltages that produce coil currents in Inom vector
%      Inomstates = vector of nominal coil, vessel, and plasma current states
% 
%  RESTRICTIONS:  Only works right now for outputs that are modeled as linear relative to a
%	plasma equilibrium.
%
%  METHOD:  See Prop. 3 of Walker & Humphreys, Valid Coordinate Systems for Linearized Plasma
%	Shape Response Models in Tokamaks (where it is shown that the only valid choice for 
%	nominal vessel currents is 0).  The Icnom vector is projected onto the range of the 
%	steady-state gain to compute Inom; the pseudo-inverse of steady-state gain is applied 
%	to compute Vnom (this is true inverse for vectors in range of steady-state gain).
 
%  WRITTEN BY:  Mike Walker 	ON 	4/10/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @(#)calc_ssop_from_I.m	1.2 04/29/11

% TO DO:
%   Add constrained calculation (to stay inside coil limits)

if(nargin>1)
   struct_to_ws(ssop_inputs);
end

Pcc = tok_system.Pcc;
if(~exist('Ipnom'))
   Ipnom = 0;
end
if(~exist('Icnom'))
   Icnom = zeros(size(Pcc,1),1);
end

cmat = tok_system.cmat;
ncx = tok_system.ncx;
yeq = tok_system.yeq;
Ieq = tok_system.Ieq;
iplcirc = tok_system.flags.iplcirc;

% Handle models with or without plasma circuit (including vacuum models):
if(iplcirc==0)
   if(Ipnom ~= 0)
      wait('ERROR calc_ssop_from_I: Cannot specify nonzero plasma current with iplcirc=0')
      return;
   end
end

% Project Icnom onto range of steady-state gain matrix.
ic1 = tok_system.output_map.cc.start;
ic2 = tok_system.output_map.cc.end;
ssgain = dcgain(tok_system.amat,tok_system.bmat,tok_system.cmat(ic1:ic2,:),tok_system.dmat(ic1:ic2,:));
range_ssgain = orth(ssgain);
Inom = range_ssgain*range_ssgain'*Icnom;	% steady-state currents
pinv_ssgain = pinv(ssgain);
Vnom = pinv_ssgain*Inom;			% steady-state voltages
Inomstates = pinv(Pcc)*Inom;

if(iplcirc) 	% Append vessel and plasma nominal currents
   nv = size(tok_system.Pxx,1) - length(Inom) - 1;
   Inom = [Inom; zeros(nv,1); Ipnom];
   Inomstates = [Inomstates; zeros(tok_system.nvx,1); Ipnom];
else
   nv = size(tok_system.Pxx,1) - length(Inom);
   Inom = [Inom; zeros(nv,1)];
   Inomstates = [Inomstates; zeros(tok_system.nvx,1)];
end

% This is calculation for outputs that are defined as relative to plasma equilibrium:
ynom = yeq + cmat*(Inomstates-Ieq);

% Calculations for outputs that are originally defined in absolute coordinates.

% All currents constant => all loop voltages are zero:
ynom(tok_system.output_map.lv.start:tok_system.output_map.lv.end) = 0;

% Netlist models only the dynamic equation => all outputs are absolute:
ynom(tok_system.output_map.nl.start:tok_system.output_map.nl.end) =  ...
	tok_system.netlist_model.Chat * Inomstates + tok_system.netlist_model.Dhat * Vnom;

ssop = struct( ...
'ynom',ynom, ...
'Inom',Inom, ...
'Vnom',Vnom, ...
'Inomstates',Inomstates, ...
'ssgain',ssgain, ...
'pinv_ssgain',pinv_ssgain);
