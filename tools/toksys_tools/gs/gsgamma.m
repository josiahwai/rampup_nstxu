function gamma = gsgamma(init,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gamma = gsgamma(eq,config)
%
%  PURPOSE: Calculate vertical growth rate for an equilibrium
%
%  INPUTS:     eq, equilibrium,
%          config, tokamak description, most relevant fields are:
%                  mpc, mpv, mcc, mcv, mvv = mutual inductances
%                  resc, resv = resistances in coils and vessel elements
%    Coils can be connected in circuits by one of 3 methods:
%    config.icci, config.cccirc, config.buscode
%    These are explained below where ic is coil and ci circuit currents
%      ic = config.icci*ci,
%      ic = sign(config.cccirc).*ci(abs(cccirc))
%      i = find(config.buscode), icci = eye(nc), icci(i(1),i(2:end)) = -1
%          icci = icci(:,[1:i(1)-1 i(1)+1:nc]), ic = icci*ci
%
%  OUTPUTS: gamma, vertical growth rate [rad/s]
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  WRITTEN BY:  Anders Welander ON 2016-04-28
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

config.outputs = 'gamma';
config.pres0 = init.pres;
config.fpol0 = init.fpol;
config.no_edge_gradient = 0;
config.no_edge_current = 0;
config.evolve_option = 1;
gs_configure
Rpla_previous = plares.values(1);
gs_initialize
gs_eq_analysis
gs_response
xs = [];
gs_dynamics
plot(pprime)
