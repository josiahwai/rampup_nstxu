%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gsevolve_init
%
%  PURPOSE: Template for initialization of gsevolve
%           All initialization options are included
%
%  INPUTS:  none
%
%  OUTPUTS: init, the variable used to initialize gsevolve
	
%  VERSION @(#)gsevolve_init.m	1.1 02/23/15
%
%  WRITTEN BY:  Anders Welander  ON	1/27/15
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PLEASE RENAME AND SAVE A LOCAL COPY OF THIS TEMPLATE

% The initial equilibrium can be on a different grid than config

% INITIALIZE FROM EXISTING EQUILIBRIUM
% Initialization can be very easy, simply load an equilibrium:
% init = read_eq(161086,4,'EFIT02','D3D');
% That's it! No need to read more unless you are interested.

% The rest of the template describes other ways to initialize.
% This includes modifying the initial state or starting from no plasma
 
% The state of the equilibrium is given by the vector 
%   xs = [circuit currents; parameters for internal profiles]
% Exact format of circuits and profile parameters depends on the configuration

% Initialization of gsevolve is about initializing the vector xs.
% However, if only the vector xs would be given on a first call then gsevolve  
% would need to find the corresponding equilibrium from scratch.
% For some xs vectors there are more than one equilibrium solution and
% for some there is none. Therefore an at least approximate solution
% psizr [Wb] should be given if the initial state has a plasma.
% The initial state can be created with the code gsdesign.m

% Hence init should include:
%  1. Initial values of parameters that are related to xs
%  2. An equilibrium (psizr) that at least roughly 
%     corresponds to initial parameters
%     init.rmaxis, init.zmaxis are also used if available

% Initialize all of xs with init.xs (highest precedence)
% This requires knowledge of what the states are
% If init.xs exists then all parameters in the following are ignored

% Initialize coil currents in order of precedence:
% init.ci, circuit currents
% init.ic, coil currents that make flux = config.mpc*init.ic [Wb]
% init.cc, EFIT coil currents (are converted to ic by cc_efit_to_tok.m)

% Initialize vessel currents in order of precedence:
% init.iv, vessel currents that make flux = config.mpv*init.iv [Wb]
% init.cc, EFIT vessel currents (converted to iv by cc_efit_to_tok.m)

% By default all profile parameters are given at values of normalized
% poloidal flux values = linspace(0,1,length(profile_parameter))

% Initialize pressure profile in order of precedence:
% init.betap (if this is a state)
% init.Wth (if this is a state) (total thermal energy)
% init.W (if these are states) (thermal energy within flux surfaces)
% init.pres, pressure values [Pa]
% init.pprime [2*pi*Pa/Wb], with init.psimag [Wb], init.psibry [Wb]

% Initialize current profile in order of precedence:
% init.cpasma, init.li (if these are states)
% init.I (if these are states) (toroidal current within flux surfaces)
% init.fpol, [Tm] [Tesla*meters]
% init.ffprim [2*pi*Tm/Wb], with init.psimag [Wb], init.psibry [Wb]
% NOTE! Unless init.fpol is given, init.rzero, init.bzero are required
% where init.bzero is applied toroidal field measured at init.rzero


% VACUUM INITIALIZATION
% For vacuum initialization pressure and current profiles are set to zero
% init.pres = zeros(65,1);
% init.fpol = zeros(65,1)+init.rzero*init.bzero;
% The applied toroidal field init.bzero should still be non-zero
% init.rmaxis, init.zmaxis, if supplied dictates point of plasma initiation
% If these are missing or nan the initiation point is found by gs_breakdown

% REINITIALIZATION
% The state can at any time be changed by calling gsevolve with new init
