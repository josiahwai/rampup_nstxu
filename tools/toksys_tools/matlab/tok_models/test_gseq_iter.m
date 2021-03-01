%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	
%  USAGE:   test_gseq_iter
%
%  PURPOSE: Examples on how to run gseq
%
%  INPUTS:  None (just run the code)
%
%  OUTPUTS: Demonstration of gseq, including plots of evolving equilibrium
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  VERSION @(#)test_gseq_iter.m	1.2 10/25/13
%
%  WRITTEN BY:  Anders Welander  ON	10/07/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid in simulation need not match grid in original equilibrium
load /m/GAtools/tokamaks/iter/make/2010v3p3/iter_obj_2010v3p3_1733
%load /m/GAtools/tokamaks/iter/make/2010v3p3/iter_obj_2010v3p3_3365
%load /m/GAtools/tokamaks/iter/make/2010v3p3/iter_obj_2010v3p3_65129

% The config structure is essentially parts of tok_data_struct, with some additions
config = tok_data_struct;


% Define output (or diagnostic) matrix, see help gseq
config.Cmat = [[config.mcc config.mcv config.mpc']; [config.mcv' config.mvv config.mpv']];
config.Cmat = [[zeros(config.nc,config.nc+config.nv) config.mpc']; ...
               [zeros(config.nv,config.nc+config.nv) config.mpv']];


% Defaults are constraints = 1 and nkn = 1 (knots for splines)
config.constraints = 1; 
config.nkn = 4;

config.plotit = 1; % Will cause gseq to plot progress by default


% gseq can be configured separately
gseq([], [], config) % y is not returned here since it has no meaning

% Read equilibrium that will be used to initialize the simulation
init = read_corsica_flat_files('/m/GAtools/tokamaks/iter/corsica/v3.3/sob');

% Initialize gseq to equilibrium 'init'
gseq([], init); % Now y is defined and therefore returned (prevented from echo with ';')

x = cc_efit_to_tok(config,init); % creates cc0t, vc0t
ic6_0 = x.cc0t(6); % Remember this original current since we will play with it later


% Define initial targets, we will converge eq0 with these targets as constraints
x.ip = 15e6;
x.li = 1;
x.betap = 0.7;
x.plotit = 1; % This flag has precedence over the default flag, config.plotit

% Initialize gseq with x in the call
gseq(x, init);
% In this case 'init' is only used as default for information that is missing from x


% Alternatively, configure and initialize, all at once
[y, e, eq, eqx] = gseq(x, init, config);
% This time also storing all possible outputs


% Remember this equilibrium as the one with exactly 15 MA
eq15 = eq;


% Test changing the total plasma current with constraints = 1
for ip = [15000000:-100000:7500000 7500000:100000:15000000]
  x.ip = ip;
  y = gseq(x);
end

% Test changing the coil current in PF6
% It's not possible to keep all of ip, li, betap fixed when the plasma limits
% so reconfigure to constraints = 0
config.constraints = 0;

% By reinitializing to eq15, we will know that eq.cplasma is exactly 15 MA
[y, e, eq, eqx] = gseq(x, eq15, config);


for ic6 = [ic6_0 35000:5000:150000 150000:-5000:35000 ic6_0]
  x.cc0t(6) = ic6; 
  y = gseq(x);
end
% At this point the equilibrium is not converged
% This can be seen in the psi error and Ip, li, betap values in the plot


% The flag 'converge' will make gseq converge to a flux error < 1e-6  
% on every call, and not just during initialization
x.converge = true;
y = gseq(x);
% This makes cpasma, li, betap very closely the same as for eq15


% Reconfigure again to ip, li, betap constraints
config.constraints = 1;
gseq([], init, config);


% Make a funny shape
x.li = 0.5;
x.betap = 1;
y = gseq(x);
% If there is a solution that fits inside the limiter then it will be found
% For large plasmas there is at the most 1 solution for a given input vector
% Note that the freaky shape is what you get if the coil currents aren't changed


% Return to the normal shape
x.li = 1;
x.betap = 0.7;
y = gseq(x);


% Reset the converge flag to false for faster execution
x.converge = false;
