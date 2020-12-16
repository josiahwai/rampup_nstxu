%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gsdesign_demo_east_LSNF
%
%  PURPOSE: DEMO of gsdesign showing design of lower-snowflake for EAST
%
%  INPUTS:  none
%
%  OUTPUTS: eq, a lower-snowflake equilibrium
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%
%  WRITTEN BY:  Anders Welander  ON	3/12/19
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load EAST tokamak configuration
if exist('east_obj_filename','var')
  try
    clear tok_data_struct
    load(east_obj_filename)
    tok_data_struct.rg;
  catch
    error('east_obj_filename does not hold name of a matfile containing tok_data_struct')
  end
else
  east_obj_filename = [getenv('GATOOLS_ROOT'), ...
    '/tokamaks/east/make/east_obj_2014_6565.mat'];
  if exist(east_obj_filename,'file')
    load(east_obj_filename)
  else
    east_obj_filename = ...
      '/m/GAtools/tokamaks/east/make/east_obj_2014_6565.mat';
    if exist(east_obj_filename,'file')
      load(east_obj_filename)
    else
      disp('Could not find east objects. Please set the variable:')
      disp('east_obj_filename')
      disp('to name of matfile containing tok_data_struct for east')
    end
  end
end

config = tok_data_struct;
config.constraints = 1;
config.psikn = [0 0.40 0.70 1];
config.no_edge_current = true;
config.no_edge_gradient = true;
config.plot_settings.SOL.n = 9;
config.plot_settings.SOL.d = 1e-3;
init = [];

clear spec gsdesign

% Specify snowflake
spec.targets.rsnf =  1.49;
spec.targets.zsnf = -0.85;
spec.targets.tsnf = 20;
spec.weights.snf = 2000;

% Specify bdef point
spec.targets.rbdef = spec.targets.rsnf;
spec.targets.zbdef = spec.targets.zsnf;

% Specify where the flux should equal the boundary flux
spec.targets.rsep = ...
  [ 1.49, 1.47, 2.31, 1.65, 2.20, 1.90, 1.45, 2.302, 2.302, ...
    1.75, 1.50, 1.80, 1.55, 1.60, 1.46, 1.44, 1.700];
spec.targets.zsep = ...
  [-0.85, 0.00, 0.00, 0.82, 0.40,-0.67, 0.40, 0.100,-0.100, ...
   -0.76, 0.70, 0.80, 0.77, 0.80, 0.59, 0.20, 0.825];
spec.weights.sep = ...
  [  500, 1.00, 10.0, 99.0, 1.00, 1.00, 1.00, 10.00, 10.00, ...
    1.00, 10.0, 10.0, 10.0, 10.0, 10.0, 00.0, 10.00];

% Finally, specify plasma current, li, betap, boundary flux, connections
spec.targets.cpasma = 250e3;
spec.weights.cpasma = 1e-3;
spec.targets.li = 1;
spec.weights.li = 10;
spec.targets.betap = 1;
spec.weights.betap = 5;
spec.targets.psibry = 0;
spec.weights.psibry = 5;
spec.cccirc =  [1 2 3 4 5 6 7 8 7 8 9 10 11 12 13 -13];
spec.locks.ic = nan(16,1);
spec.locks.ic(15:16) = 0;

% Call gsdesign with these specs
eq = gsdesign(spec, init, config);
