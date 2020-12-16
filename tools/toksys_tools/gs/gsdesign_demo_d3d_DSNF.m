%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gsdesign_demo_d3d_DSNF
%
%  PURPOSE: DEMO of gsdesign showing design of double-snowflake for DIII-D
%
%  INPUTS:  none
%
%  OUTPUTS: eq, a double snowflake equilibrium
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%
%  WRITTEN BY:  Anders Welander  ON	3/12/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load D3D tokamak configuration
if exist('d3d_obj_filename','var')
  try
    clear tok_data_struct
    load(d3d_obj_filename)
    tok_data_struct.rg;
  catch
    error('d3d_obj_filename does not hold name of a matfile containing tok_data_struct')
  end
else
  d3d_obj_filename = [getenv('GATOOLS_ROOT'), ...
    '/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_6565.mat'];
  if exist(d3d_obj_filename,'file')
    load(d3d_obj_filename)
  else
    d3d_obj_filename = ...
      '/m/GAtools/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_6565.mat';
    if exist(d3d_obj_filename,'file')
      load(d3d_obj_filename)
    else
      disp('Could not find d3d objects. Please set the variable:')
      disp('d3d_obj_filename')
      disp('to name of matfile containing tok_data_struct for d3d')
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

% Specify snowflakes
spec.targets.rsnf = 1.25*[+1 +1];
spec.targets.zsnf = 1.00*[-1 +1];
spec.targets.tsnf = [-120 120];
spec.weights.snf = 1*[1 1];

% Specify bdef point
spec.targets.rbdef = spec.targets.rsnf(1);
spec.targets.zbdef = spec.targets.zsnf(1);
spec.limits.bdef_dpsibar = 0;

% Specify where the flux should equal the boundary flux
r1 = [1.1, 2.15, 1.1,  1.1, 1.11,  1.11, 1.9,  1.9];
z1 = [0.0, 0.00, 0.1, -0.1, 0.20, -0.20, 0.5, -0.5];
r2 = [1.65,  1.65, 2.0,  2.0, 2.1,  2.1];
z2 = [0.75, -0.75, 0.4, -0.4, 0.2, -0.2];
r3 = [1.13,  1.13, 1.15 , 1.15, 1.165,  1.165, 1.18,  1.18, 1.19,  1.19];
z3 = [0.30, -0.30, 0.40, -0.40, 0.500, -0.500, 0.60, -0.60, 0.70, -0.70];
spec.targets.rsep = [spec.targets.rsnf r1 r2 r3];
spec.targets.zsep = [spec.targets.zsnf z1 z2 z3];

% Finally, specify plasma current, li, betap, boundary flux, connections
spec.targets.cpasma = 1e6;
spec.weights.cpasma = 1e-3;
spec.targets.li = 1;
spec.weights.li = 100;
spec.targets.betap = 1;
spec.weights.betap = 100;
spec.targets.psibry = 0;
spec.weights.psibry = 100;
spec.buscode = [0 0 1 1 1 1 1 0 0 1 1 1 1 1 1 1 0 0 1 1];

% Call gsdesign with these specs
eq = gsdesign(spec, init, config);
