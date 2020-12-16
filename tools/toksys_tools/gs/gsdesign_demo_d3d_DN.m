%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gsdesign_demo_d3d_DN
%
%  PURPOSE: DEMO of gsdesign showing design of EAST-like double-x-point for DIII-D
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
if exist('tok_data_struct','var')
  if ~strcmp(upper(tok_data_struct.tokamak),'D3D')
    disp('Warning: tok_data_struct should contain description of D3D')
  end
else % try to load tok_data_struct, don't want to bother user if possible since this is a demo
  if ~exist('gatools_root','var')
    gatools_root = getenv('GATOOLS_ROOT');
  end
  tokdir = [gatools_root,'/tokamaks/d3d/make/with_helicon2017/'];
  if ~exist(tokdir,'dir')
    tokdir = '/fusion/pillar-archive/data/m/GAtools_sccs/tokamaks/d3d/make/with_helicon2017/';
  end
  if exist([tokdir,'d3d_obj_6565.mat'],'file')
    load([tokdir,'d3d_obj_6565.mat'])
    disp(['Loaded tok_data_struct from: ',tokdir,'d3d_obj_6565.mat'])
  elseif exist([tokdir,'d3d_obj_tmp_mks.mat'],'file')
    load ([tokdir,'d3d_obj_tmp_mks.mat'])
    disp(['Loaded tok_data_struct from: ',tokdir,'d3d_obj_tmp_mks.mat'])
  else
    error('No tok_data_struct in workspace, please load one for D3D')
  end
end

config = regrid(65,65,tok_data_struct);
config.constraints = 1;
config.psikn = [0 0.40 0.70 1];
config.plot_settings.SOL.n = 9;
config.plot_settings.SOL.d = 1e-3;
buscode = [0 0 1 1 1 1 1 0 0 1 1 1 1 1 1 1 0 0 1 1];
ind = find(buscode);
Pcc = [1 zeros(1,config.nc-2); eye(config.nc-1)]; % icc(1) is current in both EcoilA and B
Pcc(ind(1),ind-1) = -1;
Pcc = Pcc(:,[1:ind(1)-2 ind(1):config.nc-1]);
config.Pcc = Pcc;

init = [];

clear spec gsdesign

% EAST shape
eqeast.rbbbs = [1.3981 1.3976 1.4002 1.4057 1.4138 1.4187 1.4244 1.4376 1.4536 1.4625 1.4724 1.4942 1.5062 1.5193 1.5475 1.55 1.5787 1.5938 1.608 1.6375 1.6812 1.725 1.7319 1.7687 1.8125 1.8522 1.8562 1.9 1.9437 1.9515 1.9875 2.0312 2.0339 2.075 2.1025 2.1187 2.1592 2.1625 2.2057 2.2062 2.2429 2.25 2.2715 2.2915 2.2937 2.3037 2.308 2.3044 2.2937 2.2929 2.2735 2.25 2.2457 2.2093 2.2062 2.1635 2.1625 2.1187 2.1075 2.075 2.0396 2.0312 1.9875 1.9577 1.9437 1.9 1.8582 1.8562 1.8125 1.7687 1.7363 1.725 1.6812 1.6375 1.608 1.5938 1.5801 1.55 1.5496 1.5211 1.5062 1.495 1.4716 1.4625 1.4511 1.4337 1.4197 1.4187 1.409 1.4018 1.3981];
eqeast.zbbbs = [-0.075 0 0.075 0.15 0.225 0.2622 0.3 0.375 0.45 0.487 0.525 0.6 0.63718 0.675 0.75 0.75619 0.825 0.85917 0.88999 0.87559 0.85293 0.82891 0.825 0.80368 0.77656 0.75 0.74721 0.71583 0.68143 0.675 0.64406 0.60266 0.6 0.55692 0.525 0.50498 0.45 0.44519 0.375 0.37403 0.3 0.28323 0.225 0.15 0.13909 0.075 0 -0.075 -0.14559 -0.15 -0.225 -0.28964 -0.3 -0.375 -0.3805 -0.45 -0.45146 -0.51102 -0.525 -0.56279 -0.6 -0.60838 -0.64942 -0.675 -0.6865 -0.72044 -0.75 -0.75136 -0.78019 -0.80668 -0.825 -0.83123 -0.85444 -0.87618 -0.88997 -0.85716 -0.825 -0.75105 -0.75 -0.675 -0.63337 -0.6 -0.525 -0.49297 -0.45 -0.375 -0.3 -0.29433 -0.225 -0.15 -0.075];

% Transform EAST shape into DIII-D
scale_factor = 1.25;
r_offset = -0.7;
z_offset = 0;

% Specify points rsep, zsep where flux should equal the boundary flux
spec.targets.rsep = r_offset+scale_factor*eqeast.rbbbs;
spec.targets.zsep = z_offset+scale_factor*eqeast.zbbbs;

% Specify points rx, zx where the poloidal field should vanish
[~, ix1] = min(spec.targets.zsep);
[~, ix2] = max(spec.targets.zsep);
spec.targets.rx = spec.targets.rsep([ix1 ix2]);
spec.targets.zx = spec.targets.zsep([ix1 ix2]);
spec.weights.x = 1000*[1 1];

% The DN can be more balanced with higher weights on the x-points
spec.weights.sep = ones(1,length(spec.targets.rsep));
spec.weights.sep(ix1) = 1;
spec.weights.sep(ix2) = 1;
spec.locks.rsep = spec.targets.rx;
spec.locks.zsep = spec.targets.zx;

% Finally, specify plasma current, li, betap, boundary flux, connections
spec.targets.cpasma = 250e3;
spec.weights.cpasma = 1;
spec.targets.li = 1;
spec.weights.li = 10;
spec.targets.betap = 1;
spec.weights.betap = 10;
spec.targets.psibry = 0;
spec.weights.psibry = 10;

% Call gsdesign with these specs
eq = gsdesign(spec, init, config);
