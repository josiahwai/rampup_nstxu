% Use gsdesign to determine what the coil currents were for given geqdsk
% file. Compare this to actual coil & vessel currents.

clear; clc; close all

% ========
% SETTINGS
% ========
shot = 204660;
time_ms = 100;
% shotdir = '/u/jwai/rampup_nstxu2/eq/geqdsk/';
shotdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk/';


load('nstxu_obj_config2016_6565.mat')
specify_vessel_currents = 1;

% ====================================
% Define spec,init,config for gsdesign
% ====================================

efit_eq = read_eq(shot, time_ms/1000, shotdir, 'NSTX');
init = efit_eq.gdata;


% NSTXU geometry
fcnturn = tok_data_struct.fcnturn';
Ifrac = zeros(size(fcnturn));
cccirc = [2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 11 12 13] - 1;
for icirc = 1:max(cccirc)
  icoils = find(cccirc == icirc);
  Ifrac(icoils) = fcnturn(icoils) / sum(fcnturn(icoils));
end
init.turnfc = tok_data_struct.fcnturn';
init.fcturn = Ifrac;
init.fcid = cccirc;
init.ecid = [1 1 1 1 1 1 1 1];
init.ecturn = [112 110 109.5 108.5 108.5 109.5 110 112];


% Define coil circuits and limits
spec.cccirc = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 ...
    11 12 13];

spec.limits.ic(1:tok_data_struct.nc,1) = [-20 0 0 -8 0 0 -13 -13 -13 ...
    -13 0 0 0 0 0 0 -24 -24 -24 -24 -13 -13 -13 -13 0 0 -8 0 0]'*1000;

spec.limits.ic(1:tok_data_struct.nc,2) = [20 15 0 13.5 15 15 8 8 8 8 ...
    13 13 13 13 13 13 0 0 0 0 8 8 8 8 15 15 13.5 0 15]'*1000;

% PF1BU/L, PF1CU/L, and PF4 were not used this campaign (shot 204660)
% see tok_data_struct.ccnames coil names/indices
ilock = [3 4 11:16 27 28]; 
spec.locks.ic = nan(size(spec.cccirc));
spec.locks.ic(ilock) = 0;


config = tok_data_struct;
config.max_iterations = 20;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

% Constrain profiles to three degrees of freedom
config.constraints = 1;

% Take profile details from init
config.pres0 = init.pres;
config.fpol0 = init.fpol;

% A generous number of knots ensures the profile can be made in gsdesign
config.psikn = (0:config.nr-1)/(config.nr-1);


% To uniquely specify the equilibrium it will now suffice to specify:
% boundary, boundary flux, total current, li, betap

spec.targets.rsep = init.rbbbs;
spec.targets.zsep = init.zbbbs;
spec.weights.sep = 10 * ones(length(spec.targets.rsep),1);
   
spec.targets.psibry = init.psibry;
spec.weights.psibry = 1;

spec.targets.cpasma = init.cpasma;
spec.weights.cpasma = 1;

gs_configure
gs_initialize
gs_eq_analysis

spec.targets.li = li;
spec.weights.li = 10;

spec.targets.betap = betap;
spec.weights.betap = 10;

if specify_vessel_currents
  coils = load('cc204660_100.mat');
  
  % vacuum vessel degrouping matrix
  vvgroup = [1  1  2  2  3  4  5  5  5  6  6  6  6  6  6  6  6  6  6  6 ...
       6  6  6  6  7  7  7  7  7  8  8  8  8  9  9  9  9  9  9  9  9  9 ...
       9  9  9  9  9  9  9  9  9  9  9  9  9  9  9 10 10 10 11 11 11 11 ...
      11 12 12 13 13 13 14 15 16 17 18 18 18 19 19 20 20 20 20 20 21 21 ...
      21 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 ...
      22 22 22 23 23 23 23 24 24 24 24 24 25 25 25 25 25 25 25 25 25 25 ...
      25 25 25 25 25 26 26 26 27 28 29 29 30 30 31 31 32 32 33 34 35 36 ...
      37 38 39 40];

  Pvv = zeros(length(vvgroup), max(vvgroup));  
  for i = 1:max(vvgroup)
    Pvv(:,i) = (vvgroup==i) / sum(vvgroup==i);
  end
  spec.locks.iv = Pvv*coils.iv_true';
end


eq = gsdesign(spec, init, config);
eq.cc = eq.ic;

[~,i] = unique(spec.cccirc);
ic = eq.ic(i);
coils = load(['cc' num2str(shot) '_' num2str(time_ms) '.mat']);

figure
bar([coils.ic_true' ic])
title('Coil current Comparison', 'fontsize', 18)
legend('True', 'gsdesign reconstructed', 'fontsize', 16)








































