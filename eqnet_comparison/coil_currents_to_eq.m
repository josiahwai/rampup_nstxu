function eq = coil_currents_to_eq(init, ic, iv, ip)

load('nstxu_obj_config2016_6565.mat')

spec.cccirc = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 ...
    11 12 13];

spec.limits.ic(1:tok_data_struct.nc,1) = [-20 0 0 -8 0 0 -13 -13 -13 ...
    -13 0 0 0 0 0 0 -24 -24 -24 -24 -13 -13 -13 -13 0 0 -8 0 0]'*1000;

spec.limits.ic(1:tok_data_struct.nc,2) = [20 15 0 13.5 15 15 8 8 8 8 ...
    13 13 13 13 13 13 0 0 0 0 8 8 8 8 15 15 13.5 0 15]'*1000;

config = tok_data_struct;
config.max_iterations = 20;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

% Constrain profiles to three degrees of freedom
config.constraints = 1;

% Take profile details from init
% config.pres0 = init.pres * .1389 / .1746;
config.pres0 = init.pres;
config.fpol0 = init.fpol;

% A generous number of knots ensures the profile can be made in gsdesign
config.psikn = (0:config.nr-1)/(config.nr-1);

% To uniquely specify the equilibrium it will now suffice to specify:
% boundary, boundary flux, total current, li, betap

spec.targets.rsep = init.rbbbs;
spec.targets.zsep = init.zbbbs;
spec.weights.sep = 1e-3 * ones(length(spec.targets.rsep),1);
   
spec.targets.psibry = init.psibry;
spec.weights.psibry = 1e-3;

% gs_configure
% gs_initialize
% gs_eq_analysis

% spec.targets.li = init.li;
spec.targets.li = .5854;
spec.weights.li = 1;

% spec.targets.betap = init.betap;
spec.targets.betap = .1389;
spec.weights.betap = 1;

% Lock coil, vessel, plasma currents
spec.locks.cpasma = ip;

spec.locks.ic = ic(spec.cccirc);

% PF1BU/L, PF1CU/L, and PF4 were not used this campaign (shot 204660)
% see tok_data_struct.ccnames coil names/indices
ilock = [3 4 11:16 27 28]; 
spec.locks.ic(ilock) = 0;

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
spec.locks.iv = Pvv*iv;

eq = gsdesign(spec, init, config);









































