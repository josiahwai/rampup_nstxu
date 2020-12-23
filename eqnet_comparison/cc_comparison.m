% Compare equilibrium using true coil currents vs. equilibrium from coil currents
% predicted by eqnet/optimization

clear; clc; close all

shot = 204660;
time_ms = 120;
% shotdir = '/u/jwai/rampup_nstxu2/eq/geqdsk/';
shotdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk/';
load('nstxu_obj_config2016_6565.mat')

% ============================
% DEFINE SPEC, INIT, CONFIG
% ============================

efit_eq = read_eq(shot, time_ms/1000, shotdir);
init = efit_eq.gdata;


config = tok_data_struct;
config.max_iterations = 20;
config.constraints = 1;
config.nkn = 4;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

spec.targets.rsep = init.rbbbs;
spec.targets.zsep = init.zbbbs;
spec.weights.sep = 0.1*ones(length(spec.targets.rsep),1);
   
spec.targets.cpasma = init.cpasma;
spec.weights.cpasma = 1;

spec.targets.psibry = init.psibry;
spec.weights.psibry = 1;

config.constraints = 1;  % allow for scaling/peaking of profiles
config.pres0 = init.pres;
config.fpol0 = init.fpol;


gs_configure
gs_initialize
gs_eq_analysis
spec.targets.li = li;
spec.targets.betap = betap;
spec.weights.li = 10;
spec.weights.betap = 10;


spec.cccirc = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 ...
    11 12 13];

spec.limits.ic(1:tok_data_struct.nc,1) = [-20 0 0 -8 0 0 -13 -13 -13 ...
    -13 0 0 0 0 0 0 -24 -24 -24 -24 -13 -13 -13 -13 0 0 -8 0 0]'*1000;

spec.limits.ic(1:tok_data_struct.nc,2) = [20 15 0 13.5 15 15 8 8 8 8 ...
    13 13 13 13 13 13 0 0 0 0 8 8 8 8 15 15 13.5 0 15]'*1000;



%%
% ============================
% MAP COIL CURRENTS TO TOKSYS
% ============================
load('nstxu_obj_config2016_6565.mat');
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

% spec.locks.ic = coils.ic_true(spec.cccirc);
% eq0 = gsdesign(spec, init, config);

spec.locks.ic = coils.ic_pred(spec.cccirc);
eq0 = gsdesign(spec, init, config);



% =======
% PLOT IT
% =======
figure
hold on
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;
% plot_nstxu_geo(tok_data_struct)

[~,cs] = contour(rg,zg,init.psizr,30,'--r', 'linewidth', 0.5);
contour(rg,zg,eq0.psizr,cs.LevelList,'k', 'linewidth', 0.5);
set(gcf,'Position',[204 38 312 533])

















