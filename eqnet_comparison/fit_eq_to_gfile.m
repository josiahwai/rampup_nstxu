% Loads NSTX gfile, coil currents etc and fits a toksys equilibrium to the gfile
clear all; clc; close all


shot = 204653;
time = 0.5;
tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';
opts.plotit = 0;
opts.cache_dir = '/Users/jwai/Research/rampup_nstxu/eq/cache_eqs/';
opts.cache_it = 1;

% import gfile
init = fetch_eq_nstxu(shot, time, tree, tokamak, server, opts);

%%
% ====================================
% Define spec,init,config for gsdesign
% ====================================
clear spec config

load('nstxu_obj_config2016_6565.mat')

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

% load circuit connection info
circ = nstxu2016_circ(tok_data_struct);
spec.cccirc = circ.cccirc;

% apply limits
spec.limits.ic = circ.limits.ic;

% lock coil currents
spec.locks.ic = init.ic;
spec.locks.iv = init.iv;
   
% lock Ip
spec.locks.cpasma = init.cpasma;

%%
eq1 = gsdesign(spec, init, config);


%%
% Perform linearizations for each coil



icoil = 1;
dI = 1000;

% perturb coil current
spec2 = spec;
init2 = init;
icx = init.icx;
icx(icoil) = icx(icoil) + dI;
ic = circ.Pcc * icx;
init2.ic = ic;
init2.icx = icx;
spec2.locks.ic = ic;
dic = init2.ic - init.ic;

eq2 = gsdesign(spec2, init2, config);

dpsi_pla = eq2.psizr_pla - eq1.psizr_pla;

figure
contourf(eq1.rg, eq1.zg, dpsi_pla / dI, 20)
colorbar

















