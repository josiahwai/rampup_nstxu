clear all; clc; close all
warning('off','MATLAB:polyshape:repairedBySimplify')
warning('off', 'curvefit:fit:nonDoubleYData')
warning('off', 'stats:pca:ColRankDefX')

ROOT = getenv('RAMPROOT');

% we will try to match shape of this shot at the shape_time
shot = 204660;  
shape_time = 0.7; 
enforce_stability = 0;

t0 = 0.07;
tf = 0.9;
N = 200;
t = linspace(t0, tf, N);

% load vacuum geometry, inductances, resistances
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);
sys = load([ROOT 'sysid/fit_coils_vessels/coil_vessel_fit.mat']).coil_vessel_fit;
res = load([ROOT 'sysid/fit_plasma_resistance/fits_all/res' num2str(shot) '.mat']).res;


% fetch efits - will use efit profiles and shape as targets
[targets, targets_array, efit01_eqs] = read_target(shot, t, tok_data_struct);
[~,i] = min(abs(targets.time - shape_time));
target = targets_array(i);


% Use gsdesign to import/fit the equilibrium at t=t0
% For this initial eq, we use higher weights on matching equilibrium flux
% and dont lock the coil currents. Then we use the fit coil currents as the
% starting state vector. 
efit_eq0 = efit01_eqs.gdata(1);

x = [efit_eq0.icx; efit_eq0.ivx; efit_eq0.cpasma];

profiles.pres = efit_eq0.pres;
profiles.fpol = efit_eq0.fpol;

[spec,init,config] = make_gs_inputs(x, profiles, efit_eq0, tok_data_struct);

spec.locks.iv = nan(size(spec.locks.iv));
spec.locks.ic = nan(size(spec.locks.ic));

spec.targets.iv = efit_eq0.iv;
spec.weights.iv(1:length(spec.targets.iv)) = 1e-4;

spec.targets.ic = efit_eq0.ic;
spec.weights.ic(1:length(spec.targets.ic)) = 1e-4;

spec.weights.sep = ones(size(spec.weights.sep)) * 10;

config.max_iterations = 10;

eq0 = gsdesign(spec, init, config);
x0 = pinv(circ.Pxx) * [eq0.ic; eq0.iv; eq0.cpasma];  % initial state

eq = eq0;
x = x0;

% Load dynamics matrices


% MAIN CONTROL SIMULATION LOOP
% for i = 1:N
      
  i = 1; 
  y = read_isoflux(eq, target, tok_data_struct);

  
  % determine dynamics
  [Lp, Li, Le, li, mcIp, mvIp] = inductance(eq, tok_data_struct);
  M = sys.Mxx
  
  % error vector: target - actual
  e.ip = targets.ip(i) - y.ip;
  e.psicp = y.target_bdef_psi - y.psicp;
  e.psix_r = 0 - y.target_bdef_psi_r;
  e.psix_z = 0 - y.target_bdef_psi_z;
  e.zcur = 0 - y.zcur;
  
    
  % Ip control
  kp.ip = 1;
  u.oh = kp.ip * e.ip;
  
  % Shape control
  
  
  
  
% end




%%
% ref_eq = efit01_eqs.gdata(i);
% y = read_isoflux(ref_eq.psizr, ref_eq.psibry, target, tok_data_struct);
% response = vacuum_response2(target, tok_data_struct);













































