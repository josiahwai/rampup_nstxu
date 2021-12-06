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
dt = mean(diff(t));

% load vacuum geometry, inductances, resistances
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
struct_to_ws(tok_data_struct);
circ = nstxu2016_circ(tok_data_struct);
sys = load([ROOT 'sysid/fit_coils_vessels/coil_vessel_fit.mat']).coil_vessel_fit;
res = load([ROOT 'sysid/fit_plasma_resistance/fits_all/res' num2str(shot) '.mat']).res;


% fetch efits - will use efit profiles and shape as targets
[targets, targets_array, efit01_eqs] = read_target(shot, t, tok_data_struct);
[~,i] = min(abs(targets.time - shape_time));
target = targets_array(i);
response = vacuum_response3(target, tok_data_struct);
eq_target = efit01_eqs.gdata(i);
targets.icx_smooth = targets.icx;
targets.icx_smooth(:,[6 8 9]) = smoothdata(targets.icx(:,[6 8 9]), 'lowess', 70);


% fetch voltages
include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
  'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
vobjcsignals = get_vobjcsignals(shot, t, [], include_coils);
vobjcsignals.sigs = smoothdata(vobjcsignals.sigs, 'lowess', 50);
vsignals = interp1(vobjcsignals.times, vobjcsignals.sigs, t);


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
istart = 1;

%%
%
% clear all; clc; close all
% load('matlab.mat')

%%
% close all
% istart = 10;
% eq = eqs{istart};
% x = xall(istart,:)';

%%
% MAIN CONTROL SIMULATION LOOP
for i = istart:N    
  
  % solve for new equilibrium
  target = targets_array(i);
  x = [target.icx target.ivx target.ip]';  
  profiles.pres = efit01_eqs.gdata(i).pres;
  profiles.fpol = efit01_eqs.gdata(i).fpol;
  
  [spec,init,config] = make_gs_inputs(x, profiles, efit01_eqs.gdata(i), tok_data_struct);  
  % [spec,init,config] = make_gs_inputs(x, profiles, eq, tok_data_struct);
  config.max_iterations = 8;
  
  %   init = efit01_eqs.gdata(i);
  spec.weights.sep = ones(size(spec.weights.sep)) * 10;
%   spec.targets.ic = spec.locks.ic;
%   spec.weights.ic = ones(size(spec.targets.ic)) * 3e-3;
%   spec.locks.ic = nan(size(spec.locks.ic));
  
  
  eq = gsdesign(spec, init, config);
  
  % save stuff
  eqs{i} = eq;
  xall(i,:) = x;
  
  
  figure(2)
  clf
  hold on
  plot(targets.time, targets.icx, '--r')
  plot(t(1:i), xall(1:i, circ.iicx), 'b')
  % plot(t(1:i), xtarg(1:i, circ.iicx), 'g')
  xlim([0 t(i)+0.1])
  
end




























