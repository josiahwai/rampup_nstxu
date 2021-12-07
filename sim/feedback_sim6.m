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
clear all; clc; close all
load('matlab.mat')


%%
% clear all; clc; close all
% load('matlab2.mat')
% istart = 27;
% eq = eqs{istart};
% x = xall(istart,:)';


%%
di_err_integral = 0;

% MAIN CONTROL SIMULATION LOOP
for i = istart:N
  
  % SHAPE CONTROL
  gap_opts.plotit = 0;   
  gaps = get_nstxu_gaps(efit01_eqs.gdata(i+1), gap_opts);
  
  target = targets_array(i+1);
  target.icx = targets.icx(i+1,:)';
  target.rcp = [gaps(:).r];
  target.zcp = [gaps(:).z];
  target.psicp_err = 0*target.rcp(:);
  target.psi_r = 0;
  target.psi_z = 0;
  target.gap_dist = gaps.dist;
  [target.rx, target.zx] = isoflux_xpFinder(efit01_eqs.gdata(i).psizr, 0.6, -1.1, rg, zg);
  target.r_ingap = 0.35;
  target.z_ingap = 0;
  
  y = read_isoflux(eq,target,tok_data_struct);
    
  gaps = get_nstxu_gaps(eq, gap_opts);
  e = target.gap_dist - gaps.dist;
  e([5 7 9])
  
  targetvec = [target.icx; target.psicp_err; target.psi_r; target.psi_z; 0; 0];
  yvec = [y.icx; y.psicp - y.target_bdef_psi; y.psix_r; y.psix_z; y.psix - y.target_bdef_psi; y.psix - y.psi_in_gap];
  
  response = vacuum_response3(target, tok_data_struct);
  C = [response.disdis(circ.iicx,:); response.dpsicpdis - response.dpsibrydis; 
       response.dpsixdis_r; response.dpsixdis_z; response.dpsixdis - response.dpsibrydis; 
       response.dpsixdis - response.dpsi_ingapdis];     
  C = C(:,circ.iicx);
  
  idx = [1 2 5 10 13];
  
  if t(i) < 0.21
    wt.icx = ones(circ.ncx,1) * 1e-5;
    wt.icx(idx) = 1e8;
    wt.icx(circ.iicx_remove) = 1e8;  
    wt.psicp = ones(size(target.rcp(:))) * 0;
    wt.psicp([5 7 9]) = 1;  
    wt.psix_r = 0;
    wt.psix_z = 0;
    wt.psix_bdef = 0;
    wt.psix_ingap = 0;
  else
    wt.icx = ones(circ.ncx,1) * 1e-5;
    wt.icx(idx) = 1e8;
    wt.icx(circ.iicx_remove) = 1e8; 
    wt.icx([2 13]) = 5e-5; % PF1AU/L    
    wt.psicp = ones(size(target.rcp(:))) * 1;
    wt.psicp([5 7 9]) = 1;      
    wt.psix_r = 0;
    wt.psix_z = 0;
    wt.psix_bdef = 0;
    wt.psix_ingap = 1;
  end      
  
  W = diag([wt.icx; wt.psicp; wt.psix_r; wt.psix_z; wt.psix_bdef; wt.psix_ingap]);
  
  dy_target = targetvec - yvec;
  di_target = pinv(W*C)*W*dy_target;    
  di_target(circ.iicx_remove) = 0;  
  
%   di_max = max(di_target([6 8 9]));  
%   if di_max > 200  % trust region
%     di_target([6 8 9]) = 200/di_max * di_target([6 8 9]);  
%   end
  
  I_target = y.icx + di_target;
  
  
  % determine dynamics and integrate
  [Lp, Li, Le, li, mcIp, mvIp] = inductance(eq, tok_data_struct);
  Rp = interp1(res.t, res.Rp, t(i));
  
  M = [sys.mvv mvIp; mvIp' Lp];
  Minv = inv(M);
  A = -Minv * diag([sys.rv; Rp]);
  B = -Minv * [sys.mcv'; mcIp'];
  A = numerically_stabilize(A,100);
  
  [Ad,Bd] = c2d(A,B,dt);
  
  ic = I_target;  
  icdot = di_target / dt;
  
  ivp = x([circ.iivx circ.iipx]);
  ivp = Ad*ivp + Bd*icdot;
  
  x = [ic; ivp];
  
  
  % solve for new equilibrium
  profiles.pres = efit01_eqs.gdata(i+1).pres;
  profiles.fpol = efit01_eqs.gdata(i+1).fpol;
  
  % [spec,init,config] = make_gs_inputs(x, profiles, efit01_eqs.gdata(i), tok_data_struct);
  [spec,init,config] = make_gs_inputs(x, profiles, eq, tok_data_struct);
  config.max_iterations = 8;
  
  %   init = efit01_eqs.gdata(i);
  spec.weights.sep = ones(size(spec.weights.sep)) * 1;
  spec.targets.ic = spec.locks.ic;
  spec.weights.ic = ones(size(spec.targets.ic)) * 3e-3;
  spec.locks.ic = nan(size(spec.locks.ic));
  
  j = [find(circ.Pcc(:,6)); find(circ.Pcc(:,9))];
  spec.weights.ic(j) = 1e-3;
  
  
  %   spec.targets.iv = spec.locks.iv;
  %   spec.weights.iv = ones(size(spec.targets.iv)) * 3e-3;
  %   spec.locks.iv = nan(size(spec.locks.iv));
  %   spec.weights.fpol(1:length(spec.targets.fpol)) = 1e-4;
  %   spec.weights.pres(1:length(spec.weights.pres)) = 1e-4;
  spec.weights.fpol(1:length(spec.targets.fpol)) = 100;
  spec.weights.pres(1:length(spec.weights.pres)) = 1e-1;
  
  spec.locks.cpasma = efit01_eqs.gdata(i).cpasma;
  
%   spec.targets.li = target.li;
%   spec.weights.li = 10;
  
  
  pcurrt = efit01_eqs.gdata(i).pcurrt;
  zcur = sum(sum(pcurrt.*zgg)) / sum(pcurrt(:));
  spec.targets.zcur = zcur;
  spec.weights.zcur = 100;
      
%   spec.targets.cpasma = spec.locks.cpasma;
%   spec.weights.cpasma = 1;
%   spec.locks.cpasma = nan;
  
  eq = gsdesign(spec, init, config);
  
  x = pinv(circ.Pxx) * [eq.ic; eq.iv; eq.cpasma];
  
  % save stuff
  eqs{i} = eq;
  xall(i,:) = x;
  yall(i,:) = y;
  xtarg(i,:) = I_target;
  
  
  figure(2)
  clf
  hold on
  plot(targets.time, targets.icx, '--r')
  plot(t(1:i), xall(1:i, circ.iicx), 'b')
  plot(t(1:i), xtarg(1:i, circ.iicx), 'g')
  xlim([0 t(i)+0.1])
  for j = 1:circ.ncx
    text(t(i)+0.001, xall(i,j), circ.ccnames{j})
  end
  
end


%%
x = pinv(circ.Pxx) * [eq.ic; eq.iv; eq.cpasma];
if i == 1
  eq_prev = eq0;
  dx = x - x0;
else
  eq_prev = eqs{i-1};
  dx = x - xall(i-1,:)';
end
mpcx = tok_data_struct.mpc * circ.Pcc;
mpvx = tok_data_struct.mpv * circ.Pvv;
dcphidip = eq_prev.pcurrt(:) / sum(eq_prev.pcurrt(:));
mp_ip = mpp_full * dcphidip;
dpsizr_app = [mpcx mpvx mp_ip] * dx;
dpsizr_app = reshape(dpsizr_app, nz, nr);
psizr_pred = eq_prev.psizr + dpsizr_app;
eq_pred = eq_params(psizr_pred, tok_data_struct, 0);


dpsizr_app = reshape(mpcx * di_target, nz, nr);
psizr_targ = eq_prev.psizr + dpsizr_app;
eq_targ = eq_params(psizr_targ, tok_data_struct, 0);


figure
plot_nstxu_geo(tok_data_struct)
scatter(target.rcp, target.zcp, 'k', 'filled')
contour(rg,zg,eq_prev.psizr,[eq_prev.psibry eq_prev.psibry], 'g')   % Previous
contour(rg,zg,eq_pred.psizr,[eq_pred.psibry eq_pred.psibry], 'b')   % Expected based on: vacuum model + actual coil currents
contour(rg,zg,eq_targ.psizr,[eq_targ.psibry eq_targ.psibry], '--b')  % Expected based on: vacuum model + target coil currents
contour(rg,zg,eq.psizr,[eq.psibry eq.psibry], 'r')  % Actual fitted
set(gcf, 'Position', [640 115 564 1062])



































