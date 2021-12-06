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
clear all; clc; close all
load('matlab.mat')

%%
% close all
% istart = 10;
% eq = eqs{istart};
% x = xall(istart,:)';

%%
di_err_integral = 0;

% MAIN CONTROL SIMULATION LOOP
for i = istart:N
  
  
  % determine target coil currents
  target = targets_array(i);
  target.icx = targets.icx_smooth(i,:)';
  target.psicp_err = 0*target.rcp(:);
  target.psi_r = 0;
  target.psi_z = 0;
  
  y = read_isoflux(eq, target, tok_data_struct);
  
  targetvec = [target.icx(:); target.psicp_err; target.psi_r; target.psi_z];
  yvec = [y.icx; y.psicp - y.target_bdef_psi; y.target_bdef_psi_r; y.target_bdef_psi_z];
  
  C = [response.disdis(circ.iicx,:); response.dpsicpdis - response.dpsibrydis; response.dpsibrydis_r; response.dpsibrydis_z];
  C = C(:,circ.iicx);
  
  wt.icx = ones(circ.ncx,1) * 1e-6;
  wt.icx([2 5 10 13]) = 1e3;
  wt.icx(1) = 1e3;  
  wt.icx(circ.iicx_remove) = 1e8;  % high weight that unused coils stay at zero
  wt.psicp = ones(size(target.rcp(:))) * 1;
  wt.psicp(1) = 3;
  wt.psi_r = 0;
  wt.psi_z = 0;
  
  W = diag([wt.icx; wt.psicp; wt.psi_r; wt.psi_z]);
  
  dy_target = targetvec - yvec;
  di_target = pinv(W*C)*W*dy_target;
  I_target = y.icx + di_target;
  di_err_integral = di_err_integral + di_target*dt;
  
  %   mpcx = tok_data_struct.mpc * circ.Pcc;
  %   dpsizr_app = reshape(mpcx*di_target, nz, nr);
  %   psizr_new = eq.psizr + dpsizr_app;
  %   [rx,zx,psix] = isoflux_xpFinder(psizr_new,0.6,-1,rg,zg);
  %   figure
  %   plot_eq(eq_target)
  %   contour(rg,zg,psizr_new,[psix psix], 'b')
  %   contour(rg,zg,eq.psizr,[eq.psibry eq.psibry], 'g')
  %   figure
  %   bar([target.icx I_target])
  
  
  % determine dynamics
  [Lp, Li, Le, li, mcIp, mvIp] = inductance(eq, tok_data_struct);
  Rp = interp1(res.t, res.Rp, t(i));
  M = [sys.Mxx [mcIp; mvIp]];
  M = [M; [mcIp' mvIp' Lp]];
  R = diag([sys.rc; sys.rv; Rp]);
  
  Minv = inv(M);
  A = -Minv*R;
  B = Minv(:,circ.iicx);
  A = numerically_stabilize(A,100);
  
  
  % shape control
  %   wt.vc = ones(circ.ncx_keep,1) * 500;
  %   wt.vc(1) = 1e4; % OH coil
  %   wt.ic_tracking = ones(circ.ncx_keep,1) * 1;
  %   Abar = A(circ.ikeep, circ.ikeep);
  %   Bbar = B(circ.ikeep, circ.ikeep);
  %   Q = diag(wt.ic_tracking);
  %   R = diag(wt.vc);
  %   Klqr = lqr(Abar, Bbar, Q, R);
  %   dv = zeros(circ.ncx,1);
  %   dv(circ.ikeep) = Klqr*di_target(circ.ikeep);
  %   dv(circ.ikeep) = sys.rc(circ.ikeep) .* di_target(circ.ikeep);
   
  
  dv = kps .* di_target + kis .* di_err_integral;  
  dv(1) = 0;
  
  v0 = vsignals(i,:)';    
  v = v0 + dv;
  
  
  v(circ.iicx_remove) = 0;
  xdot = A*x + B*v;
  x = x + xdot*dt;
  x(circ.ii_unipolar) = max(x(circ.ii_unipolar), 0);
  x(circ.iicx_remove) = 0;
  
  % solve for new equilibrium
  profiles.pres = efit01_eqs.gdata(i).pres;
  profiles.fpol = efit01_eqs.gdata(i).fpol;
  
  [spec,init,config] = make_gs_inputs(x, profiles, efit01_eqs.gdata(i), tok_data_struct);  
  % [spec,init,config] = make_gs_inputs(x, profiles, eq, tok_data_struct);
  config.max_iterations = 8;
  
  %   init = efit01_eqs.gdata(i);
  spec.weights.sep = ones(size(spec.weights.sep)) * 10;
  spec.targets.ic = spec.locks.ic;
  spec.weights.ic = ones(size(spec.targets.ic)) * 3e-3;
  spec.locks.ic = nan(size(spec.locks.ic));
  
  
  eq = gsdesign(spec, init, config);
  
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
  % plot(t(1:i), xtarg(1:i, circ.iicx), 'g')
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
contour(rg,zg,eq_prev.psizr,[eq_prev.psibry eq_prev.psibry], 'g')
contour(rg,zg,eq_pred.psizr,[eq_pred.psibry eq_pred.psibry], 'b')
contour(rg,zg,eq_targ.psizr,[eq_targ.psibry eq_targ.psibry], '--b')
contour(rg,zg,eq.psizr,[eq.psibry eq.psibry], 'r')
set(gcf, 'Position', [640 115 564 1062])

  
  






%%

close all
clear tf

for i = 1:circ.ncx
  
  lc = sys.lc(i);
  rc = sys.rc(i);
  tau = lc/rc;
  K = 1/rc;
  s = tf('s');
  P = K / (tau*s + 1);
  
  Tr = 0.06;
  w0 = 2.2/Tr;
  zeta = 0.7;
  kp = (2*zeta*w0*tau - 1) / K;
  ki = tau*w0^2 / K;
  C = kp + ki/s;
  kr = rc;
  
  
  u_CL = (C+kr) / (1+P*C);
  I_CL = P*(C+kr) / (1 + P*C);
  
  Iref = 1000;
  
  %   figure
  %   sgtitle(circ.ccnames{i})
  %   subplot(211)
  %   step(Iref * I_CL)
  %   subplot(212)
  %   step(Iref * u_CL)
  %   drawnow
  
  kps(i,1) = kp;
  kis(i,1) = ki;
  krs(i,1) = kr;
  
end
kis(circ.iicx_remove) = 0;
kps(circ.iicx_remove) = 0;
krs(circ.iicx_remove) = 0;


























