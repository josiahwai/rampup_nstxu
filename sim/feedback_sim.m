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


% MAIN CONTROL SIMULATION LOOP
for i = istart:N
         
  y = read_isoflux(eq, target, tok_data_struct);

  % determine target coil currents
  target.icx = targets.icx(i,:)'; 
  target.psicp_err = 0*target.rcp(:);
  target.psi_r = 0;
  target.psi_z = 0;
  
  targetvec = [target.icx; target.psicp_err; target.psi_r; target.psi_z];
  yvec = [y.icx; y.psicp - y.target_bdef_psi; y.target_bdef_psi_r; y.target_bdef_psi_z];
  
  C = [response.disdis(circ.iicx,:); response.dpsicpdis - ones(12,1)*response.dpsibrydis; response.dpsibrydis_r; response.dpsibrydis_z];
  C = C(:,circ.iicx);
  
  wt.icx = ones(circ.ncx,1) * 1e-4;
  wt.icx(1) = 1e-4;  % OH coil
  if t(i) < 0.35, wt.icx([5 10]) = 1e3; end  
  wt.icx(circ.iicx_remove) = 1e8;  % high weight that unused coils stay at zero
  wt.psicp = ones(size(target.rcp(:))) * 1;
  wt.psi_r = 0.1;
  wt.psi_z = 0.1;
  
  W = diag([wt.icx; wt.psicp; wt.psi_r; wt.psi_z]);
  
  dy_target = targetvec - yvec;
  di_target = pinv(W*C)*W*dy_target;
  I_target = y.icx + di_target;
  
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
     
%     dv(circ.ikeep) = sys.rc(circ.ikeep) .* di_target(circ.ikeep);
    
  % Ip control
%   kp = 1e-2;
%   dv(1) = dv(1) + kp * (targets.ip(i) - y.ip);
    
  
  % integrate
  v0 = vsignals(i,:)';  
  
  kp = 3e-2;
  dv = kp * (targets.icx_smooth(i,:)' - y.icx); 
  
  v = v0 + dv;
  
%   v0 = sys.rc .* y.icx; 
%   v0(circ.iicx_remove) = 0;  
%   v0(1) = vsignals(i,1);
%   
%   v = v0 + dv;
  


%   I_target(1) = targets.icx(i,1);  
%   v = 1 * sys.rc .* (I_target - y.icx);
%   v(circ.iicx_remove) = 0;
  
  xdot = A*x + B*v;
  x = x + xdot*dt;
  x(circ.ii_unipolar) = max(x(circ.ii_unipolar), 0);
  
  % solve for new equilibrium
  profiles.pres = efit01_eqs.gdata(i).pres;
  profiles.fpol = efit01_eqs.gdata(i).fpol;
  [spec,init,config] = make_gs_inputs(x, profiles, eq, tok_data_struct);
  
%   init = efit01_eqs.gdata(i);
  spec.weights.sep = ones(size(spec.weights.sep)) * 1;
  spec.targets.ic = spec.locks.ic;
  spec.weights.ic = ones(size(spec.targets.ic)) * 3e-3;
  spec.locks.ic = nan(size(spec.locks.ic));
  
  
  eq = gsdesign(spec, init, config);
    
  % save stuff
  eqs{i} = eq;
  xall(i,:) = x;
  yall(i,:) = y;
  
  
  figure(2)
  clf
  hold on
  plot(targets.time, targets.icx, '--r')
  plot(t(1:i), xall(1:i, circ.iicx), 'b')
  xlim([0 t(i)+0.1])
  
end




figure
hold on
plot(t(1:size(xall,1)), xall(:,circ.iicx), 'b')
% plot(t(1:size(xall,1)), xall(:,[2 13]), 'g')
plot(targets.time, targets.icx, '--r')
xlim([0 t(i)+0.1])











































