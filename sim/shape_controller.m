clear all; clc; close all
warning('off','MATLAB:polyshape:repairedBySimplify')
warning('off', 'curvefit:fit:nonDoubleYData')
warning('off', 'stats:pca:ColRankDefX')

ROOT = getenv('RAMPROOT');

% we will try to match shape of this shot at the shape_time
shot = 204660;  
shape_time = 0.7; 
% t = [shape_time shape_time+0.1];
enforce_stability = 0;

t0 = 0.07;
tf = 0.9;
N = 50;
t = linspace(t0, tf, N);

% load vacuum geometry, inductances, resistances
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
struct_to_ws(tok_data_struct);
circ = nstxu2016_circ(tok_data_struct);
sys = load([ROOT 'sysid/fit_coils_vessels/coil_vessel_fit.mat']).coil_vessel_fit;
res = load([ROOT 'sysid/fit_plasma_resistance/fits_all/res' num2str(shot) '.mat']).res;

% fetch efits - will use efit profiles and shape as targets
[targets, targets_array, efit01_eqs] = read_target(shot, t, tok_data_struct);

%%
[~,i] = min(abs(targets.time - shape_time));
target = targets_array(i);


% weights: coil currents, flux errs, grad psi err, coil voltages
% wt.vc = ones(1,circ.ncx_keep) * 1;
% wt.vc(1) = 1e6;  % don't use OH coil for shaping
% wt.ic = ones(1,circ.ncx_keep) * 1;
% wt.psi_cp = ones(size(target.rcp)) * 1;
% wt.xp_psi_r = 1;
% wt.xp_psi_z = 1;

sys = load('coil_vessel_fit.mat').coil_vessel_fit;
response = vacuum_response3(target, tok_data_struct);

% Mxxi = inv(sys.Mxx);
% A = -Mxxi*diag(sys.Rxx);
% B = Mxxi(:,circ.ikeep);
% C = [response.disdis(circ.ikeep,:); response.dpsicpdis; response.dpsibrydis_r; response.dpsibrydis_z];
% D = 0;
% 
% W = diag([wt.ic wt.psi_cp wt.xp_psi_r wt.xp_psi_z]);
% Q = C'*W*C;
% R = diag(wt.vc);
% 
% Klqr = lqr(A,B,Q,R);

%%
eq = efit01_eqs.gdata(i-10);
eqT = efit01_eqs.gdata(i);

y = read_isoflux(eq, target, tok_data_struct);

wt.icx = ones(circ.ncx,1) * 0;
wt.psicp = ones(size(target.rcp(:))) * 1;
wt.psi_r = 1;
wt.psi_z = 1;

target.icx = y.icx;
% target.icx = zeros(circ.ncx,1);
target.psicp_err = 0*target.rcp(:);
target.psi_r = 0;
target.psi_z = 0;

targetvec = [target.icx(circ.ikeep); 
             target.psicp_err; 
             target.psi_r; 
             target.psi_z];

yvec = [y.icx(circ.ikeep); 
        y.psicp - y.target_bdef_psi; 
        y.target_bdef_psi_r; 
        y.target_bdef_psi_z];

C = [response.disdis(circ.ikeep,:); 
     response.dpsicpdis - response.dpsibrydis; 
     response.dpsibrydis_r; response.dpsibrydis_z];

C = C(:, [circ.iicx_keep(1:end)]);
   
W = diag([wt.icx(circ.ikeep); wt.psicp; wt.psi_r; wt.psi_z]);
   
dy_target = targetvec - yvec;
di_target = pinv(W*C)*W*dy_target;


psizr = eq.psizr;
% psizr_app = reshape(mpc*eq.ic + mpv*eq.iv, nz, nr);
% psizr_pla = eq.psizr - psizr_app;

mpcx = mpc * circ.Pcc;

dpsizr_app = mpcx(:,circ.ikeep(1:end)) * di_target;
dpsizr_app = reshape(dpsizr_app, nz, nr);

psizr_new = psizr + dpsizr_app;

[rx,zx,psix] = isoflux_xpFinder(psizr_new,0.6,-1,rg,zg);


figure
plot_eq(eqT)
contour(rg,zg,psizr_new,[psix psix], 'b')
contour(rg,zg,eq.psizr,[eq.psibry eq.psibry], 'g')

%%
% 
% opts.use_vessel_currents = 0;
% target.islimited = 0;
% 
% r = vacuum_response3(target, tok_data_struct);
% 
% if opts.use_vessel_currents
%   i = [circ.iicx circ.iivx];
%   n = circ.ncx + circ.nvx;
% else
%   i = circ.iicx;
%   n = circ.ncx;
% end
% 
% % form the equations A*currents = b
% % b = flux error at control points, flux gradient at x-point
% % A = response of errors to the currents
% 
% A = [eye(n);  
%   r.dpsicpdis(:,i) - r.dpsibrydis(:,i);
%   r.dpsibrydis_r(i);
%   r.dpsibrydis_z(i)];
% 
% [psibry, psibry_r, psibry_z] = bicubicHermite(...
%   rg,zg,eq.psizr,target.rbdef, target.zbdef);
% 
% psi_cp_err = bicubicHermite(rg,zg,eq.psizr,target.rcp, target.zcp)';
% 
% b = -[zeros(n,1);
%   psi_cp_err - psibry;
%   psibry_r;
%   psibry_z];
% 
% % weights
% wt.ic = ones(1,circ.ncx) * 1e-6;
% wt.ic(circ.iicx_remove) = 1e8;
% wt.iv = [];
% if opts.use_vessel_currents
%   wt.iv = ones(1,circ.nvx) * 1e-6;
% end
% wt.cp = ones(1,length(target.rcp));
% if target.islimited
%   wt.x = [0 0];
% else
%   wt.x = [1 1];
% end
% weights = diag([wt.ic wt.iv wt.cp wt.x]);
% 
% % solve weighted least squares problem to obtain currents
% currents = pinv(weights*A)*(weights*b);
% 
% 
% 
% dpsizr_app = mpc * currents;
% dpsizr_app = reshape(dpsizr_app, nz, nr);
% 
% psizr_new = psizr + dpsizr_app;
% 
% [rx,zx,psix] = isoflux_xpFinder(psizr_new,0.6,-1,rg,zg);
% 
% 
% figure
% plot_eq(eqT)
% contour(rg,zg,psizr_new,[psix psix], 'b')
% contour(rg,zg,eq.psizr,[eq.psibry eq.psibry], 'g')











%%
% Use gsdesign to import/fit the equilibrium at t=t0
% For this initial eq, we use higher weights on matching equilibrium flux
% and dont lock the coil currents. Then we use the fit coil currents as the
% starting state vector. 
% efit_eq0 = efit01_eqs.gdata(1);
% 
% x = [efit_eq0.icx; efit_eq0.ivx; efit_eq0.cpasma];
% 
% profiles.pres = efit_eq0.pres;
% profiles.fpol = efit_eq0.fpol;
% 
% [spec,init,config] = make_gs_inputs(x, profiles, efit_eq0, tok_data_struct);
% 
% spec.locks.iv = nan(size(spec.locks.iv));
% spec.locks.ic = nan(size(spec.locks.ic));
% 
% spec.targets.iv = efit_eq0.iv;
% spec.weights.iv(1:length(spec.targets.iv)) = 1e-4;
% 
% spec.targets.ic = efit_eq0.ic;
% spec.weights.ic(1:length(spec.targets.ic)) = 1e-4;
% 
% spec.weights.sep = ones(size(spec.weights.sep)) * 10;
% 
% config.max_iterations = 10;
% 
% eq0 = gsdesign(spec, init, config);
% x0 = pinv(circ.Pxx) * [eq0.ic; eq0.iv; eq0.cpasma];  % initial state
% 
% eq = eq0;
% x = x0;
% 
% % Load dynamics matrices
% 
% 
% % MAIN CONTROL SIMULATION LOOP
% % for i = 1:N
%       
%   i = 1; 
%   y = read_isoflux(eq, target, tok_data_struct);
% 
%   
%   % determine dynamics
%   [Lp, Li, Le, li, mcIp, mvIp] = inductance(eq, tok_data_struct);
%   M = sys.Mxx
%   
%   % error vector: target - actual
%   e.ip = targets.ip(i) - y.ip;
%   e.psicp = y.target_bdef_psi - y.psicp;
%   e.psix_r = 0 - y.target_bdef_psi_r;
%   e.psix_z = 0 - y.target_bdef_psi_z;
%   e.zcur = 0 - y.zcur;
%   
%     
%   % Ip control
%   kp.ip = 1;
%   u.oh = kp.ip * e.ip;
%   
  % Shape control
  
  
  
  
% end




%%
% ref_eq = efit01_eqs.gdata(i);
% y = read_isoflux(ref_eq.psizr, ref_eq.psibry, target, tok_data_struct);
% response = vacuum_response2(target, tok_data_struct);













































