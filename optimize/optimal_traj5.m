ccc
RAMPROOT = getenv('RAMPROOT');

% ========
% Settings 
% ========
shot = 204660;
enforce_stability = 0;

% simulation timing
Ts = .01;
tstart = 0;
tend = 0.4;
tsample = tstart:Ts:tend;
N = length(tsample);
t = tsample;


% constraints.t = tstart:0.01:0.9;
constraints.t = [tstart 0.4];
constraints.n = length(constraints.t);

coil_opts.plotit = 0;
coil_constraints = fetch_coilcurrents_nstxu(shot, constraints.t, coil_opts);

coil_targs = fetch_coilcurrents_nstxu(shot, t);

vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;

tok_data_struct = vac_sys.build_inputs.tok_data_struct;
tok_data_struct.imks = 1;

circ = nstxu2016_circ(tok_data_struct);

init = fetch_coilcurrents_nstxu(shot, tstart, coil_opts);
ic0 = init.icx;
iv0 = init.ivx;
ip0 = init.ip;

wt.icx = ones(1,circ.ncx) * 1e-1;
wt.ivx = ones(1,circ.nvx) * 1e-6;
wt.ip = ones(1,circ.np) * 1e-4;
wt.cp = ones(1,10) * 1;

wt.dicx = ones(1,circ.ncx) / Ts^2 * 1e-5;
wt.divx = ones(1,circ.nvx) / Ts^2 * 1e-5;
wt.dip = ones(1,circ.np)  / Ts^2 * 0;
wt.dcp = ones(1,10) / Ts^2 * 0;

% ==============================
% Load time dependent parameters
% ==============================

% in the future, should iteratively solve for the time-dep params
% for now, just load.

% Load fitted vacuum model parameters
Mxx = vac_sys.sysid_fits.Mxx;
Rxx = vac_sys.sysid_fits.Rxx;
Mvv = Mxx(circ.iivx, circ.iivx);
Mcc = Mxx(circ.iicx, circ.iicx);
Mvc = Mxx(circ.iivx, circ.iicx);
Mcv = Mvc';
Rv = Rxx(circ.iivx);
Rc = Rxx(circ.iicx);

% load plasma mutuals, Rp, Lp
modeldir = [RAMPROOT '/buildmodel/built_models/mcc/'];
traj.t = tsample;
[traj.Mpc, traj.Mpv] = load_Mpc_Mpv(modeldir, tsample);

load('Rp_ipfit2.mat');
traj.Rp_t = interp1(Rp_ipfit.t, Rp_ipfit.value, t, 'linear', 'extrap');

fit_Lp = load('fit_Lp.mat').fit_Lp;
traj.Lp_t = fit_Lp(t);

% form the state-space dynamics matrices
npv = circ.nvx + circ.np;
Apow  = eye(npv);
F  = [];
F_row = zeros(npv, N*circ.ncx);
F_ext = [];
F_row_ext = zeros(npv, N);
E  = [];

for i = 1:N  
  Rp = traj.Rp_t(i);
  Lp = traj.Lp_t(i);
  
  Mpc = squeeze( interp1(traj.t, traj.Mpc, t(i), 'linear', 'extrap'));
  Mpv = squeeze( interp1(traj.t, traj.Mpv, t(i), 'linear', 'extrap'));
  Mvp = Mpv';   
  
  M = [Mvv Mvp; Mpv Lp];
  R = diag([Rv; Rp]);
  M_inv = inv(M);  
  Ac_vp = -M_inv * R;
  Bc_vp = -M_inv * [Mvc; Mpc];
  
  if enforce_stability
    Ac_vp = numerically_stabilize(Ac_vp, 1e3);
  end
  
  [Ad,Bd] = c2d(Ac_vp, Bc_vp, Ts);
  
  Bd_ext = zeros(circ.nvx+circ.np, 1);
  Bd_ext(end) = 1;
  
  Ad_list{i} = Ad;
  Bd_list{i} = Bd;
  Bd_ext_list{i} = Bd_ext;
  
  idx = (circ.ncx*(i-1)+1):(circ.ncx*i);
  F_row = Ad * F_row;
  F_row(:,idx) = Bd;
  F = [F; F_row];
  
  Apow = Ad*Apow;
  E = [E; Apow];
    
  F_row_ext = Ad * F_row_ext;
  F_row_ext(:,i) = Bd_ext;
  F_ext = [F_ext; F_row_ext];
  
end
F = F / Ts;

%%
% =========================
% Form the prediction model
% =========================

% fetch eq
eq0 = fetch_eq_nstxu(shot, tstart);
eq_ref = fetch_eq_nstxu(shot, 0.4);

ibad = (eq_ref.rbbbs == 0 & eq_ref.zbbbs == 0);
eq_ref.rbbbs(ibad) = [];
eq_ref.zbbbs(ibad) = [];

targ_geo.cp.n = 10;
[targ_geo.cp.r, targ_geo.cp.z] = interparc(eq_ref.rbbbs, eq_ref.zbbbs, targ_geo.cp.n, 1, 0);


% measure current state, y := [ic iv ip cp_err]
ny = circ.nx + targ_geo.cp.n;

cp_err0 = bicubicHermite(eq0.rg, eq0.zg, eq0.psizr, targ_geo.cp.r, targ_geo.cp.z) - eq0.psibry;

y0 = [ic0; iv0; ip0; cp_err0];
x0 = [iv0; ip0];

xprev = [x0; zeros((N-1)*(circ.nvx+1),1)];
ic_prev = [ic0; zeros((N-1)*circ.ncx,1)];
yprev = [y0; zeros((N-1)*ny, 1)];

x0hat = repmat(x0, N, 1);
y0hat = repmat(y0, N, 1);
ic0hat = repmat(ic0, N, 1);

% velocity conversion matrices
Svp = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.nvx+circ.np));
Sc = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.ncx));   
Sy = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(ny));   

% get response model
cdata = build_cmat_data(eq_ref, circ, tok_data_struct, targ_geo);

DC = [cdata.x; cdata.dpsicpdix];

D = DC(:, 1:circ.ncx);
C = DC(:, circ.ncx+1:end);
Chat = kron(eye(N), C);
Dhat = kron(eye(N), D);

xk = x0;

z = y0hat + Chat * (E*xk - F*ic_prev - x0hat) - Dhat*ic0hat;
M = Chat*F*Sc + Dhat;

%%
% weights and targets
for i = 1:N
  Q{i} = diag([wt.icx wt.ivx wt.ip wt.cp]);
  Qv{i} = diag([wt.dicx wt.divx wt.dip wt.dcp]);
end

Qhat = blkdiag(Q{:});
Qvhat = blkdiag(Qv{:});
Qbar = Qhat + Sy'*Qvhat*Sy;

targs.t = tsample;
targs.ic = interp1(coil_targs.times, coil_targs.icx', targs.t, 'pchip', 'extrap');
targs.iv = interp1(coil_targs.times, coil_targs.ivx', targs.t, 'pchip', 'extrap');
targs.ip = interp1(coil_targs.times, coil_targs.ip, targs.t, 'pchip', 'extrap')';
targs.cp_err = zeros(N, targ_geo.cp.n);

% targs_ic = ones(N,1) * targs.ic(1,:);
% targs_iv = targs.iv * 0;
% targs_ic = interp1(coil_targs.times([1 end]), coil_targs.icx(:,[1 end])', targs.t, 'pchip', 'extrap');
% targs_ip = interp1(coil_targs.times([1 end]), coil_targs.ip([1 end]), targs.t, 'pchip', 'extrap')';


rhat = reshape([targs.ic targs.iv targs.ip targs.cp_err]', [], 1);

% rhat = reshape([targs.ic targs.iv targs.ip targs.cp_err]', [], 1);


% cost function
H = M'*Qbar*M;
ft = (z'*Qbar - rhat'*Qhat - yprev'*Qvhat*Sy) * M;

npv = size(H,1);


% Equality constraints:
% ---------------------
clear Aeq beq
ieq = 1;

% coil currents at the specified times
P_constraints = zeros(constraints.n, N);
for i = 1:length(constraints.t)
  [~,j] = min(abs(constraints.t(i) - tsample));
  P_constraints(i,j) = 1;
end
Aeq{ieq} = kron(P_constraints, eye(circ.ncx));
beq{ieq} = reshape(coil_constraints.icx, [], 1);


% some coils turned off -- constrain to zero
ieq = ieq + 1;
noff = length(circ.iicx_remove);
off = zeros(noff, circ.ncx);
for i = 1:noff
  off(i,circ.iicx_remove(i)) = 1;
end
Aeq{ieq} = kron(eye(N), off);
beq{ieq} = zeros(N*noff,1);


Aeq = cat(1,Aeq{:});
beq = cat(1,beq{:});

% enforce linear independence
[u,s,v] = svd_tol(Aeq, 0.9999);
Aeq = u'*Aeq;
beq = u'*beq;



% Inequality constraints:
% ----------------------
Aineq = zeros(0,npv);
bineq = zeros(0,1);


% Solve quadratic program
opts = mpcActiveSetOptions;
iA0 = false(length(bineq), 1);

[ic_hat,exitflag,iA,lambda] = mpcActiveSetSolver(H, ft', Aineq, bineq, Aeq, beq, iA0, opts);


% unpack solution
yhat = M*ic_hat + z;

yhat = reshape(yhat,[],N);
ic_hat = reshape(ic_hat,[],N);

xhat = yhat(circ.iisx, :);

cphat = yhat(circ.nx+1:end,:);


%%
close all

figure
hold on
plot(nan,nan,'-r')
plot(nan,nan, '--k')
scatter(nan,nan, 'k', 'filled')
plot(t,xhat(circ.iicx,:), 'r')
plot(coil_targs.times, coil_targs.icx, '--k')
for i = 1:length(coil_constraints.times)-1
  scatter(ones(circ.ncx,1) * coil_constraints.times(i), coil_constraints.icx(:,i), 'k', 'filled')
end
% plot(targs.t, targs.ic, '--r')
xlabel('Time [s]', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Current [A]', 'fontweight', 'bold', 'fontsize', 14)
title(['OH/PF Coil Currents ' num2str(shot)], 'fontsize', 14, 'fontweight', 'bold')
legend('Trajectory', 'True currents', 'Specified targets', 'fontsize', 14)
box on

figure
hold on
plot(t, xhat(circ.iivx,:), 'b')
plot(targs.t, targs.iv, '--r')

figure 
hold on
plot(t, xhat(circ.iipx,:), 'b', 'linewidth', 1.5)
plot(targs.t, targs.ip, '-r', 'linewidth', 1.5)
plot(coil_targs.times, coil_targs.ip, '--k', 'linewidth', 1.5)
box on
legend('Simulated', 'Simulation Target', 'True', 'fontsize', 14)
title(['Ip ' num2str(shot)], 'fontsize', 14, 'fontweight', 'bold')
xlabel('Time [s]', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Current [A]', 'fontweight', 'bold', 'fontsize', 14)

figure
plot(t, cphat)






%%
x0 = xhat(:,1);
xf = xhat(:,end);


[spec, init, config] = make_gsdesign_inputs2(xf, tok_data_struct, eq_ref, circ);
% [spec, init, config] = make_gsdesign_inputs2(xf, tok_data_struct, eq0, circ);

eq = gsdesign(spec,init,config)


figure
plot_eq(eq)






















































