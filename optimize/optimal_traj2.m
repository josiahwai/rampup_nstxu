
ccc
RAMPROOT = getenv('RAMPROOT');

% ========
% Settings 
% ========
shot = 204660;
constraints.t = [0 0.1 0.4 0.86 0.88]; 
constraints.t = 0:0.01:0.9;
constraints.n = length(constraints.t);

coil_opts.plotit = 0;
coil_constraints = fetch_coilcurrents_nstxu(shot, constraints.t, coil_opts);
% coils_true = fetch_coilcurrents_nstxu(shot, t);
coils_true = coil_constraints;

enforce_stability = 0;

vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);

init = fetch_coilcurrents_nstxu(shot, 0, coil_opts);
ic0 = init.icx;
iv0 = init.ivx;
ip0 = init.ip;
x0 = [ic0; iv0; ip0];

% simulation timing
Ts = .01;
tstart = 0;
tend = 0.9;
tsample = tstart:Ts:tend;
N = length(tsample);
t = tsample;

wt.icx = ones(1,circ.ncx) * 1;
wt.ivx = ones(1,circ.nvx) * 0;
wt.ip = ones(1,circ.np) * 1e-1;

wt.dicx = ones(1,circ.ncx) / Ts^2 * 0;
wt.divx = ones(1,circ.nvx) / Ts^2 * 0;
wt.dip = ones(1,circ.np)  / Ts^2 * 0;
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
%   Mpc = traj.Mpc(i,:); 
%   Mpv = traj.Mpv(i,:);
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
  
%   Ad = inv(eye(size(Ac_vp)) - Ts*Ac_vp);
%   Bd = Ad * Bc_vp * Ts;
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
  
  
  
%   F  = F  + kron(diag(ones(N-i+1,1), -i+1), Apow*Bd);
%   F_ext = F_ext + kron(diag(ones(N-i+1,1), -i+1), Apow*Bd_ext);
%   Apow  = Ad*Apow;
%   E = [E; Apow];
end
F = F / Ts;

%%
% ==============================
% Form the weights and cost fun
% ==============================

% velocity conversion matrices (see derivations): dxhat = Sx*xhat - xprev,
% dichat = Sic*ichat - icprev, ichat = Picx * xhat
Svp = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.nvx+circ.np));
Sc = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.ncx));   

ivp0 = [iv0; ip0];
ic_prev = [ic0; zeros((N-1)*circ.ncx,1)];
ivp_prev = [ivp0; zeros((N-1)*(circ.nvx+1),1)];

u_ext = zeros(N,1);

R = {};
Rv = {};
for i = 1:N
  Q{i} = diag([wt.ivx wt.ip]);
  Qv{i} = diag([wt.divx wt.dip]);
  R{i} = diag(wt.icx);
  Rv{i} = diag(wt.dicx);
end

Qhat = blkdiag(Q{:});
Qvhat = blkdiag(Qv{:});
Rhat = blkdiag(R{:});
Rvhat = blkdiag(Rv{:});

targs.t = tsample;
targs.ic = interp1(coil_constraints.times, coil_constraints.icx', targs.t, 'pchip', 'extrap');
targs.iv = interp1(coil_constraints.times, coil_constraints.ivx', targs.t, 'pchip', 'extrap');
targs.ip = interp1(coil_constraints.times, coil_constraints.ip, targs.t, 'pchip', 'extrap')';
rvp_hat = reshape([targs.iv targs.ip]', [], 1);
rc_hat = reshape(targs.ic', [], 1);

% costfun
Qbar = Qhat + Svp'*Qvhat*Svp;

H = Rhat + Sc'*Rvhat*Sc + Sc'*F'*Qbar*F*Sc;
f1 = -( (rvp_hat'*Qhat + ivp_prev'*Qvhat*Svp)*F*Sc + rc_hat'*Rhat + ic_prev'*Rvhat*Sc);
z = -F*ic_prev + E*[iv0;ip0] + F_ext*u_ext;
f2 = z'*Qbar*F*Sc;
f = (f1 + f2)';

% constraint on matching snapshots
clear Aeq beq
ieq = 1;

P_constraints = zeros(constraints.n, N);
for i = 1:length(constraints.t)
  [~,j] = min(abs(constraints.t(i) - tsample));
  P_constraints(i,j) = 1;
end
Aeq{ieq} = kron(P_constraints, eye(circ.ncx));
beq{ieq} = reshape(coil_constraints.icx, [], 1);


% Solve quadratic program
npv = size(H,1);

Aineq = zeros(0,npv);
bineq = zeros(0,1);

Aeq = cat(1,Aeq{:});
beq = cat(1,beq{:});

[u,s,v] = svd(Aeq);
s = diag(s);
s(s < sqrt(eps)) = [];
n = length(s);
u = u(:,1:n);
v = v(:,1:n);
Aeq = diag(s)*v';
beq = u'*beq;

% Aeq = Aineq;
% beq = bineq;

opts = mpcActiveSetOptions;

iA0 = false(length(bineq), 1);

[ic_hat,exitflag,iA,lambda] = mpcActiveSetSolver(H, f, Aineq, bineq, Aeq, beq, iA0, opts);

ivp_hat = E*ivp0 + F*(Sc*ic_hat - ic_prev); % F*Sc*ic_hat + z;
ivp_hat = reshape(ivp_hat, [], N);
ic_hat = reshape(ic_hat,[],N);
xhat = [ic_hat; ivp_hat];    




figure
hold on
plot(nan,nan,'-r')
plot(nan,nan, '--k')
scatter(nan,nan, 'k', 'filled')
plot(t,xhat(circ.iicx,:), 'r')
plot(coils_true.times, coils_true.icx, '--k')
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
plot(coils_true.times, coils_true.ip, '-k', 'linewidth', 1.5)
box on
legend('Simulated', 'Simulation Target', 'True', 'fontsize', 14)
title(['Ip ' num2str(shot)], 'fontsize', 14, 'fontweight', 'bold')
xlabel('Time [s]', 'fontweight', 'bold', 'fontsize', 14)
ylabel('Current [A]', 'fontweight', 'bold', 'fontsize', 14)







%%

% w = [iv0; ip0];
% 
% % ic = [ic0 xhat(circ.iicx,:)];
% % uhat = diff(ic, 1, 2);
% uhat = gradient(xhat(circ.iicx,:));
% for i = 1:N
%   wsim(:,i) = w;
%   u = uhat(:,i);
%   w = Ad_list{i}*w + Bd_list{i} / Ts * u;
% %   wsim(:,i) = w;
% end
% 
% figure
% subplot(211)
% plot(t,wsim(1:circ.nvx,:))
% subplot(212)
% plot(t,wsim(end,:))

% close all
% for k = 1:circ.nvx
%   figure
%   hold on
%   plot(targs.t, targs.iv, 'Color', [1 1 1]*0.8)
%   plot(t, xhat(circ.iivx(k),:), 'b', 'linewidth', 2)
%   plot(targs.t, targs.iv(:,k), '--r', 'linewidth', 2)
%   title(circ.vvnames{k}, 'fontsize', 18)
% end
  









































































