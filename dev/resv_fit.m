ccc
shot = 204660;
model_dir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/resv_fit/';

% Initialize Rv and Rp(t)
times = 60:10:200;
load('fit_Rp')
Rp0 = fit_Rp(times/1e3);

load('nstxu_sys.mat')
load('Rvv_fit.mat')
Rv0 = Rvv_fit;
% R0 = diag(nstxu_sys.Pxx' * nstxu_sys.rxx * nstxu_sys.Pxx);
% Rv0 = R0(end-40:end-1);

load('/Users/jwai/Research/rampup_nstxu/dev/sim_inputs204660.mat')
traj = sim_inputs.traj;

dt = mean(diff(times));
tt = 1:length(times);
ivpdot_sim = gradient([traj.iv(tt,:) traj.ip(tt)']', dt);
ivpdot_true = gradient(sim_inputs.x_all(end-40:end,tt), dt);


e = ivpdot_sim(:) - ivpdot_true(:);

dedr = [];
for k = 1:length(times)
  Rp = Rp0(k);
  Rv = Rv0;
  % sys = build_nstxu_system_fun(times(k), Rp, Rv, 0, model_dir);
  load([model_dir num2str(shot) '_' num2str(times(k)) '_sys.mat']);
  lstari = inv(sys.lstar);
  
  Ivp = [traj.iv(k,:) traj.ip(k)]'; % vessel/plasma current at timestep
  Fvpvp = -lstari(end-40:end, end-40:end) * diag(Ivp);  
  v = zeros(size(times));
  v(k) = 1;
  mask = blkdiag(eye(40), v);
  dedr = [dedr; Fvpvp*mask];  % change in error wrt vessel/plasma resistance  
end

r = [Rv0; Rp0];
H = dedr'*dedr;
f = dedr'*e;
dr = quadprog(H,f,[],[]);

%%
r = max(r+dr,0);
Rp = r(41:end);
Rv = r(1:40);

for k = 1:length(times)
  sys = build_nstxu_system_fun(times(k), Rp(k), Rv, 1, model_dir);
end
mpc_trajectory

























