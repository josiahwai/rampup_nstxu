clear all; clc; close all
clear gsdesign spec init config
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
RAMPROOT = getenv('RAMPROOT');

% ========
% SETTINGS
% ========
mpc_timebase = [270:10:940] / 1e3;
shot = 204660;
saveit = 0;
savedir = [RAMPROOT '/dev/'];
savefn = 'sim_results_postfit.mat';
modeldir = [RAMPROOT '/buildmodel/built_models/std/'];
eqdir = [RAMPROOT '/eq/geqdsk_import/'];
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;

grey_sys = load('grey_sys_270_940.mat').grey_sys;
fittedAB = grey_sys.UserData.fittedAB;

% ===================
% define mpc weights
% ===================
t0 = mpc_timebase(1);
tf = mpc_timebase(end);
N = length(mpc_timebase) - 1;
tspan = linspace(t0,tf,N+1);
dt = mean(diff(tspan));

circ = nstxu2016_circ(tok_data_struct);

wt.icx = ones(1,circ.ncx) * 1e4; % outputs and states
wt.ivx = ones(1,circ.nvx) * 1e-4;
wt.ip = 1e2;
wt.dicx = ones(1,circ.ncx) / dt^2 * 0;   % velocity of outputs and states
wt.divx = ones(1,circ.nvx) / dt^2 * 0;
wt.dip  = 1 / dt^2 * 0;
wt.u  = ones(1,circ.nu) * 1e-5;  % input voltage
wt.du = ones(1,circ.nu) / (dt^2) * 3e2;  % velocity of input voltage

% velocity conversion matrices (see derivations): duhat = Su*uhat - uprev,
%    dxhat = Sx*xhat - xprev, delta_xhat = Tx*dxhat
Su = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.nu));
Sx = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.nx));

for i = 1:N
  Q{i} = diag([wt.icx wt.ivx wt.ip]);
  Qv{i} = diag([wt.dicx wt.divx wt.dip]);
  R{i} = diag(wt.u);
  Rv{i} = diag(wt.du);
end
Qhat = blkdiag(Q{:});
Qvhat = blkdiag(Qv{:});
Rhat = blkdiag(R{:});
Rvhat = blkdiag(Rv{:});
  

% ================
% solve mpc
% ================
traj = load('traj_270_930.mat').traj;
x0 = traj.x(1,:)';
ytarg = pinv(circ.Pxx_keep) * grey_sys.UserData.ytarg;
ps_voltages = grey_sys.UserData.ps_voltages;
u0 = ps_voltages(1,:)';
rhat = reshape(ytarg(:,2:end), [], 1);

% prediction model
C = eye(circ.nx);
Chat = kron(eye(N), C);
Apow  = eye(circ.nx);
F  = zeros(N*circ.nx, N*circ.nu);
E  = [];
for i = 1:N
  A = squeeze(fittedAB.A(i,:,:));
  B = squeeze(fittedAB.B(i,:,:));  
  F  = F  + kron(diag(ones(N-i+1,1), -i+1), Apow*B);
  Apow  = A*Apow;
  E  = [E; Apow];
end

y0 = C * x0;
y0hat = repmat(y0, N, 1);
x0hat = repmat(x0, N, 1);

xprev = [x0; zeros((N-1)*circ.nx, 1)];
uprev = [u0; zeros((N-1)*circ.nu, 1)];

% Write costfun J to depend only on the decision variable, u. See derivation.
% J=u'*H*u + 2*f'*u

Rbar = Su'*Rvhat*Su + Rhat;
fu = -Su'*Rvhat*uprev;
Hu = Rbar;

Qbar = Chat'*Qhat*Chat + Sx'*Chat'*Qvhat*Chat*Sx;
fx = -Sx'*Chat'*Qvhat*Chat*xprev + Chat'*Qhat*(y0hat - rhat - Chat*x0hat);
fy = F'*Qbar*E*x0 + F'*fx;
Hy = F'*Qbar*F;
  
H = Hu + Hy;
f = fu + fy;

Uhat = -H\f;
Uhat = reshape(Uhat, circ.nu, []);


%%
% ==================
% Integrate solution
% ==================
% Uhat = ps_voltages(2:end,:)';

x = x0;
xall = zeros(circ.nx, N+1);
xall(:,1) = x0;

for i = 1:N
  A = squeeze(fittedAB.A(i,:,:));
  B = squeeze(fittedAB.B(i,:,:));
  u = Uhat(:,i);
  x = A*x + B*u;
  xall(:,i+1) = x;
end

%%
% ======================
% Save and Plot results
% ======================
if saveit
  sim_inputs.Uhat = Uhat;
  sim_inputs.tspan = tspan;
  sim_inputs.x0 = x_t0;
  sim_inputs.eq0 = eqs{1};
  sim_inputs.traj = traj;
  sim_inputs.x_all = xall;
  save([savedir savefn], 'sim_inputs');
end


figure
sgtitle('Shot 204660')
ax(1) = subplot(221);
title('Coil currents', 'fontsize', 14, 'fontweight', 'bold')
hold on
plot(tspan, xall(circ.iicx,:), 'b')
plot(tspan, traj.x(:,circ.iicx)', '--r')
ylabel('A')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

ax(2) = subplot(222);
title('Vessel currents', 'fontsize', 14, 'fontweight', 'bold')
hold on
plot(tspan, xall(circ.iivx,:), 'b')
plot(tspan, traj.x(:,circ.iivx)', '--r')
ylabel('A')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

ax(3) = subplot(223);
title('Ip', 'fontsize', 14, 'fontweight', 'bold')
hold on
plot(tspan, xall(circ.iipx,:), 'b')
plot(tspan, traj.x(:,circ.iipx)', '--r')
ylabel('A')
xlabel('Time')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

ax(4) = subplot(224);
hold on
title('Coil voltages', 'fontsize', 14, 'fontweight', 'bold')
stairs(tspan(1:end-1), Uhat')
ylabel('Voltage')
xlabel('Time')

linkaxes(ax,'x')
set(gcf,'Position',[28 275 669 505])

figure
hold on
i = 11;
plot(tspan, traj.x(:,circ.iivx), 'color', [1 1 1] * 0.8)
plot(tspan, traj.x(:, circ.iivx(i)), '--r')
plot(tspan, xall(circ.iivx(i),:), 'b')











































