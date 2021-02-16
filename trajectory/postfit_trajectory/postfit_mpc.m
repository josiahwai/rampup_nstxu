clear all; clc; close all
clear gsdesign spec init config
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
RAMPROOT = getenv('RAMPROOT');


% ============================
% LOAD PARAMETERS FROM FITTING
% ============================
% load('grey_sys2.mat')
% rcc = grey_sys.Parameters(1).Value;
% rvv = grey_sys.Parameters(2).Value;  
% rp_t = grey_sys.Parameters(3).Value;
% lp_t = grey_sys.Parameters(4).Value;
% voltage_scale = grey_sys.Parameters(5).Value;
% fitted_resistances = variables2struct(rcc, rvv, rp_t, lp_t, voltage_scale);

fittedAB = load('fittedAB_270_930.mat').fittedAB;


% ========
% Settings 
% ========
saveit = 0;
savedir = [RAMPROOT '/dev/'];
savefn = 'sim_results_postfit.mat';
modeldir = [RAMPROOT '/buildmodel/built_models/std/'];
eqdir = [RAMPROOT '/eq/geqdsk_import/'];
shot = 204660;
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;

mpc_timebase = [270:10:940] / 1e3;

if length(mpc_timebase) ~= size(fittedAB.A,1)
  error('Fitted parameters does not match simulation timing')
end

% ================
% define mpc costs 
% ================
t0 = mpc_timebase(1);
tf = mpc_timebase(end);
N = length(mpc_timebase) - 1;
tspan = linspace(t0,tf,N+1);
dt = mean(diff(tspan));

circ = nstxu2016_circ(tok_data_struct);

wt.icx = zeros(1,circ.ncx);
wt.icx(circ.ikeep) = 1e4;
wt.ivx = ones(1,circ.nvx) * 1e-4;
wt.ip = 10000;
wt.u  = ones(1,circ.nu) * 1e-5;
wt.du = ones(1,circ.nu) / (dt^2) * 0;

% velocity conversion matrix: du = Sv*u
Sv = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.nu));

% weighting matrices
for i = 1:N
  Q{i} = diag([wt.icx wt.ivx wt.ip]);
  R{i} = diag(wt.u);
  Rv{i} = diag(wt.du);
end
Qhat = blkdiag(Q{:});
Rhat = blkdiag(R{:}) + Sv' * blkdiag(Rv{:}) * Sv;


% load shot trajectory
% [traj, eqs] = load_interp_trajectory(shot, circ, mpc_timebase, tspan, eqdir, modeldir);
% [traj, eqs] = load_interp_trajectory(shot, circ, mpc_timebase, eqdir, modeldir);
traj = load('traj_270_930.mat').traj;

load('eqs.mat')
x_t0 = traj.x(1,:)';

% future state targets

rhat = reshape(traj.x', [], 1);
rhat(1:circ.nx) = [];

% prediction model
Chat = kron(eye(N), eye(circ.nx));
Apow  = eye(circ.nx);
F  = zeros(N*circ.nx, N*circ.nu);
OB  = [];
for i = 1:N
  A = squeeze(fittedAB.A(i,:,:));
  B = squeeze(fittedAB.B(i,:,:));  
  F  = F  + kron(diag(ones(N-i+1,1), -i+1), Apow*B);
  Apow  = A*Apow;
  OB  = [OB; Apow];
end

H = F'*Chat'*Qhat*Chat*F + Rhat;
if max(max(H-H')) > 1e-10, H = (H+H')/2; end

f = F' * Chat' * Qhat * (Chat * OB * x_t0 - rhat);

% Uhat = -H\f;
% Uhat = reshape(Uhat, circ.nu, []);  

%%
coils = load('coils_greybox.mat').coils;
Uhat = zeros(circ.nu, N);
grey_sys = load('grey_sys_270_930.mat').grey_sys;
% voltage_scale = grey_sys.Parameters(5).Value(circ.ikeep);
% v = interp1(coils.t, coils.v', mpc_timebase(2:end), 'next')';
v = grey_sys.UserData.ps_voltages';
Uhat = v(:,2:N+1);
% Uhat(circ.ikeep,:) = v;
x_t0 = load('x.mat').x;

[x_all, x_tf] = state_dynamics(fittedAB, x_t0, N, Uhat);

if saveit
  sim_inputs.Uhat = Uhat;
  sim_inputs.tspan = tspan;
  sim_inputs.x0 = x_t0;
  sim_inputs.eq0 = eqs{1};
  sim_inputs.traj = traj;
  sim_inputs.x_all = x_all;
  save([savedir savefn], 'sim_inputs');
end


% =============
% Plot results
% =============
figure
sgtitle('Shot 204660')
ax(1) = subplot(221);
title('Coil currents', 'fontsize', 14, 'fontweight', 'bold')
hold on
plot(tspan, x_all(circ.iicx,:), 'b')
plot(tspan, traj.x(:,circ.iicx)', '--r')
ylabel('A')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

ax(2) = subplot(222);
title('Vessel currents', 'fontsize', 14, 'fontweight', 'bold')
hold on
plot(tspan, x_all(circ.iivx,:), 'b')
plot(tspan, traj.x(:,circ.iivx)', '--r')
ylabel('A')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

ax(3) = subplot(223);
title('Ip', 'fontsize', 14, 'fontweight', 'bold')
hold on
plot(tspan, x_all(circ.iipx,:), 'b')
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
plot(tspan, x_all(circ.iivx(i),:), 'b')

% figure
% sgtitle('204660', 'fontweight', 'bold', 'fontsize', 18)
% subplot(211)
% plot(tspan, traj.li)
% ylabel('Li', 'fontweight', 'bold', 'fontsize', 16)
% subplot(212)
% plot(tspan, traj.betap)
% ylabel('Betap', 'fontweight', 'bold', 'fontsize', 14)
% xlabel('Time [s]')

% ==============
% State dynamics
% ==============
function [x_all, x_tf] = state_dynamics(fittedAB, x_t0, N, u_all)
  x = x_t0;
  x_all = zeros(length(x_t0), N+1);
  x_all(:,1) = x_t0;
  
  for i = 1:N
    A = squeeze(fittedAB.A(i,:,:));
    B = squeeze(fittedAB.B(i,:,:));
    x = A*x + B*u_all(:,i);
    x_all(:,i+1) = x;      
  end
  
  x_tf = x_all(:,end);
end





































