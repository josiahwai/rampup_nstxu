clear all; clc; close all
clear gsdesign spec init config
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

% ================
% define mpc costs 
% ================
saveit = 0;
savedir = '/Users/jwai/Research/rampup_nstxu/sim/';
modeldir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/original/';
shot = 204660;
t_snapshots = [60:10:300] / 1000;

t0 = t_snapshots(1);
tf = t_snapshots(end);
N = length(t_snapshots) - 1;
tspan = linspace(t0,tf,N+1);
dt = mean(diff(tspan));

nx = 49;
nu = 8;
iic = 1:8;
iiv = 9:48;
iip = 49;

wt.ic = ones(1,8) * 100;
wt.iv = ones(1,40) * 1;
wt.ip = 1e-4;
wt.u  = ones(1,8) * 1e-3;
wt.du = 1 / (dt^2) * ones(1,8) * 3e-3;

% velocity conversion matrix: du = Sv*u
Sv = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(nu));

% weighting matrices
Q = diag([wt.ic wt.iv wt.ip]);   
Qf = diag([wt.ic wt.iv wt.ip]);
Qhat = blkdiag(kron(eye(N-1),Q), Qf);  

R = diag(wt.u); 
Rv = Sv' * kron(eye(N), diag(wt.du)) * Sv;
Rhat = kron(eye(N), R) + Rv; 

% load shot trajectory
[traj, eqs] = make_traj(shot, t_snapshots, tspan, modeldir);
x_t0 = traj.x(1,:)';

% future state targets
rhat = reshape(traj.x', [], 1);
rhat(1:nx) = [];

% prediction model
Chat = kron(eye(N), eye(nx));
Apow  = eye(nx);
F  = zeros(N*nx, N*nu);
OB  = [];
for i = 1:N
  A = squeeze(traj.A(i,:,:));
  B = squeeze(traj.B(i,:,:));
  F  = F  + kron(diag(ones(N-i+1,1), -i+1), Apow*B);
  Apow  = A*Apow;
  OB  = [OB; Apow];
end

H = F'*Chat'*Qhat*Chat*F + Rhat;
if max(max(H-H')) > 1e-10, H = (H+H')/2; end

f = F' * Chat' * Qhat * (Chat * OB * x_t0 - rhat);

Uhat = -H\f;
Uhat = reshape(Uhat, nu, []);  

[x_all, x_tf] = state_dynamics(traj.A, traj.B, x_t0, N, Uhat);

if saveit
  sim_inputs.Uhat = Uhat;
  sim_inputs.tspan = tspan;
  sim_inputs.x0 = x_t0;
  sim_inputs.eq0 = eqs{1};
  sim_inputs.traj = traj;
  sim_inputs.x_all = x_all;
  save([savedir 'sim_inputs204660.mat'], 'sim_inputs');
end

%%
% =============
% Plot results
% =============
figure
sgtitle('Shot 204660')
subplot(221)
title('Coil currents', 'fontsize', 14, 'fontweight', 'bold')
icoil = 1:8;
hold on
plot(tspan, x_all(icoil,:), 'b')
plot(tspan, traj.x(:,icoil)', '--r')
ylabel('A')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

subplot(222)
title('Vessel currents', 'fontsize', 14, 'fontweight', 'bold')
icoil = 9:48;
hold on
plot(tspan, x_all(icoil,:), 'b')
plot(tspan, traj.x(:,icoil)', '--r')
ylabel('A')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

subplot(223)
title('Ip', 'fontsize', 14, 'fontweight', 'bold')
hold on
plot(tspan, x_all(end,:), 'b')
plot(tspan, traj.x(:,end)', '--r')
ylabel('A')
xlabel('Time')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

subplot(224)
hold on
title('Coil voltages', 'fontsize', 14, 'fontweight', 'bold')
stairs(tspan(1:end-1), Uhat')
ylabel('Voltage')
xlabel('Time')

set(gcf,'Position',[28 275 669 505])



figure
hold on
i = 11;
ivess = 9:48;
plot(tspan, traj.x(:,ivess), 'color', [1 1 1] * 0.8)
plot(tspan, traj.x(:, ivess(i)), '--r')
plot(tspan, x_all(ivess(i),:), 'b')


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
function [x_all, x_tf] = state_dynamics(Alist, Blist, x_t0, N, u_all)
  x = x_t0;
  x_all = zeros(length(x_t0), N+1);
  x_all(:,1) = x_t0;
  
  for i = 1:N
    A = squeeze(Alist(i,:,:));
    B = squeeze(Blist(i,:,:));    
    x = A*x + B*u_all(:,i);
    x_all(:,i+1) = x;      
  end
  
  x_tf = x_all(:,end);
end

% =========================================
% Load eq snapshots to make shot trajectory
% =========================================
function [traj, eqs] = make_traj(shot, t_snapshots, tspan, modeldir)
  dt = mean(diff(tspan));
  
  for i = 1:length(t_snapshots)
    
    t = t_snapshots(i) * 1000;
    load([modeldir num2str(shot) '_' num2str(t) '_sys.mat']);
    eq = load(['eq' num2str(shot) '_' num2str(t) '.mat']);
    eq = eq.eq;    
    eqs{i} = eq;
    
    [Ai, Bi] = c2d(sys.As, sys.B, dt);      
    A(i,:,:) = Ai;
    B(i,:,:) = Bi;
        
    psibry(i) = eq.psibry;
    li(i) = eq.li;
    betap(i) = eq.betap;
    pres(i,:) = eq.pres;
    fpol(i,:) = eq.fpol;
    zcur(i) = eq.zcur;
    rcur(i) = eq.rcur;
    
    % Load coil currents
    load(['coils' num2str(shot) '.mat'])
    
    ccnames = {'OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', ...
      'PF5', 'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL'};
    
    icoil = [];
    for k = 1:length(sys.inputs)
      icoil(k) = find(strcmp(ccnames, sys.inputs{k}));
    end
    
    itime = find( floor(coils.t*1000) == t);
    ic(i,:) = coils.ic(itime, icoil);
    iv(i,:) = coils.iv(itime,:);
    ip(i,:) = eq.cpasma;        
  end
  
  % finer interpolation
  method = 'pchip';
  traj.A      = interp1(t_snapshots, A,      tspan, 'previous');
  traj.B      = interp1(t_snapshots, B,      tspan, 'previous');
  traj.psibry = interp1(t_snapshots, psibry, tspan, method);
  traj.li     = interp1(t_snapshots, li,     tspan, method);
  traj.betap  = interp1(t_snapshots, betap,  tspan, method);
  traj.pres   = interp1(t_snapshots, pres,   tspan, method);
  traj.fpol   = interp1(t_snapshots, fpol,   tspan, method);
  traj.zcur   = interp1(t_snapshots, zcur,   tspan, method);
  traj.rcur   = interp1(t_snapshots, rcur,   tspan, method);
  traj.ic     = interp1(t_snapshots, ic,     tspan, method);
  traj.iv     = interp1(t_snapshots, iv,     tspan, method);
  traj.ip     = interp1(t_snapshots, ip,     tspan, method);
  traj.x      = [traj.ic traj.iv traj.ip'];
end








































