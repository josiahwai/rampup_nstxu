clear all; clc; close all
clear gsdesign spec init config
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

% ==============
% Define problem
% ==============
shot = 204660;
t_snapshots = [60:10:300] / 1000;

t0 = t_snapshots(1);
tf = t_snapshots(end);
N = 50;
tspan = linspace(t0,tf,N+1);
dt = mean(diff(tspan));

[traj, eqs] = make_traj(shot, t_snapshots, tspan);
x_t0 = traj.x(1,:)';

nx = 54;
nu = 13;
iic = 1:13;
iiv = 14:53;
iip = 54;

wt.ic = ones(1,13) * 100;
wt.iv = ones(1,40) * 1e-4;
wt.ip = 1e-2;
wt.u  = ones(1,13) * 1e-3;

wt.ic([3 4 7 11 12]) = 0;

Q = diag([wt.ic wt.iv wt.ip]);
Qf = diag([wt.ic wt.iv wt.ip]);
Qhat = blkdiag(kron(eye(N-1),Q), Qf);

R = diag(wt.u);
Rhat = kron(eye(N),R);

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

% future targets
rhat = reshape(traj.x', [], 1);
rhat(1:nx) = [];

f = F' * Chat' * Qhat * (Chat * OB * x_t0 - rhat);

Uhat = -H\f;
Uhat = reshape(Uhat, nu, []);  

[x_all, x_tf] = state_dynamics(traj.A, traj.B, x_t0, N, Uhat);

%%
icoil = 1:13;
% icoil = 14:53;
figure
subplot(311)
hold on
plot(tspan, x_all(icoil,:), 'b')
plot(tspan, traj.x(:,icoil)', '--r')
subplot(312)
hold on
plot(tspan, x_all(end,:), 'b')
plot(tspan, traj.x(:,end)', '--r')
subplot(313)
plot(tspan(1:end-1), Uhat)
set(gcf,'Position',[287 418 560 420])

% ==============
% State dynamics
% ==============
function [x_all, x_tf] = state_dynamics(Alist, Blist, x_t0, N, u_all)
  x = x_t0;
  x_all = zeros(length(x_t0), N+1);
  x_all(:,1) = x_t0;
  
  for i = 1:N
    A = squeeze(Alist(i,:,:));
    max(abs(eig(A)))
    B = squeeze(Blist(i,:,:));    
    x = A*x + B*u_all(:,i);
    x_all(:,i+1) = x;      
  end
  
  x_tf = x_all(:,end);
end








function [traj, eqs] = make_traj(shot, t_snapshots, tspan)
  dt = mean(diff(tspan));
  
  for i = 1:length(t_snapshots)
    t = t_snapshots(i) * 1000;
    load([num2str(shot) '_' num2str(t) '_sys.mat']);
    eq = load(['eq' num2str(shot) '_' num2str(t) '.mat']);
    eq = eq.eq;
    
    eqs{i} = eq;
    
    [A(i,:,:), B(i,:,:)] = c2d(sys.As, sys.B, dt);    
    
    psibry(i) = eq.psibry;
    li(i) = eq.li;
    betap(i) = eq.betap;
    pres(i,:) = eq.pres;
    fpol(i,:) = eq.fpol;
    
    load(['coils' num2str(shot) '.mat'])
    k = find( floor(coils.t*1000) == t);
    ic(i,:) = coils.ic(k,:);
    iv(i,:) = coils.iv(k,:);
    ip(i,:) = eq.cpasma;
    
%     traj.W(i,:) = eq.W;
%     traj.WV(i) = sum(eq.W .* eq.V);
%     traj.PV(i) = sum(eq.pres .* eq.V);
    
  end
  
  method = 'pchip';
  traj.A      = interp1(t_snapshots, A,      tspan, 'previous');
  traj.B      = interp1(t_snapshots, B,      tspan, 'previous');
  traj.psibry = interp1(t_snapshots, psibry, tspan, method);
  traj.li     = interp1(t_snapshots, li,     tspan, method);
  traj.betap  = interp1(t_snapshots, betap,  tspan, method);
  traj.pres   = interp1(t_snapshots, pres,   tspan, method);
  traj.fpol   = interp1(t_snapshots, fpol,   tspan, method);
  traj.ic     = interp1(t_snapshots, ic,     tspan, method);
  traj.iv     = interp1(t_snapshots, iv,     tspan, method);
  traj.ip     = interp1(t_snapshots, ip,     tspan, method);
  traj.x      = [traj.ic traj.iv traj.ip'];

end








































