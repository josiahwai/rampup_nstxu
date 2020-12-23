clear all; clc; close all
clear gsdesign spec init config
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

% ==============
% Define problem
% ==============
shot = 204660;
[Ac,Bc,x_t0,eq60] = read_x(shot, 60);
[~,~,xd_tf,eq120] = read_x(shot, 120);

t0 = 0.060;
tf = 0.120;
N = 2;
tspan = linspace(t0,tf,N+1);
dt = mean(diff(tspan));
nx = length(x_t0);
nu = 13;
iic = 1:13;
iiv = 14:53;
iip = 54;

[A,B] = c2d(Ac,Bc,dt);

wt.ic = ones(1,13);
wt.iv = ones(1,40) * 1e-3;
wt.Ip = 1e-3;
wt.u  = ones(1,13) * 1e-3;

wt.ic(1) = 10;

Q = diag([wt.ic wt.iv wt.Ip]);
Qf = diag([wt.ic wt.iv wt.Ip]);
Qhat = blkdiag(kron(eye(N-1),Q), Qf);

R = diag(wt.u);
Rhat = kron(eye(N),R);

% prediction model
Chat = kron(eye(N), eye(nx));
Apow  = eye(nx);
F  = zeros(N*size(B));
OB  = [];
for i = 1:N
  F  = F  + kron(diag(ones(N-i+1,1), -i+1), Apow*B);
  Apow  = A*Apow;
  OB  = [OB; Apow];
end

H = F'*Chat'*Qhat*Chat*F + Rhat;

% manipulate for numerical stability
if max(max(H-H')) < 1e-10, H = (H+H')/2; end
Linv = inv(chol(H,'lower'));
Linv(abs(Linv) < 1e-10) = 0; 

% future targets
rhat = [];
for k = 1:N
  r = x_t0 + k/N * (xd_tf - x_t0);  
  rhat = [rhat; r];
end

f = F' * Chat' * Qhat * (Chat * OB * x_t0 - rhat);

% mpcqpsolver options
Aeq = [];
beq = zeros(0,1);
Aineq = [];
bineq = zeros(0,1);
if ~exist('iA','var'), iA = false(size(bineq)); end
opt = mpcqpsolverOptions;

[Uhat,status,iA] = mpcqpsolver(Linv, f, Aineq, bineq, Aeq, beq, iA, opt);
Uhat = reshape(Uhat, nu, []);     

[x_all, x_tf] = state_dynamics(A, B, x_t0, N, Uhat);




% ==========
% PLOT STUFF
% ==========
% eq_tf = coil_currents_to_eq(eq60, x_tf(iic), x_tf(iiv), x_tf(iip));
% close all
% figure
% plot_eq(eq120)
% figure
% plot_eq(eq_tf)


% figure
% subplot(211)
% plot(tspan, x_all')
% subplot(212)
% plot(tspan(1:end-1), Uhat')

icoil = 1:5;
figure
subplot(211)
traj.x = [x_t0 xd_tf];
hold on
plot(tspan, x_all(icoil,:), 'b')
plot([t0 tf], traj.x(icoil,:), '--r')
subplot(212)
plot(tspan(1:end-1), Uhat)


% figure
% subplot(311)
% bar([x_t0(iic) x_tf(iic) xd_tf(iic)])
% subplot(312)
% bar([x_t0(iiv) x_tf(iiv) xd_tf(iiv)])
% subplot(313)
% bar([x_t0(iip) x_tf(iip) xd_tf(iip)])


% ==============
% State dynamics
% ==============
function [x_all, x_tf] = state_dynamics(A, B, x_t0, N, u_all)
  x = x_t0;
  x_all = zeros(length(x_t0), N+1);
  x_all(:,1) = x_t0;
  for i = 1:N
    x = A*x + B*u_all(:,i);
    x_all(:,i+1) = x;      
  end
  x_tf = x_all(:,end);
end


% ==============
% Read A,B,x
% ==============
function [A,B,x,eq] = read_x(shot, time_ms)

  load([num2str(shot) '_' num2str(time_ms) '_sys.mat']);
  load(['eq' num2str(shot) '_' num2str(time_ms) '.mat']);
  Ip = eq.cpasma;
  load(['coils' num2str(shot) '.mat'])
  i = find( floor(coils.t*1000) == time_ms);
  if isempty(i)
    warning('Could not import coil currents')
  end
  ic = coils.ic(i,:)';
  iv = coils.iv(i,:)';

  x = [ic; iv; Ip];
  A = sys.As;
  B = sys.B;
end































