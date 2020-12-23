clear; clc; close all
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

% ==============
% Define problem
% ==============
shot = 204660;
[A,B,x_t0] = read_x(shot, 60);
[~,~,xd_tf] = read_x(shot, 120);

t0 = 0.060;
tf = 0.100;
dt = 1e-5;
tspan = t0:dt:tf;
nt = length(tspan);

% slow down the ultra-fast poles for numerical stability
[V,D] = eig(A);
D = diag(D);
D(abs(D) > 1e4) = -1e4;
A = real(V*diag(D)/V);

wt.ic = ones(1,13);
wt.iv = ones(1,40) * 0;
wt.Ip = 1e-4;
wt.u  = ones(1,13);

wt.Ip(1) = 1e-4;

Q = 1 / (tf-t0) * diag([wt.ic wt.iv wt.Ip]);
R = 1 / (tf-t0) * diag(wt.u);
Qf = diag([wt.ic wt.iv wt.Ip]);


% ===============
% Adjoint looping
% ===============
uk = zeros(nt,13);

for i = 1:10
  eps = 1e-5;
  [J1, grad] = adjoint_loop( uk             , x_t0, xd_tf, A, B, R, Q, Qf, tspan);
  [J2, ~]    = adjoint_loop( uk - eps*grad  , x_t0, xd_tf, A, B, R, Q, Qf, tspan);
  [J3, ~]    = adjoint_loop( uk - 2*eps*grad, x_t0, xd_tf, A, B, R, Q, Qf, tspan);
  
  M = [1 0 0; 1 eps eps^2; 1 2*eps 4*eps^2];
  J = [J1 J2 J3]';
  a = M \ J;
  eps = -a(2) / (2*a(3));

  uk = uk - eps*grad;
end

%%
close all
% Plot results
[x_tf, x_all] = state_dynamics(A,B,x_t0,uk,tspan);

icoil = [54];
figure
hold on
plot(tspan, x_all(:,icoil), 'b')
plot([t0 tf], [x_t0(icoil) xd_tf(icoil)]', '--r')
xlabel('Time')
legend('Ip', 'Ip target', 'fontsize', 18)
ylabel('A')

figure
subplot(211)
plot(tspan,x_all)
title('Coil currents')
ylim([-1 1] * 3e4)
subplot(212)
i = [1:13];
bar([x_t0(i) x_tf(i) xd_tf(i)])
title('Coil currents at final time')
legend('x0', 'xf', 'xf target', 'fontsize', 14)





% ============
% Adjoint Loop
% ============
function [grad_cost, grad] = adjoint_loop( u, x_t0, xd_tf, A, B, R, Q, Qf, tspan)
  
  [x_tf, x_all] = state_dynamics(A,B,x_t0,u,tspan);
  xadj_tf = Qf * (x_tf - xd_tf);
  xadj_all = adjoint_dynamics(xadj_tf, x_all, x_t0, xd_tf, A, Q, tspan);
  grad = u*R + xadj_all * B;
  
  grad_cost_integrand = zeros(size(tspan));
  for i = 1:length(tspan)
    grad_cost_integrand(i) = 0.5 * norm(grad(i,:))^2;
  end
  grad_cost = trapz(tspan, grad_cost_integrand) + 0.5 * (x_tf-xd_tf)'*Qf*(x_tf-xd_tf);  
end




% ==============
% State dynamics
% ==============
function [x_tf, x_all] = state_dynamics(A,B,x0,u_all,tspan)

function xdot = statefun(t,x,A,B,u_all,tspan)  
  i = min(find(t < [tspan inf], 1), length(tspan)) - 1;
  u = linear_interp(t, tspan(i), tspan(i+1), u_all(i,:)', u_all(i+1,:)');  
  xdot = A*x + B*u;  
end

[~,x_all] = ode45(@(t,x) statefun(t,x,A,B,u_all,tspan), tspan, x0);

x_tf = x_all(end,:)';
end


  
% ================
% Adjoint dynamics
% ================
function xadj_all = adjoint_dynamics(xadj_tf, x_all, x_t0, xd_tf, A, Q, tspan)

function xadj_dot = adjointfun(t, xadj, x_all, x_t0, xd_tf, A, Q, tspan)
  i = min(find(t < [tspan inf], 1), length(tspan)) - 1;  
  x = linear_interp(t, tspan(i), tspan(i+1), x_all(i,:)', x_all(i+1,:)');  
  xd = linear_interp(t, tspan(1), tspan(end), x_t0, xd_tf);  
  xadj_dot = -A'*xadj - Q * (x - xd);  
end

[~,xadj_all] = ode45(@(t,xadj) adjointfun(t, xadj, x_all, x_t0, xd_tf, A, Q, tspan),...
  flip(tspan), xadj_tf);

end



% ==============
% Read A,B,x
% ==============
function [A,B,x] = read_x(shot, time_ms)

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


function z = linear_interp(t,t0,t1,z0,z1)
  c = (t-t0)/(t1-t0);
  z = (1-c)*z0 + c*z1;
end
































