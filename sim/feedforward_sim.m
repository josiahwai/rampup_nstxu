% load x0, Uhat
% define: spec, init, config from: x, traj.li/betap/rzbbbs
% eq = gsdesign(spec,init,config)
% model = build_nstxu_system_fun()
% [A,B] = c2d(A,B)
% x = Ax + Bu
% save: eq.psizr, x

clear all; clc; close all
clear gsdesign spec init config

% ========
% Settings
% ========
load('nstxu_obj_config2016_6565.mat')
load('sim_inputs204660.mat');
struct_to_ws(sim_inputs);

% interpolate to different timebase
tspan0 = tspan;
dt0 = mean(diff(tspan0));
Uhat0 = Uhat;

N = 40;
tspan = linspace(tspan0(1), tspan0(end), N+1);
dt = mean(diff(tspan));
Uhat = interp1(tspan0(1:end-1), Uhat0', tspan(1:end-1), 'linear', 'extrap')';


% initialize
dt = mean(diff(tspan));
N = length(tspan) - 1;
sim.eq0 = eq0;
sim.x = x0;
sim.psizr(1,:,:) = eq0.psizr;
eq = eq0;
x  = x0;

% simulate
for i = 1:N
  t = tspan(i);  
  t_ms = t*1000
  opts.islimited = t_ms < 220;
  opts.Te_res = min( t_ms/0.3, 4000);  
  sys = build_nstxu_system_fun(eq, tok_data_struct, opts);      
  [A,B] = c2d(sys.As, sys.B, dt);
  u = Uhat(:,i);
  x = A*x + B*u;  
  [spec, config] = make_gsdesign_inputs(x, tok_data_struct, eq, traj, t, tspan0);
  eq = gsdesign(spec, eq, config);
    
  % write to sim
  sim.x(:,i+1) = x;
  sim.psizr(i+1,:,:)  = eq.psizr;
  sim.psibry(i+1) = eq.psibry;
end


%%
% =============
% Plot results
% =============
x_all = sim.x;

figure
sgtitle('Shot 204660')
subplot(221)
title('Coil currents', 'fontsize', 14, 'fontweight', 'bold')
icoil = 1:8;
hold on
plot(tspan, x_all(icoil,:), 'b')
plot(tspan0, traj.x(:,icoil)', '--r')
ylabel('A')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

subplot(222)
title('Vessel currents', 'fontsize', 14, 'fontweight', 'bold')
icoil = 9:48;
hold on
plot(tspan, x_all(icoil,:), 'b')
plot(tspan0, traj.x(:,icoil)', '--r')
ylabel('A')
mylegend({'Simulation', 'Experiment'}, {'-', '--'}, 1.5, {'b', 'r'});

subplot(223)
title('Ip', 'fontsize', 14, 'fontweight', 'bold')
hold on
plot(tspan, x_all(end,:), 'b')
plot(tspan0, traj.x(:,end)', '--r')
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














