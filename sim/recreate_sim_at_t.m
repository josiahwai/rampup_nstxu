ccc
time_ms = 180;

t = time_ms/1000;
load('sim.mat')
load(['eq204660_' num2str(time_ms) '.mat'])
load('nstxu_obj_config2016_6565.mat')
load('sim_inputs204660.mat')

% True currents
traj = sim_inputs.traj;
[~,i] = min(abs(t - sim_inputs.tspan));
xt = traj.x(i,:)';

% Predicted simulated currents
load('xsim_pred.mat')
xs = xsim_pred(:,i);

% Simulated currents
% [~,i] = min(abs(t - sim.tspan));
% xs = sim.x(:,i);

% xs(9:48) = 0;
xs(9:48) = xt(9:48);
% xs(end) = xt(end);


[spec, config] = make_gsdesign_inputs(xs, tok_data_struct, eq, traj, t, sim_inputs.tspan);

eq = gsdesign(spec, eq, config);








