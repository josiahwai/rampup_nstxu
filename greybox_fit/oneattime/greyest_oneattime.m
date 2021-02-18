ccc
RAMPROOT = getenv('RAMPROOT');

% ========
% SETTINGS
% ========
icoil_to_estimate = 1;
saveit = 0;
shot = 204660;
grey_timebase = (360:10:940) / 1000;

eqdir = [RAMPROOT '/eq/geqdsk_import/'];
savedir = [RAMPROOT '/greybox_fit/fitted_models/'];

modeldir = [RAMPROOT '/buildmodel/built_models/std/'];
traj = load('sim_inputs204660_smoothed.mat').sim_inputs.traj;
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;

% modeldir = [RAMPROOT '/buildmodel/built_models/std_coneqt/'];
% traj = load('sim_inputs204660_smoothed_coneqt.mat').sim_inputs.traj;
% tok_data_struct = load('coneqt_nstxu_obj_config2016_6565.mat').tok_data_struct;


% tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
% circ = nstxu2016_circ(tok_data_struct);
% modeldir = [RAMPROOT '/buildmodel/built_models/mcc/'];
% traj = load_interp_trajectory(shot, circ, grey_timebase, eqdir, modeldir);


% =============================
% Initialize Greybox estimator
% =============================
N = length(grey_timebase);
circ = nstxu2016_circ(tok_data_struct);

fit_Rp = load('fit_Rp5.mat').fit_Rp;
coils = load('coils_greybox.mat').coils;

sim_timebase = traj.tspan;
ps_voltages  = pinv(circ.Pcc_keep) * coils.v;
ps_voltages = interp1(coils.t, ps_voltages', grey_timebase);
ps_voltages = smoothdata(ps_voltages, 1, 'movmean', 13);

xtarg = traj.x(:,icoil_to_estimate);
xtarg = interp1(sim_timebase, xtarg, grey_timebase)';

% Initialize estimates for parameters
rss = [tok_data_struct.resc; tok_data_struct.resv; 1];
rxx = diag(circ.Pxx' * diag(rss) * circ.Pxx);

r0 = rxx(icoil_to_estimate);
l0 = tok_data_struct.mcc(icoil_to_estimate, icoil_to_estimate);
m_plasma_fac0 = 1;
voltage_fac0 = 1;
coil_voltage0 = ps_voltages(:,icoil_to_estimate);
parameters = {r0, l0, m_plasma_fac0, voltage_fac0, coil_voltage0};

fn = 'greymodel_oneattime';
order = [1 circ.nu circ.nx];
[~,k] = min(abs(sim_timebase - grey_timebase(1)));
x0 = traj.x(k,:)';
Ts = mean(diff(sim_timebase));

% load file args
grey_timebase_shifted = grey_timebase - grey_timebase(1) + Ts;
lstar_invs = interp1(sim_timebase, traj.lstari, grey_timebase);
file_args = {Ts, circ, lstar_invs, grey_timebase_shifted, traj, icoil_to_estimate, rxx};

grey_init_sys = idnlgrey(fn,order,parameters,x0,Ts, 'FileArgument', file_args);
    
grey_init_sys.Parameters(1).Fixed(1:end) = true;  % r0 
grey_init_sys.Parameters(2).Fixed(1:end) = true ; % l0
grey_init_sys.Parameters(3).Fixed(1:end) = true ; % m_plasma_fac
grey_init_sys.Parameters(4).Fixed(1:end) = true;  % voltage_fac
grey_init_sys.Parameters(5).Fixed(1:end) = false; % voltage_fac

grey_init_sys.Parameters(1).Minimum(1:end) = 0;  % r0
grey_init_sys.Parameters(2).Minimum(1:end) = 0;  % l0 
grey_init_sys.Parameters(3).Minimum(1:end) = 0;  % m_plasma_fac
grey_init_sys.Parameters(4).Minimum(1:end) = 0;  % voltage_fac

grey_data = iddata(xtarg, ps_voltages, Ts);   

%%
grey_sys = nlgreyest(grey_data, grey_init_sys);

%%
r = grey_sys.Parameters(1).Value;
l = grey_sys.Parameters(2).Value;
m_plasma_fac = grey_sys.Parameters(3).Value;
voltage_fac = grey_sys.Parameters(4).Value;
% coil_voltage = grey_sys.Parameters(5).Value;

x = x0;
for i = 1:N
  t = grey_timebase_shifted(i);
  u = ps_voltages(i,:)';
  [x,y] = greymodel_oneattime(t, x, u, r, l, m_plasma_fac, voltage_fac, coil_voltage, file_args);      
  yall(i) = y;      
end

figure
hold on
title(['Coil ' num2str(icoil_to_estimate)])
plot(grey_timebase, yall, '-b')
plot(grey_timebase, xtarg, '--r')
legend('Fit', 'Target', 'fontsize', 14)










































