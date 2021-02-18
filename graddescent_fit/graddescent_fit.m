ccc
RAMPROOT = getenv('RAMPROOT');

% ========
% SETTINGS
% ========
saveit = 1;
shot = 204660;
fit_timebase = (270:10:940) / 1000;
load('nstxu_obj_config2016_6565.mat')
eqdir = [RAMPROOT '/eq/geqdsk_import/'];
modeldir = [RAMPROOT '/buildmodel/built_models/std/'];
savedir = [RAMPROOT '/greybox_fit/fitted_models/'];
traj = load('sim_inputs204660_smoothed.mat').sim_inputs.traj;

% ==============
% Define targets
% ==============
circ = nstxu2016_circ(tok_data_struct);
sim_timebase = traj.tspan;

% coils = mds_fetch_current_voltage(shot,0);
coils = load('coils_greybox.mat').coils;

ps_voltages  = pinv(circ.Pcc_keep) * coils.v;
ps_voltages = interp1(coils.t, ps_voltages', fit_timebase);
ps_voltages = smoothdata(ps_voltages, 1, 'movmean', 13);

xtarg = traj.x;
xtarg = interp1(sim_timebase, xtarg, fit_timebase)';
xtarg = smoothdata(xtarg, 2, 'movmean', 13);

% =====================
% Initial guess for fit
% =====================
rc0 = diag(circ.Pcc' * diag(tok_data_struct.resc) * circ.Pcc);
rv0 = diag(circ.Pvv' * diag(tok_data_struct.resv) * circ.Pvv);
fit_Rp = load('fit_Rp5.mat').fit_Rp;
rp0_time = fit_Rp(fit_timebase);


% =====================
% Integrate and compare
% =====================
[~,k] = min(abs(sim_timebase - fit_timebase(1)));
x0 = traj.x(k,:)';

% arguments needed for integrator
Ts = mean(diff(fit_timebase));
fit_timebase_shifted = fit_timebase - fit_timebase(1) + Ts;
lstar_invs = interp1(sim_timebase, traj.lstari, fit_timebase);
file_args = {Ts, circ, lstar_invs, fit_timebase_shifted};
voltage_scale = ones(circ.nu,1);

rcc = rc0;
rvv = rv0;
rp_t = rp0_time;

N = length(fit_timebase);
x = x0;
for i = 1:N
  xall(:,i) = x;
  t = fit_timebase_shifted(i);
  u = ps_voltages(i,:)';
  [x,y,A,B] = nl_grey_nstxu_modelAB(t, x, u, rcc, rvv, rp_t, [], voltage_scale, file_args);        
end




% ================
% Plot differences
% ================
figure
subplot(311)
hold on
plot(fit_timebase, xall(circ.iipx, :), 'b')
plot(fit_timebase, xtarg(circ.iipx,:), '--r')
set(gcf, 'Position', [680 248 467 730])
xlabel('Time [s]', 'fontsize', 14)
ylabel('Ip', 'fontsize', 14)
legend('Fit', 'Target', 'fontsize', 14)

subplot(312)
hold on
plot(fit_timebase, xall(circ.ikeep, :), 'b')
plot(fit_timebase, xtarg(circ.ikeep,:), '--r')
ylabel('Coil Currents', 'fontsize', 14)
legend('Fit', 'Target', 'fontsize', 14)

subplot(313)
hold on
plot(fit_timebase, xall(circ.iivx, :), 'b')
plot(fit_timebase, xtarg(circ.iivx,:), '--r')
ylabel('Vessel Currents', 'fontsize', 14)
legend('Fit', 'Target', 'fontsize', 14)












































