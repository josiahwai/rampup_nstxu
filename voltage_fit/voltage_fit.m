ccc
RAMPROOT = getenv('RAMPROOT');

% ========
% SETTINGS
% ========
saveit = 0;
shot = 204660;
fit_timebase = (30:10:940) / 1000;
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
N = length(fit_timebase);

% coils = mds_fetch_current_voltage(shot,0);
coils = load('coils_greybox.mat').coils;

ps_voltages  = pinv(circ.Pcc_keep) * coils.v;
ps_voltages = interp1(coils.t, ps_voltages', fit_timebase);
ps_voltages = smoothdata(ps_voltages, 1, 'movmean', 13);

xtarg = traj.x;
xtarg = interp1(sim_timebase, xtarg, fit_timebase)';
xtarg = smoothdata(xtarg, 2, 'movmean', 13);

[~,k] = min(abs(sim_timebase - fit_timebase(1)));
x0 = traj.x(k,:)';
% ==============
% Find voltages
% ==============
rcc = diag(circ.Pcc' * diag(tok_data_struct.resc) * circ.Pcc);
rvv = diag(circ.Pvv' * diag(tok_data_struct.resv) * circ.Pvv);
rp_t = load('rp_t').rp_t;
voltage_scale = ones(circ.nu,1);
Ts = mean(diff(fit_timebase));
fit_timebase_shifted = fit_timebase - fit_timebase(1) + Ts;
lstar_invs = interp1(sim_timebase, traj.lstari, fit_timebase);
file_args = {Ts, circ, lstar_invs, fit_timebase_shifted};


for i = 1:N-1 
  t = fit_timebase(i);
  udum = ps_voltages(i,:)';
  x = xtarg(:,i);
  [~,~,A,B] = nl_grey_nstxu_modelAB(t, x, udum, rcc, rvv, rp_t, [], voltage_scale, file_args);  
  
  x_next = xtarg(:,i+1);
  u = pinv(B) * (x_next - A*x);
  ufit(i,:) = u;
end


% ufit = smoothdata(ufit,1);

figure
plot(fit_timebase(1:end-1), ufit)




x = x0;
xfit(i,:) = x;
for i = 1:N-1 
  t = fit_timebase(i);
  u = ufit(i,:)';
  [x,y,A,B] = nl_grey_nstxu_modelAB(t, x, u, rcc, rvv, rp_t, [], voltage_scale, file_args);    
  xfit(i+1,:) = x;    
end
xfit = xfit';

%%
figure
subplot(311)
hold on
plot(fit_timebase, xfit(circ.iipx, :), 'b')
plot(fit_timebase, xtarg(circ.iipx,:), '--r')
set(gcf, 'Position', [680 248 467 730])
xlabel('Time [s]', 'fontsize', 14)
ylabel('Ip', 'fontsize', 14)
legend('Fit', 'Target', 'fontsize', 14)

subplot(312)
hold on
plot(fit_timebase, xfit(circ.iicx, :), 'b')
plot(fit_timebase, xtarg(circ.iicx,:), '--r')
ylabel('Coil Currents', 'fontsize', 14)
legend('Fit', 'Target', 'fontsize', 14)

subplot(313)
hold on
plot(fit_timebase, xfit(circ.iivx, :), 'b')
plot(fit_timebase, xtarg(circ.iivx,:), '--r')
ylabel('Vessel Currents', 'fontsize', 14)
legend('Fit', 'Target', 'fontsize', 14)



































