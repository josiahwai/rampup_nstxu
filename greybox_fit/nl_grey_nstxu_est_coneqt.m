ccc
RAMPROOT = getenv('RAMPROOT');

saveit = 1;
shot = 204660;
% grey_timebase = (30:10:360) / 1000;
% grey_timebase = (360:10:940) / 1000;
% grey_timebase = (30:10:940) / 1000;
grey_timebase = (270:10:940) / 1000;

ts = mean(diff(grey_timebase));

eqdir = [RAMPROOT '/eq/geqdsk_import/'];
modeldir = [RAMPROOT '/buildmodel/built_models/std_coneqt/'];
savedir = [RAMPROOT '/greybox_fit/fitted_models/'];

load('sim_inputs204660_smoothed_coneqt.mat')
load('nstxu_obj_config2016_6565.mat')
fit_Rp = load('fit_Rp5.mat').fit_Rp;

% coils = mds_fetch_current_voltage(shot,0);
coils = load('coils_greybox.mat').coils;

circ = nstxu2016_circ(tok_data_struct);

sim_timebase = sim_inputs.traj.tspan;
ps_voltages  = pinv(circ.Pcc_keep) * coils.v;
ps_voltages = interp1(coils.t, ps_voltages', grey_timebase);
ps_voltages = smoothdata(ps_voltages, 1, 'movmean', 13);

ytarg = sim_inputs.traj.x * circ.Pxx_keep';
ytarg = interp1(sim_timebase, ytarg, grey_timebase)';
ytarg = smoothdata(ytarg, 2, 'movmean', 13);
% close; plot(grey_timebase,ytarg); ylim([-1 1]*1.5e4)

% Initialize estimates for resistances
rc0 = diag(circ.Pcc' * diag(tok_data_struct.resc) * circ.Pcc);
% rv0 = diag(circ.Pvv' * diag(tok_data_struct.resv) * circ.Pvv);
rv0 = load('Rvv_fit.mat').Rvv_fit;
rp0_time = fit_Rp(grey_timebase);
[lp0, tdum] = read_Lp;
lp0_time = interp1(tdum, lp0, grey_timebase)';
voltage_scale = ones(circ.nu,1);
parameters = {rc0, rv0, rp0_time, lp0_time, voltage_scale};

fn = 'nl_grey_nstxu_model';
order = [circ.nxx_keep circ.nu circ.nx];
[~,k] = min(abs(sim_timebase - grey_timebase(1)));
x0 = sim_inputs.traj.x(k,:)';
Ts = mean(diff(sim_timebase));

% load file args
grey_timebase_shifted = grey_timebase - grey_timebase(1) + Ts;
lstar_invs = interp1(sim_timebase, sim_inputs.traj.lstari, grey_timebase);
file_args = {Ts, circ, lstar_invs, grey_timebase_shifted};

grey_init_sys = idnlgrey(fn,order,parameters,x0,Ts, 'FileArgument', file_args);
    
grey_init_sys.Parameters(1).Fixed(1:end) = false;  % rc 
grey_init_sys.Parameters(2).Fixed(1:end) = false;  % rv 
grey_init_sys.Parameters(3).Fixed(1:end) = false;  % rp_time 
grey_init_sys.Parameters(4).Fixed(1:end) = true;   % lp_time
grey_init_sys.Parameters(5).Fixed(1:end) = false;  % voltage_scale

grey_init_sys.Parameters(1).Minimum(1:end) = 0;  % rc 
grey_init_sys.Parameters(2).Minimum(1:end) = 0;  % rv 
grey_init_sys.Parameters(3).Minimum(1:end) = 0;  % rp_time 
grey_init_sys.Parameters(4).Minimum(1:end) = 0;  % lp_time
grey_init_sys.Parameters(5).Minimum(1:end) = 0;  % voltage_scale


opt = nlgreyestOptions;
wt.icx(1:circ.ncx_keep) = 10;
wt.ivx(1:circ.nvx) = 1;
wt.ip = 1;
opt.OutputWeight = diag([wt.icx wt.ivx wt.ip]);

grey_data = iddata(ytarg', ps_voltages, Ts);   

%%
grey_sys = nlgreyest(grey_data, grey_init_sys, opt);

%%

rcc = grey_sys.Parameters(1).Value;
rvv = grey_sys.Parameters(2).Value;
rp_t = grey_sys.Parameters(3).Value;
lp_t = grey_sys.Parameters(4).Value;
voltage_scale = grey_sys.Parameters(5).Value;

N = length(grey_timebase);
x = x0;
for i = 1:N
  yall(:,i) = circ.Pxx_keep * x;  
  xall(:,i) = x;
  t = grey_timebase_shifted(i);
  u = ps_voltages(i,:)';
  [x,y,A,B] = nl_grey_nstxu_modelAB(t, x, u, rcc, rvv, rp_t, lp_t, voltage_scale, file_args);        
  fittedAB.A(i,:,:) = A;
  fittedAB.B(i,:,:) = B;
end


figure
subplot(311)
hold on
plot(grey_timebase, yall(end, :), 'b')
plot(grey_timebase, ytarg(end,:), '--r')
set(gcf, 'Position', [680 248 467 730])
xlabel('Time [s]', 'fontsize', 14)
ylabel('Ip', 'fontsize', 14)
legend('Fit', 'Target', 'fontsize', 14)

subplot(312)
hold on
plot(grey_timebase, yall(1:8, :), 'b')
plot(grey_timebase, ytarg(1:8,:), '--r')
ylabel('Coil Currents', 'fontsize', 14)

subplot(313)
hold on
plot(grey_timebase, yall(9:48, :), 'b')
plot(grey_timebase, ytarg(9:48,:), '--r')
ylabel('Vessel Currents', 'fontsize', 14)


t = grey_timebase;
UserData = variables2struct(t, ytarg, ps_voltages, fittedAB);
grey_sys.UserData = UserData;

if saveit
  fn = ['grey_sys_' num2str(t(1)*1e3) '_' num2str(t(end)*1e3) '_coneqt.mat'];
  save([savedir fn], 'grey_sys');
end











































