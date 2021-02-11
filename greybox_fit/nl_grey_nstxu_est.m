ccc
shot = 204660;
times = [30:10:960];
ts = mean(diff(times))/1000;

eqdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk_import/';
modeldir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/mcc/';
load('sim_inputs204660_mcc.mat')
load('nstxu_obj_config2016_6565.mat')
fit_Rp = load('fit_Rp5.mat').fit_Rp;

% coils = mds_fetch_current_voltage(shot,0);
coils = load('coils_greybox.mat').coils;

circ = nstxu2016_circ(tok_data_struct);

sim_timebase = sim_inputs.tspan;
ps_voltages  = pinv(circ.Pcc_keep) * coils.v;

ps_voltages = interp1(coils.t, ps_voltages', sim_timebase);
ytarg = circ.Pxx_keep * sim_inputs.x_all;

sys = load([modeldir num2str(shot) '_100_sys.mat']).sys;

% Initialize estimates for resistances
rc0 = diag(circ.Pcc' * diag(tok_data_struct.resc) * circ.Pcc);
rv0 = diag(circ.Pvv' * diag(tok_data_struct.resv) * circ.Pvv);
rp0_time = fit_Rp(sim_timebase);
[lp0, tdum] = read_Lp;
lp0_time = interp1(tdum, lp0, sim_timebase)';
voltage_scale = ones(circ.nu,1);
parameters = {rc0, rv0, rp0_time, lp0_time, voltage_scale};

fn = 'nl_grey_nstxu_model';
order = [circ.nxx_keep circ.nu circ.nx];
x0 = sim_inputs.x0;
Ts = mean(diff(sim_timebase));

sim_timebase_shifted = sim_timebase - sim_timebase(1) + Ts;
file_args = {Ts, circ, sys.lstar, sim_timebase_shifted};

grey_init_sys = idnlgrey(fn,order,parameters,x0,Ts, 'FileArgument', file_args);
    
grey_init_sys.Parameters(1).Fixed(1:end) = false;  % rc 
grey_init_sys.Parameters(2).Fixed(1:end) = false;  % rv 
grey_init_sys.Parameters(3).Fixed(1:end) = false;  % rp_time 
grey_init_sys.Parameters(4).Fixed(1:end) = true;   % lp_time
grey_init_sys.Parameters(5).Fixed(1:end) = false;  % voltage_scale


opt = nlgreyestOptions;

grey_data = iddata(ytarg', ps_voltages, Ts);   

grey_sys = nlgreyest(grey_data, grey_init_sys);


%%

rcc = grey_sys.Parameters(1).Value;
rvv = grey_sys.Parameters(2).Value;
rp_t = grey_sys.Parameters(3).Value;
lp_t = grey_sys.Parameters(4).Value;
voltage_scale = grey_sys.Parameters(5).Value;

N = length(sim_timebase);
x = sim_inputs.x0;
for i = 1:N
  t = sim_timebase_shifted(i);
  u = ps_voltages(i,:)';
  [x,y] = nl_grey_nstxu_model(t, x, u, rcc, rvv, rp_t, lp_t, voltage_scale, file_args);      
  yall(:,i) = y;  
  xall(:,i) = x;
end

figure
hold on
plot(sim_timebase, yall(1:8, :), 'b')
plot(sim_timebase, ytarg(1:8,:), '--r')

figure
hold on
plot(sim_timebase, yall(9:48, :), 'b')
plot(sim_timebase, ytarg(9:48,:), '--r')


figure
hold on
plot(sim_timebase, yall(end, :), 'b')
plot(sim_timebase, ytarg(end,:), '--r')















































