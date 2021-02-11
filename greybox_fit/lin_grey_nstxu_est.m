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
time_ms = 100;
sys = load([modeldir num2str(shot) '_' num2str(time_ms) '_sys.mat']).sys;

sim_timebase = sim_inputs.tspan;
ps_voltages  = pinv(circ.Pcc_keep) * coils.v;

ps_voltages = interp1(coils.t, ps_voltages', sim_timebase);
ytarg = circ.Pxx_keep * sim_inputs.x_all;
Ts = mean(diff(sim_timebase));

% Initialize estimates for resistances
rc0 = diag(circ.Pcc' * diag(tok_data_struct.resc) * circ.Pcc);
rv0 = diag(circ.Pvv' * diag(tok_data_struct.resv) * circ.Pvv);
rp0 = median(fit_Rp(sim_timebase));
[lp0, tdum] = read_Lp;
lp0 = median(lp0);
parameters = {'rc', rc0; 'rv', rv0; 'rp0', rp0; 'lp0', lp0};

grey_init_sys = idgrey('lin_greybox_model_res', parameters, 'd', {circ, sys}, Ts);

grey_init_sys.Structure.Parameters(1).Free(1:end) = true;  % rc 
grey_init_sys.Structure.Parameters(2).Free(1:end) = false;  % rv 
grey_init_sys.Structure.Parameters(3).Free(1:end) = true;  % rp
grey_init_sys.Structure.Parameters(4).Free(1:end) = false;  % lp

grey_init_sys.Structure.Parameters(1).Minimum(1:end) = 0;  % rc 
grey_init_sys.Structure.Parameters(2).Minimum(1:end) = 0;  % rv 
grey_init_sys.Structure.Parameters(3).Minimum(1:end) = 1e-8;  % rp
grey_init_sys.Structure.Parameters(4).Minimum(1:end) = 0;  % lp


grey_data = iddata(ytarg', ps_voltages, Ts);   

opt = greyestOptions;
wt.icx(1:circ.ncx_keep) = 1e3;
wt.ivx(1:circ.nvx) = 1e-6;
wt.ip = 1e-6;
opt.OutputWeight = diag([wt.icx wt.ivx wt.ip]);

grey_sys = greyest(grey_data, grey_init_sys, opt);


rc = grey_sys.Structure.Parameters(1).Value;
rv = grey_sys.Structure.Parameters(2).Value;
rp = grey_sys.Structure.Parameters(3).Value;
lp = grey_sys.Structure.Parameters(4).Value;



%%

% create some data
N = length(sim_timebase);
x = sim_inputs.x0;
for i = 1:N
%   [A,B,C,D] = lin_greybox_model_res(rc0, rv0, rp0, lp0, Ts, circ, sys);
  [A,B,C,D] = lin_greybox_model_res(rc, rv, rp, lp, Ts, circ, sys);
  u = ps_voltages(i,:)';
  x = A*x + B*u;
  y(:,i) = C*x;  
  xall(:,i) = x;
end

figure
hold on
plot(sim_timebase, xall([circ.iicx], :), 'r')
plot(sim_timebase, ytarg(1:8,:), 'b')














































