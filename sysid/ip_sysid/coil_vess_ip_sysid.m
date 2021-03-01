% Vessel Resistance Fitting Algorithm
ccc
RAMPROOT = getenv('RAMPROOT');
load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);
option = 1;
enforce_stability = 0;

shotlist = [204660];
starttimes = [0];
endtimes = [0.96];

Ts = 0.01;
shot = shotlist(1);
tstart = starttimes(1);
tend = endtimes(1);
tsample = tstart:Ts:tend;

parameter_times = linspace(tstart, tend, 30);

%%
% =========
% LOAD DATA
% =========

% Coil Currents
include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
        'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
        
icsignals = get_icsignals(shot, [], [], include_coils);
ivsignals = get_ivsignals(shot);
ipsignals = mds_fetch_signal(shot, 'efit01', '.RESULTS.AEQDSK:IPMEAS');

icts = timeseries(icsignals.sigs,icsignals.times);      
ivts = timeseries(ivsignals.sigs, ivsignals.times);
ipts = timeseries(ipsignals.sigs,ipsignals.times);      

icts = resample(icts,tsample);   
ivts = resample(ivts,tsample);
ipts = resample(ipts,tsample);   

% obtain derivatives
Tsmooth = 10;  % [ms]
nsmooth = floor(Tsmooth/1000/Ts);

ic = smoothdata(icts.Data, 1, 'movmean', nsmooth);
iv = smoothdata(ivts.Data, 1, 'movmean', nsmooth);
ip = smoothdata(ipts.Data, 1, 'movmean', nsmooth);

icdot = gradient(ic', Ts)';
ivdot = gradient(iv', Ts)';
ipdot = gradient(ip', Ts)';

icdot = smoothdata(icdot,1,'movmean',nsmooth);
ivdot = smoothdata(ivdot,1,'movmean',nsmooth);
ipdot = smoothdata(ipdot,1,'movmean',nsmooth);

% DO NOT use the filtered values for y here
y = double(ipts.Data);  

% DO NOT filter the voltages used in u here
u = double([icdot ivdot]);

shotdata = iddata(y, u, Ts);

x0 = y(1);

%%
% ================
% Initialize Model
% ================

fit_Lp = load('fit_Lp.mat').fit_Lp;
Lp_t = fit_Lp(parameter_times);

fit_Rp = load('fit_Rp5.mat').fit_Rp;
Rp_t = fit_Rp(parameter_times);

modeldir = [RAMPROOT '/buildmodel/built_models/mcc/'];
traj.t = tsample;
[traj.Mpc, traj.Mpv] = load_Mpc_Mpv(modeldir, tsample);


odefun = 'coil_vess_ip_dynamics';
order = [1 (circ.ncx+circ.nvx) 1];
file_args = {traj, parameter_times, Ts};
parameters = {Rp_t, Lp_t};

sys = idnlgrey(odefun,order,parameters,x0,Ts,'FileArgument',file_args);


sys.Parameters(1).Fixed = false;  % Rp_t
sys.Parameters(2).Fixed = true; % Lp_t

sys.Parameters(1).Minimum = 0;  % Rp_t
sys.Parameters(2).Minimum = 0;  % Lp_t



% ============
% Fit to data
% ============
opt = nlgreyestOptions('Display', 'on', 'SearchMethod', 'auto');

opt.SearchOptions.MaxIterations = 100;
opt.SearchOptions.StepTolerance = 1e-8;
opt.SearchOptions.FunctionTolerance = 1e-7;
  
%%
sys_est = nlgreyest(shotdata, sys, opt);


%%
Rp_t = sys.Parameters(1).Value;
figure
plot(parameter_times, Rp_t)

% Rp_ipfit.t = parameter_times;
% Rp_ipfit.value = Rp_t
% Rp_ipfit
% save('Rp_ipfit','Rp_ipfit')

simopt = simOptions;
simopt.InitialCondition = x0;

yest = sim(sys_est, u, simopt);

figure
hold on
plot(tsample, yest, 'b')
plot(tsample, y, '--r')
title([num2str(shot) ' Ip'], 'fontsize', 14)
ylabel('[A]', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
legend('True', 'Simulated')






























































