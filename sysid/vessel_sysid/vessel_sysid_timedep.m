% Vessel Resistance Fitting Algorithm
ccc
RAMPROOT = getenv('RAMPROOT');
load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);
option = 1;

shotlist = [204660];
starttimes = [0];
endtimes = [0.95];

% shotlist = [202376,202379,202381,202384];
% starttimes = [-2.5,-2.5,-1,-2];
% endtimes = [3.5,2,2,3];

shot = shotlist(1);
Ts = 0.01;
tstart = starttimes(1);
tend = endtimes(1);
tsample = tstart:Ts:tend;

% ================
% Initialize Model
% ================

% Option 1: identical to coneqt resvFit
if option == 1
  vacuum_system = load('NSTXU_vaccum_system.mat').NSTXU_vacuum_system;
  tok_data_struct = vacuum_system.build_inputs.tok_data_struct;
  Rxx = diag(vacuum_system.Rxx);
  Mxx = vacuum_system.Mxx;
  Mxx = blkdiag(Mxx,0);
  Rvv0 = Rxx(circ.iivx);
end

% Option 2: coneqt, lstar
if option == 2
  sys = load([RAMPROOT 'buildmodel/built_models/std_coneqt/204660_300_sys.mat']).sys;
  tok_data_struct = load('coneqt_nstxu_obj_config2016_6565.mat').tok_data_struct;
  Mxx = sys.lstar;
  Rxx = sys.rxx;
  Rvv0 = load('Rvv_fit.mat').Rvv_fit;
end

% Option 3: std, lstar
if option == 3
  sys = load([RAMPROOT 'buildmodel/built_models/std/204660_300_sys.mat']).sys;
  tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
  Mxx = sys.lstar;
  Rxx = sys.rxx;
  Rvv0 = Rxx(circ.iivx);
end

% Option 4: coneqt, mcc
if option == 4
  sys = load([RAMPROOT 'buildmodel/built_models/mcc_coneqt/204660_300_sys.mat']).sys;
  tok_data_struct = load('coneqt_nstxu_obj_config2016_6565.mat').tok_data_struct;
  Mxx = sys.lstar;
  Rxx = sys.rxx;
  Rvv0 = load('Rvv_fit.mat').Rvv_fit;
end

% Option 5: std, mcc
if option == 5
  sys = load([RAMPROOT 'buildmodel/built_models/mcc/204660_300_sys.mat']).sys;
  tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
  Mxx = sys.lstar;
  Rxx = sys.rxx;
  Rvv0 = Rxx(circ.iivx);
end

Mvc = Mxx(circ.iivx, circ.iicx);
Mvv = Mxx(circ.iivx, circ.iivx);
Mvp = zeros(circ.nvx,1);
Lvv0 = diag(Mvv);
Rvv0 = Rvv0 * 1000;

parameters = {'Mvv', Mvv; 'Mvc', Mvc; 'Mvp', Mvp; 'Rvv', Rvv0; 'Lvv', Lvv0};
odefun = 'vessel_dynamics_timedep';
sys = idgrey(odefun, parameters, 'd',{},Ts);

sys.Structure.Parameters(1).Free = false; % Mvv
sys.Structure.Parameters(2).Free = false; % Mvc
sys.Structure.Parameters(3).Free = false; % Mvp
sys.Structure.Parameters(4).Free = true;  % Rvv
sys.Structure.Parameters(5).Free = false; % Lvv

sys.Structure.Parameters(4).Minimum = Rvv0 ./ 100; 
sys.Structure.Parameters(4).Maximum = Rvv0 .* 100;

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

lpfreq = 1000; %Hz
ictspass = idealfilter(icts,[0,lpfreq],'pass');
ivtspass = idealfilter(ivts,[0,lpfreq],'pass');
iptspass = idealfilter(ipts,[0,lpfreq],'pass');

icdot = gradient(ictspass.Data', Ts)';
ipdot = gradient(iptspass.Data', Ts)';
icdot = smoothdata(icdot,1,'movmean',13);
ipdot = smoothdata(ipdot,1,'movmean',13);

y = double(ivtspass.Data);
u = double([icdot ipdot]);
shotdata = iddata(y, u, Ts);
x0 = y(1,:)';

% ============
% Fit to data
% ============
search_options.Algorithm = 'sqp';
search_options.FunctionTolerance = 5e-6;
search_options.StepTolerance = 5e-6;
search_options.MaxIterations = 30;
search_options.Advanced.TolFun = search_options.FunctionTolerance;
search_options.Advanced.TolX = search_options.StepTolerance;
search_options.Advanced.MaxIter = search_options.MaxIterations;
search_options.Advanced.Algorithm = search_options.Algorithm;
 


opt = greyestOptions('Display', 'on', 'InitialState', x0, ...
    'DisturbanceModel', 'none', 'Focus', 'simulation', ...
    'SearchMethod', 'auto','OutputWeight',eye(circ.nvx));

opt.SearchOptions.MaxIterations = 100;
opt.SearchOptions.Tolerance = 0.001;
% opt.SearchOptions.StepTolerance = 1e-6;
% opt.SearchOptions.FunctionTolerance = 1e-6;
  
% [A, B, C, D] = vessel_dynamics(Mvv, Mvc, Mvp, Rvv0, Lvv0, Ts);

%%
sys_est = greyest(shotdata, sys, opt);

%%
[yest,test,xest] = lsim(sys_est, u, tsample, x0);

figure
hold on
plot(test,yest,'b')
plot(tsample,y,'--r')

% save('sys_est_opt3','sys_est')

































