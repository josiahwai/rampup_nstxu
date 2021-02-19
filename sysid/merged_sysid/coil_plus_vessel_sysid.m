% Vessel Resistance Fitting Algorithm
ccc
RAMPROOT = getenv('RAMPROOT');
load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);
option = 1;
enforce_stability = 0;

shotlist = [204660];
starttimes = [0];
endtimes = [1];

% shotlist = [202376,202379,202381,202384];
% starttimes = [-2.5,-2.5,-1,-2];
% endtimes = [3.5,2,2,3];

Ts = 2e-4;
shot = shotlist(1);
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
  Rxx(end+1) = 0;
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

Mvp = Mxx(circ.iivx, circ.iipx);
Mcp = Mxx(circ.iicx, circ.iipx);
Mvc = Mxx(circ.iivx, circ.iicx);
Mvv = Mxx(circ.iivx, circ.iivx);
Lvv0 = diag(Mvv);
Rvv0 = Rvv0 * 1000;

load('ext_fit.mat')
fit_coils = [1     2     5     6     8     9    10    13];
for j = 1:length(fit_coils)
  k = find(fit_coils(j) == ext_fit.icoil);
  Rext_mOhm(j) = ext_fit.Rext_fit(k) * 1000;
  Lext_mH(j)   = ext_fit.Lext_fit(k) * 1000;    
end

% remove the plasma portion of Mxx, Rxx
Mxx(end,:) = [];   
Mxx(:,end) = [];
Rxx(end) = [];

% turn off the unused coils
Rxx_use = Rxx;
Rxx_use(circ.iicx) = 10;
Rxx_use(fit_coils) = Rxx(fit_coils);
Rxx = Rxx_use;

file_args = {Mxx, Rxx, fit_coils, circ, Rext_mOhm, Lext_mH, enforce_stability};
parameters = {'Mvv', Mvv; 'Mvc', Mvc; 'Mvp', Mvp; 'Rvv_mOhm', Rvv0; 'Lvv', Lvv0; 'Mcp' Mcp};
odefun = 'coil_plus_vessel_dynamics';
sys = idgrey(odefun, parameters, 'd', file_args, Ts, 'InputDelay', 3);

sys.Structure.Parameters(1).Free = false; % Mvv
sys.Structure.Parameters(2).Free = false; % Mvc
sys.Structure.Parameters(3).Free = false; % Mvp
sys.Structure.Parameters(4).Free = true;  % Rvv
sys.Structure.Parameters(5).Free = false; % Lvv
sys.Structure.Parameters(6).Free = false; % Mcp


sys.Structure.Parameters(2).Minimum = 0;    % Mvc
sys.Structure.Parameters(4).Minimum = Rvv0 ./ 100; % Rvv
sys.Structure.Parameters(4).Maximum = Rvv0 .* 100; % Rvv
sys.Structure.Parameters(6).Minimum = 0;    % Mcp
% =========
% LOAD DATA
% =========

% Coil Currents
include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
        'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
        
icsignals = get_icsignals(shot, [], [], include_coils);
ivsignals = get_ivsignals(shot);
ipsignals = mds_fetch_signal(shot, 'efit01', '.RESULTS.AEQDSK:IPMEAS');
vsignals = get_vobjcsignals(shot, [], [], include_coils);

icts = timeseries(icsignals.sigs,icsignals.times);      
ivts = timeseries(ivsignals.sigs, ivsignals.times);
ipts = timeseries(ipsignals.sigs,ipsignals.times);      
vts = timeseries(vsignals.sigs, vsignals.times);

icts = resample(icts,tsample);   
ivts = resample(ivts,tsample);
ipts = resample(ipts,tsample);   
vts = resample(vts, tsample);

% BEWARE, filtering affects the absolute magnitudes. Therefore one can 
% safely use the filtered values ONLY for derivatives.
% lpfreq = 1000; %Hz
% ictsfilt = idealfilter(icts,[0,lpfreq],'pass');
% ivtsfilt = idealfilter(ivts,[0,lpfreq],'pass');
% iptsfilt = idealfilter(ipts,[0,lpfreq],'pass');

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
y = double([icts.Data ivts.Data]);  

% DO NOT filter the voltages used in u here
u = double([vts.Data icdot ipdot]);

shotdata = iddata(y, u, Ts);

load('sim_inputs204660_smoothed.mat')
x0 = sim_inputs.traj.x(1,:)';
x0(end) = [];

% x0 = y(1,:)';

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
 
wt.ic = ones(circ.ncx,1) * 1e-6;
wt.iv = ones(circ.nvx,1) * 1; 
wt = diag([wt.ic; wt.iv]);

opt = greyestOptions('Display', 'on', 'InitialState', x0, ...
    'DisturbanceModel', 'none', 'Focus', 'simulation', ...
    'SearchMethod', 'auto','OutputWeight', wt, ...
    'EnforceStability', true);

opt.SearchOptions.MaxIterations = 100;
opt.SearchOptions.Tolerance = 0.001;
% opt.SearchOptions.StepTolerance = 1e-6;
% opt.SearchOptions.FunctionTolerance = 1e-6;
  
%%
sys_est = greyest(shotdata, sys, opt);

%%

% Debugging:
% sys = idgrey(odefun, parameters, 'd', file_args, Ts, 'InputDelay', 3);
% u = double([vts.Data icdot ipdot]);
% u = double(vts.Data);
% [yest,t,xest] = lsim(sys, u, tsample, x0);


[yest,t,xest] = lsim(sys_est, u, tsample, x0);

figure
hold on
plot(t, xest(:,circ.iicx), 'b')
plot(t, y(:,circ.iicx), '--r')
title([num2str(shot) ' Coil Currents'], 'fontsize', 14)
ylabel('[A]', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
mylegend({'True', 'Simulated'}, {'--','-'}, [], {'r','b'}, [], 'Northeast', 14);

load('coils_greybox.mat')
figure
hold on
plot(t,xest(:,circ.iivx),'b')
plot(coils.t, coils.iv,'--r')
xlim([0 0.9])
title([num2str(shot) ' Vessel Currents'], 'fontsize', 14)
ylabel('[A]', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
mylegend({'True', 'Simulated'}, {'--','-'}, [], {'r','b'}, [], 'Northeast', 14);

% save('merge_sys_est', 'sys_est')


Rvv = sys_est.Structure.Parameters(4).Value;
figure
vvlabels = categorical(circ.vvnames);
bar(vvlabels, [Rvv0 Rvv])
ylabel('Resistance [mOhm]')
legend('Original', 'Fit', 'fontsize', 14)
































































