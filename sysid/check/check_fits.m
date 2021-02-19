% Vessel Resistance Fitting Algorithm
ccc
RAMPROOT = getenv('RAMPROOT');

% shot = 204330;    
shot = 204660;
starttimes = [0];
endtimes = [0.95];

Ts = 0.01;
tstart = starttimes(1);
tend = endtimes(1);
tsample = tstart:Ts:tend;

vacuum_system = load('NSTXU_vaccum_system.mat').NSTXU_vacuum_system;
tok_data_struct = vacuum_system.build_inputs.tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);

Rxx = diag(vacuum_system.Rxx);
Mxx = vacuum_system.Mxx;
Rvv0 = Rxx(circ.iivx);

sys_est = load('sys_est_opt1.mat').sys_est;
Mvc = sys_est.Structure.Parameters(2).Value;
Rvv = sys_est.Structure.Parameters(4).Value;

ext_fit = load('ext_fit.mat').ext_fit;
fit_coils = [1     2     5     6     8     9    10    13];
for j = 1:length(fit_coils)
  k = find(fit_coils(j) == ext_fit.icoil);
  Rext_mOhm(j) = ext_fit.Rext_fit(k) * 1000;
  Lext_mH(j)   = ext_fit.Lext_fit(k) * 1000;
end

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

lpfreq = 1000; %Hz
ictspass = idealfilter(icts,[0,lpfreq],'pass');
ivtspass = idealfilter(ivts,[0,lpfreq],'pass');
iptspass = idealfilter(ipts,[0,lpfreq],'pass');
vtspass = idealfilter(vts, [0 lpfreq], 'pass');

icdot = gradient(ictspass.Data', Ts)';
ivdot = gradient(ivtspass.Data', Ts)';
ipdot = gradient(iptspass.Data', Ts)';

icdot = smoothdata(icdot,1,'movmean',13);
ivdot = smoothdata(ivdot,1,'movmean',13);
ipdot = smoothdata(ipdot,1,'movmean',13);

ic = double(ictspass.Data);
iv = double(ivtspass.Data);
ip = double(iptspass.Data);

ic0 = ic(1,:)';
iv0 = iv(1,:)';
ip0 = ip(1);
x0 = [ic0; iv0];

%%
% ================
% Integrate system
% ================
lstar = load([RAMPROOT 'buildmodel/built_models/std/204660_300_sys.mat']).sys.lstar;
Mcp = lstar(circ.iicx,end);
Mvp = lstar(circ.iivx,end);

u1 = vtspass.Data;
u2 = double([icdot ivdot ipdot]);
u = [u1 u2];

odefun = 'coil_and_vessel_dynamics';
parameters = {'Mxx' Mxx; 'Rxx' Rxx; 'fit_coils' fit_coils; ...
  'Rext_mOhm' Rext_mOhm; 'Lext_mH' Lext_mH; 'Mvc' Mvc; 'Rvv' Rvv};

sys = idgrey(odefun, parameters, 'd', {circ, Mcp, Mvp}, Ts);

[yest,t,xest] = lsim(sys, u, tsample, x0);

figure
hold on
plot(t, xest(:,circ.iicx), 'b')
plot(t, ic, '--r')

figure
hold on
plot(t, xest(:,circ.iivx), 'b')
plot(t, iv, '--r')





























