ccc
RAMPROOT = getenv('RAMPROOT');

% Timing
Ts = .01;
shot = 204660;
tstart = 0;
tend = 1;
tsample = tstart:Ts:tend;
N = length(tsample);
t = tsample;


% Load fitted vacuum model parameters
vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);
Mxx = vac_sys.sysid_fits.Mxx;
Rxx = vac_sys.sysid_fits.Rxx;


% Load the time-dep plasma parameters
modeldir = [RAMPROOT '/buildmodel/built_models/mcc/'];
traj.t = tsample;
[traj.Mpc, traj.Mpv] = load_Mpc_Mpv(modeldir, tsample);


load('Rp_ipfit2.mat');
Rp_t = interp1(Rp_ipfit.t, Rp_ipfit.value, t, 'linear', 'extrap');

fit_Lp = load('fit_Lp.mat').fit_Lp;
Lp_t = fit_Lp(t);


% Load data
include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
        'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
vsignals = get_vobjcsignals(shot, [], [], include_coils);
vts = timeseries(vsignals.sigs, vsignals.times);
vts = resample(vts, tsample);
u = double(vts.Data);


load('coils_greybox.mat') % True Data
[~,k] = min(abs(coils.t - tstart));
ic0 = zeros(circ.ncx,1);
ic0(circ.ikeep) = coils.ic(:,k);
iv0 = coils.iv(:,k);
ip0 = coils.ip(k);
x0 = [ic0; iv0; ip0];

enforce_stability = 0;


icsignals = get_icsignals(shot, [], [], include_coils);
ipsignals = mds_fetch_signal(shot, 'efit01', '.RESULTS.AEQDSK:IPMEAS');
icts = timeseries(icsignals.sigs,icsignals.times);      
ipts = timeseries(ipsignals.sigs,ipsignals.times);      
icts = resample(icts,tsample);   
ipts = resample(ipts,tsample);   
Tsmooth = 10;  % [ms]
nsmooth = floor(Tsmooth/1000/Ts);
ic = smoothdata(icts.Data, 1, 'movmean', nsmooth);
ip = smoothdata(ipts.Data, 1, 'movmean', nsmooth);
icdot = gradient(ic', Ts)';
ipdot = gradient(ip', Ts)';
icdot = smoothdata(icdot,1,'movmean',nsmooth);
ipdot = smoothdata(ipdot,1,'movmean',nsmooth);



%%
% Simulate
x = x0;
iv = iv0;
sim = zeros(N, length(x0));
iv_sim = zeros(N,length(iv0));

for i = 1:N
  ui = u(i,:)';   
  Rp = Rp_t(i) * 1.3;
  Lp = Lp_t(i);
  Mpc = squeeze( interp1(traj.t, traj.Mpc, t(i), 'linear', 'extrap'));
  Mpv = squeeze( interp1(traj.t, traj.Mpv, t(i), 'linear', 'extrap'));
  
  [x, Ad, Bd, C, D, iv] = full_dynamics_vessmod(x, ui, Mxx, Rxx, Rp, Lp, Mpc, Mpv,...
  Ts, circ, enforce_stability, iv, icdot(i,:)', ipdot(i,:)');
  xsim(i,:) = x;
  iv_sim(i,:) = iv;
end

%%

% Plots

figure
ax(1) = subplot(131);
hold on
plot(t, xsim(:,circ.iicx), 'b')
plot(coils.t, coils.ic, '--r')
title([num2str(shot) ' Coil Currents'], 'fontsize', 14)
ylabel('[A]', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
mylegend({'True', 'Simulated'}, {'--','-'}, [], {'r','b'}, [], 'Northeast', 14);


ax(2) = subplot(132);
hold on
% plot(t,xsim(:,circ.iivx),'b')
plot(coils.t, coils.iv,'--r')
plot(t, iv_sim, 'b')
xlim([0 0.9])
title([num2str(shot) ' Vessel Currents'], 'fontsize', 14)
ylabel('[A]', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
mylegend({'True', 'Simulated'}, {'--','-'}, [], {'r','b'}, [], 'Northeast', 14);


ax(3) = subplot(133);
hold on
plot(t,xsim(:,circ.iipx), 'b')
plot(coils.t, coils.ip, '--r')
title('Ip')
linkaxes(ax, 'x')
xlim([0 1])
set(gcf, 'Position', [739 508 600 587])


%%
% ipdot = gradient(xsim(:,circ.iipx)');
% iptruedot = gradient(interp1(coils.t, coils.ip, t));
% dipdot = ipdot-iptruedot;
% 
% eps = 1e-10;
% Rp_t = Rp_t + eps*dipdot;
% 
% thresh = 1e-9;
% Rp_t(Rp_t < thresh) = thresh;
% 
% Rp_t = smooth(Rp_t)

% mean(Rp_t)  
% mean(dip)


% Rp_ipfit.t = t;
% Rp_ipfit.value = Rp_t;
% save('Rp_ipfit2','Rp_ipfit')


