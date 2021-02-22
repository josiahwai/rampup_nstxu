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


load('Rp_ipfit.mat');
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

%%
% Simulate
x = x0;
xsim = zeros(N, length(x0));
for i = 1:N
  ui = u(i,:)';   
  Rp = Rp_t(i);
  Lp = Lp_t(i);
  Mpc = squeeze( interp1(traj.t, traj.Mpc, t(i), 'linear', 'extrap'));
  Mpv = squeeze( interp1(traj.t, traj.Mpv, t(i), 'linear', 'extrap'));
  
  [x, Ad, Bd, C, D] = full_dynamics(x, ui, Mxx, Rxx, Rp, Lp, Mpc, Mpv, Ts, circ, enforce_stability);
  xsim(i,:) = x;
end



% Plots

figure
ax(1) = subplot(311);
hold on
plot(t, xsim(:,circ.iicx), 'b')
plot(coils.t, coils.ic, '--r')
title([num2str(shot) ' Coil Currents'], 'fontsize', 14)
ylabel('[A]', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
mylegend({'True', 'Simulated'}, {'--','-'}, [], {'r','b'}, [], 'Northeast', 14);


ax(2) = subplot(312);
hold on
plot(t,xsim(:,circ.iivx),'b')
plot(coils.t, coils.iv,'--r')
xlim([0 0.9])
title([num2str(shot) ' Vessel Currents'], 'fontsize', 14)
ylabel('[A]', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
mylegend({'True', 'Simulated'}, {'--','-'}, [], {'r','b'}, [], 'Northeast', 14);


ax(3) = subplot(313);
hold on
plot(t,xsim(:,circ.iipx), 'b')
plot(coils.t, coils.ip, '--r')
title('Ip')
linkaxes(ax, 'x')
xlim([0 1])
set(gcf, 'Position', [739 508 600 587])













