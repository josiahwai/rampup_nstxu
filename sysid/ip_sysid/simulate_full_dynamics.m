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


% Load parameters
vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);
Mxx = vac_sys.sysid_fits.Mxx;
Rxx = vac_sys.sysid_fits.Rxx;

sys = load([RAMPROOT 'buildmodel/built_models/std_coneqt/204660_300_sys.mat']).sys;
Mpc = sys.lstar(circ.iipx,circ.iicx);
Mpv = sys.lstar(circ.iipx,circ.iivx);
Lp = sys.lstar(circ.iipx, circ.iipx);
Rp = sys.rxx(end);

modeldir = [RAMPROOT 'buildmodel/built_models/std/'];
[Mpc_all, Mpv_all] = load_Mpc_Mpv(modeldir, tsample);

fit_Rp = load('fit_Rp5.mat').fit_Rp;
fit_Lp = load('fit_Lp.mat').fit_Lp;

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
  Rp = fit_Rp(t(i));
  Lp = fit_Lp(t(i));
  Mpc = Mpc_all(i,:) * 0;
  Mpv = Mpv_all(i,:) * 0;
  
  [x, Ad, Bd, C, D] = full_dynamics(x, ui, Mxx, Rxx, Rp, Lp, Mpc, Mpv, Ts, circ, enforce_stability);
  xsim(i,:) = x;
end


%% 
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













