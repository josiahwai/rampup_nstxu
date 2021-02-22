ccc
RAMPROOT = getenv('RAMPROOT');
load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);


% Timing
traj = load('traj_fit.mat').traj;
tstart = traj.tspan(1);
tend = traj.tspan(end);
Ts = mean(diff(traj.tspan));
shot = 204660;
tsample = tstart:Ts:tend;
N = length(tsample);
t = tsample;

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


%%
% Simulate
[xsim] = state_dynamics(traj.Ad, traj.Bd, x0, N, u');
xsim = xsim(:,1:N)';


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


%%

function [x_all, x_tf] = state_dynamics(Alist, Blist, x_t0, N, u_all)
  x = x_t0;
  x_all = zeros(length(x_t0), N+1);
  x_all(:,1) = x_t0;
  
  for i = 1:N
    A = squeeze(Alist(i,:,:));
    B = squeeze(Blist(i,:,:));    
    x = A*x + B*u_all(:,i);
    x_all(:,i+1) = x;      
  end
  
  x_tf = x_all(:,end);
end










