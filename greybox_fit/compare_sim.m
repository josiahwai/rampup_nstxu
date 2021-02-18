% clear all
% traj = load('sim_inputs204660_smoothed.mat').sim_inputs.traj;
% load('nstxu_obj_config2016_6565.mat')
% circ = nstxu2016_circ(tok_data_struct);
% times = (30:10:940)/1000;
% plotit = 1;

% coils = load('coils_greybox.mat').coils;
% ps_voltages  = pinv(circ.Pcc_keep) * coils.v;
% ps_voltages = interp1(coils.t, ps_voltages', times);
% uall = smoothdata(ps_voltages, 1, 'movmean', 13);

% addpath(genpath('/usr/local/mdsplus/'))
% include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
%         'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
%       
% coils = get_vobjcsignals(204660, [], [], include_coils);
% uall = interp1(coils.times, coils.sigs, times);

function xfit = compare_sim(uall, traj, circ, times, plotit)

Ts = mean(diff(times));
lstar_invs = interp1(traj.tspan, traj.lstari, times);
file_args = {Ts, circ, lstar_invs, times};
N = length(times);

rcc = diag(circ.Pcc' * diag(tok_data_struct.resc) * circ.Pcc);
rvv = diag(circ.Pvv' * diag(tok_data_struct.resv) * circ.Pvv);
rp_t = load('rp_t').rp_t;
voltage_scale = ones(circ.nu,1);

xtarg = traj.x;
xtarg = interp1(traj.tspan, xtarg, times)';
% xtarg = smoothdata(xtarg, 2, 'movmean', 13);
x0 = xtarg(:,1);


xfit(:,1) = x0;
x = x0;
for i = 1:N-1 
  t = times(i);
  u = uall(i,:)';
  x = nl_grey_nstxu_modelAB(t, x, u, rcc, rvv, rp_t, [], voltage_scale, file_args);  
  xfit(:,i+1) = x;
end

if plotit
  figure
  subplot(411)
  hold on
  plot(times, xfit(circ.iipx, :), 'b')
  plot(times, xtarg(circ.iipx,:), '--r')
  set(gcf, 'Position', [680 248 467 730])
  xlabel('Time [s]', 'fontsize', 14)
  ylabel('Ip', 'fontsize', 14)
  legend('Fit', 'Target', 'fontsize', 14)
  
  subplot(412)
  hold on
  plot(times, xfit(circ.iicx, :), 'b')
  plot(times, xtarg(circ.iicx,:), '--r')
  ylabel('Coil Currents', 'fontsize', 14)
  legend('Fit', 'Target', 'fontsize', 14)
  
  subplot(413)
  hold on
  plot(times, xfit(circ.iivx, :), 'b')
  plot(times, xtarg(circ.iivx,:), '--r')
  ylabel('Vessel Currents', 'fontsize', 14)
  legend('Fit', 'Target', 'fontsize', 14)
  
  subplot(414)
  hold on
  plot(times, uall')
  title('Voltages')
end































