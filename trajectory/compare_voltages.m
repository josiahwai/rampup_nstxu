
sim = load('sim_inputs204660_smoothed.mat').sim_inputs;
usim = sim.Uhat';
tsim = sim.tspan(1:end-1);

include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
        'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
ps = get_vobjcsignals(204660, [], [], include_coils);
u = ps.sigs;
t = ps.times;

load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);

u = ps.sigs(:,circ.ikeep);
usim = usim(:,circ.ikeep);

figure
for i = 1:circ.ncx_keep
  subplot(2,4,i)
  hold on
  plot(t, u(:,i), '--r')
  plot(tsim, usim(:,i), 'b')
  legend('Power Supplies', 'MPC Sim')
  xlim([0 1])
  title(circ.keep_coils{i})
  
end
set(gcf,'Position', [392 93 1115 651])















