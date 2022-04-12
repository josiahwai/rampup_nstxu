

icx_efit = [efit01_eqs.gdata(:).icx];

figure
cmap = colormap('lines');
hold on

for i = [2 5 6 8 9 10 13]

  co = cmap(i,:);
  y1 = icxhat(i,:)/1e3;
  y2 = icx_efit(i,:)/1e3;

  plot(t, y1, 'linewidth', 2, 'color', co)
  plot(t, y2, '--', 'linewidth', 2, 'color', co)

  text(t(end)+.002, y1(end), circ.ccnames{i}, 'fontsize', 14, 'color', co)

end

xlabel('Time [s]', 'fontsize', 16)
ylabel('Coil currents [kA]', 'fontsize', 14)
l = mylegend({'Optimizer', 'Experiment'}, {'-', '--'}, [], {'k', 'k'});
l.FontSize = 14;
title('Coil trajectories', 'fontsize', 18)
set(gcf, 'Position', [669 326 505 314])
xlim([0 1.02])
