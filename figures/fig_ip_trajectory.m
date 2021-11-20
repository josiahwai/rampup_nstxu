close all

figure
hold on
plot(t, icxhat(2:end,:), '-b')
icx_efit = [efit01_eqs.gdata(:).icx];
plot(t, icx_efit(2:end,:) , '--r')
for i = circ.iicx_keep(2:end)  
  text(t(end)+.002, icxhat(i,end), circ.ccnames{i}, 'fontsize', 14)
end
xlabel('Time [s]', 'fontsize', 16)
ylabel('Coil currents [A]', 'fontsize', 14)
l = mylegend({'Optimizer', 'Experiment'}, {'-', '--'}, [], {'b', 'r'});
l.FontSize = 14;
title('Coil trajectories', 'fontsize', 18)
set(gcf, 'Position', [669 326 505 314])
xlim([0 1.02])

%%

figure
hold on
plot(t,iphat/1e6,'--b', t,targets.ip/1e6,'-b')
ylabel('Ip [MA]', 'fontsize', 18)
yyaxis right
hold on
plot(t,icxhat(1,:)/1e3, 'r')
plot(t, icx_efit(1,:)/1e3, '--r')
ylabel('OH [kA]', 'fontsize', 18)
title('Ip trajectory', 'fontweight', 'bold', 'fontsize', 16)
xlabel('Time [s]', 'fontsize', 16)
l = mylegend({'Optimizer', 'Experiment'}, {'-', '--'}, [], {[1 1 1]*.02, [1 1 1]*0.2});
l.FontSize = 16;
l.Location = 'best';
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
set(gcf, 'Position', [680 723 487 255])