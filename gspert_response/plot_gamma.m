figure
subplot(211)
hold on
title(['Ip ' num2str(targs.shot)], 'fontsize', 16)
plot(targs.actualtime, targs.ip/1e3)
ylabel('kA', 'fontsize', 14)
subplot(212)
hold on
title('Vertical growth Rate \gamma', 'fontsize', 16)
plot(targs.actualtime, targs.gamma)
xlabel('Time [s]', 'fontsize', 14)
ylabel('Hz', 'fontsize', 14)
ylim([-1 200])


