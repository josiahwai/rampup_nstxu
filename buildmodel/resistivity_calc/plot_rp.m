clear all; clc; close all

shots = [203172, 203708, 204069, 204660, 204960];
shot0 = 204960;

load(['res' num2str(shot0) '.mat'])

figure(1)
hold on
plot(res.t, res.A)
xlabel('Time [s]', 'fontsize', 18)
ylabel('Plasma Area [m^2]', 'fontsize', 18)
xlim([0 0.5])
title(shot0, 'fontsize', 18)


figure(2)
hold on
i = res.t > 0.07;
c = res.Rp(i) / res.Rp_spitzer(i);
plot(res.t, res.Rp, 'r')
plot(res.t, res.Rp2, 'b')
plot(res.t, c*res.Rp_spitzer, 'k')
xlabel('Time [s]', 'fontsize', 18)
ylabel('Resistance [Ohm]', 'fontsize', 18)
xlim([0 0.5])
legend('Fit', 'Fit with Vloop', 'Spitzer Rp(Te)', 'fontsize', 18)
title(shot0, 'fontsize', 18)


figure(3)
hold on
plot(res.t, res.Te_avg, 'r', 'linewidth', 3)
plot(res.t, res.Te, 'color', [1 1 1] * 0.8)
plot(res.t, res.Te_avg, 'r', 'linewidth', 3)
xlabel('Time [s]', 'fontsize', 18)
ylabel('Thomson Te [eV]', 'fontsize', 18)
legend('channel average', 'fontsize', 18)
title(shot0, 'fontsize', 18)


figure(4)
hold on
for ishot = 1:length(shots)
  
  shot = shots(ishot);
  load(['res' num2str(shot) '.mat'])    
  plot(res.t, res.eta)
  xlabel('Time [s]', 'fontsize', 18)
  ylabel('Resistivity [Ohm/m^2]', 'fontsize', 18)
  xlim([0 0.5])
  
end
legend(cellstr(num2str(shots')), 'fontsize', 18, 'location', 'best')


figure(5)
hold on
for ishot = 1:length(shots)
  
  shot = shots(ishot);
  load(['res' num2str(shot) '.mat'])    
  plot(res.t, res.Lp)
  xlabel('Time [s]', 'fontsize', 18)
  ylabel('Plasma Inductance [H]', 'fontsize', 18)
  xlim([0 0.5])
  
end
legend(cellstr(num2str(shots')), 'fontsize', 18, 'location', 'best')







