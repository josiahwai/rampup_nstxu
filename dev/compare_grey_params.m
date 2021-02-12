ccc
load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);
cclabels = categorical(circ.ccnames);
vvlabels = categorical(circ.vvnames);
t = (30:10:940) / 1000;
Nt = length(t);

i = 1;
load('grey_init_sys.mat')
[rcc(:,i), rvv(:,i), rp_t(:,i), lp_t(:,i), voltage_scale(:,i)] = unload_grey(grey_init_sys, Nt);

i = 2;
load('grey_sys1.mat')
[rcc(:,i), rvv(:,i), rp_t(:,i), lp_t(:,i), voltage_scale(:,i)] = unload_grey(grey_sys, Nt);

i = 3;
load('grey_sys2.mat')
[rcc(:,i), rvv(:,i), rp_t(:,i), lp_t(:,i), voltage_scale(:,i)] = unload_grey(grey_sys, Nt);

% i = 1;
% load('grey_init_sys.mat')
% [rcc{i}, rvv{i}, rp_t{i}, lp_t{i}, voltage_scale{i}] = unload_grey(grey_init_sys);
% 
% i = 2;
% load('grey_sys1.mat')
% [rcc{i}, rvv{i}, rp_t{i}, lp_t{i}, voltage_scale{i}] = unload_grey(grey_sys);
% 
% i = 3;
% load('grey_sys2.mat')
% [rcc{i}, rvv{i}, rp_t{i}, lp_t{i}, voltage_scale{i}] = unload_grey(grey_sys);
% 

figure
subplot(221)
bar(cclabels, rcc)
title('Coil Resistance')
legend('Original', 'Fit 1', 'Fit 2')

subplot(222)
bar(vvlabels, rvv)
ylim([0 .005])
title('Vessel Resistance')
legend('Original', 'Fit 1', 'Fit 2')

subplot(223)
plot(t, rp_t)
title('Plasma Resistance')
xlabel('Time [s]')
legend('Original', 'Fit 1', 'Fit 2')

subplot(224)
bar(cclabels, voltage_scale)
title('Voltage Scaling')
legend('Original', 'Fit 1', 'Fit 2')



function [rcc, rvv, rp_t, lp_t, voltage_scale] = unload_grey(grey_sys, Nt)  
  rcc = grey_sys.Parameters(1).Value;
  rvv = grey_sys.Parameters(2).Value;  
  voltage_scale = grey_sys.Parameters(5).Value;
  
  M = length(grey_sys.Parameters(3).Value);
  rp_t = nan(Nt, 1);
  lp_t = nan(Nt, 1);
  
  rp_t(end-M+1:end) = grey_sys.Parameters(3).Value;
  lp_t(end-M+1:end) = grey_sys.Parameters(4).Value;  
end
























