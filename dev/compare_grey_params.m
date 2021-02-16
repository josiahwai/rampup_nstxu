ccc
load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);
cclabels = categorical(circ.ccnames);
vvlabels = categorical(circ.vvnames);
t = (30:10:940) / 1000;
Nt = length(t);

i = 1;
load('grey_init_sys.mat')
[rcc(:,i), rvv(:,i), rp_t(:,i), lp_t(:,i), voltage_scale(:,i)] = unload_grey(grey_init_sys, Nt, 1);

i = 2;
load('grey_sys0.mat')
[rcc(:,i), rvv(:,i), rp_t(:,i), lp_t(:,i), voltage_scale(:,i)] = unload_grey(grey_sys, Nt, 1);

i = 3;
load('grey_sys1.mat')
[rcc(:,i), rvv(:,i), rp_t(:,i), lp_t(:,i), voltage_scale(:,i)] = unload_grey(grey_sys, Nt, -1);

figure
subplot(221)
bar(cclabels, rcc)
ylim([0 0.1])
title('Coil Resistance')
legend('Original', 'Fit (t=0-230ms)', 'Fit (t=230-950ms)', 'fontsize', 14)

subplot(222)
bar(vvlabels, rvv)
ylim([0 .005])
title('Vessel Resistance')
legend('Original', 'Fit (t=0-230ms)', 'Fit (t=230-950ms)', 'fontsize', 14)

subplot(223)
plot(t, rp_t, 'linewidth', 2)
title('Plasma Resistance')
xlabel('Time [s]')
legend('Original', 'Fit (t=0-230ms)', 'Fit (t=230-950ms)', 'fontsize', 14)

subplot(224)
bar(cclabels, voltage_scale)
ylim([0 4])
title('Voltage Scaling')
legend('Original', 'Fit (t=0-230ms)', 'Fit (t=230-950ms)', 'fontsize', 14)



function [rcc, rvv, rp_t, lp_t, voltage_scale] = unload_grey(grey_sys, Nt, side)  
  rcc = grey_sys.Parameters(1).Value;
  rvv = grey_sys.Parameters(2).Value;  
  voltage_scale = grey_sys.Parameters(5).Value;
  
  M = length(grey_sys.Parameters(3).Value);
  rp_t = nan(Nt, 1);
  lp_t = nan(Nt, 1);
  
  if side == 1
    rp_t(1:M) = grey_sys.Parameters(3).Value;
    lp_t(1:M) = grey_sys.Parameters(4).Value;
  elseif side == -1
    rp_t(end-M+1:end) = grey_sys.Parameters(3).Value;
    lp_t(end-M+1:end) = grey_sys.Parameters(4).Value;
  end
end
























