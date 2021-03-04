ccc
eqdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk_import/';
load('nstxu_obj_config2016_6565.mat')

Te = mds_fetch_signal(204660, 'activespec', [], '.MPTS.OUTPUT_DATA.BEST:FIT_TE');
[~,i] = rmoutliers(mean(Te.sigs));
Te.sigs = Te.sigs(:,~i);
Te_avg = mean(Te.sigs, 2);
Te_res_list = Te_avg * 1e3 * 0.7;

times = (30:10:960) / 1e3;
for i = 1:length(times)
  t = times(i);
  t_ms = times(i) * 1e3;
  
  Te_res = interp1(Te.times, Te_res_list, t);
  eq = load([eqdir 'eq204660_' num2str(t_ms) '.mat']).eq;
    
  Zeff_res = 1.5;
  li_res = 0.5;
  idxpl = find(eq.jphi(:)~=0);
  a0 = (max(tok_data_struct.rgg(idxpl))-min(tok_data_struct.rgg(idxpl)))/2;
  b0 = (max(tok_data_struct.zgg(idxpl))-min(tok_data_struct.zgg(idxpl)))/2;
  if(a0==b0)   % handle "point" distributions of current
    kap0=1;
  else
    kap0 =  b0/a0;
  end
  
  [taup(i),Lp(i),Resp(i),eta(i)] = calc_LR(eq.rcur,a0,kap0,li_res,Te_res,Zeff_res);
end

load('Rp_ipfit2.mat')

figure
semilogy(Rp_ipfit.t, Rp_ipfit.value, 'linewidth', 1.5)
hold on
semilogy(times, Resp, 'linewidth', 1.5)
xlim([0 0.9])   
ylim([3e-7 1e-4])
xlabel('Time [s]', 'fontsize', 14)
ylabel('Resistance [\Omega]', 'fontsize', 14)
title('Plasma resistance 204660', 'fontsize', 16)
legend('Fitted from coil dynamics', 'Calculated from Te', 'fontsize', 14)
grid on


