% Use experimental coil and vessel current data to calculate the
% time-dependent plasma resistivity

clear all; clc; close all
ROOT = getenv('RAMPROOT');

% ========
% Settings
% ========
% shot = 204960;
% shot = 204069;
% shot = 204660;
% shot = 203708;
% shot = 203172;

shotlist = load('/p/nstxusr/nstx-users/jwai/nstxu-nns/data/matlab/train_test_val_split.mat').split.valshots;

saveplot = 0;
saveit = 1;

%%
for ishot = 1:length(shotlist)
  
  try
  shot = shotlist(ishot)
  

  load('nstxu_obj_config2016_6565.mat');
  circ = nstxu2016_circ(tok_data_struct);
  
  tree = 'EFIT01';
  tokamak = 'nstxu';
  server = 'skylark.pppl.gov:8501';
  opts.cache_it = 1;
  opts.cache_dir = [ROOT '/fetch/cache/'];
  
  eqs = fetch_eqs_nstxu(shot, 'all', tree, tokamak, server, opts);
  
  i = eqs.time(:)' <= 0 | gradient([eqs.gdata(:).cpasma]) < -1e5 | [eqs.gdata(:).cpasma] < 1e3;
  eqs.adata(i) = [];
  eqs.gdata(i) = [];
  eqs.tms(i) = [];
  eqs.time(i) = [];
  
  t = eqs.time;
  dt = mean(diff(t));
  neq = length(t);
  
  icx = [eqs.gdata(:).icx];
  ivx = [eqs.gdata(:).ivx];
  ip = [eqs.gdata(:).cpasma];
  I = [icx; ivx; ip];
  Idot = gradient(I, dt);
  
  Idot = smoothdata(Idot);
  
  %%
  % close all
  
  % figure
  % k = 1;
  % hold on
  % grid on
  % % scatter(t,Idot(k,:))
  % plot(t,-Idot(k,:), '-b')
  % % plot(t,I(k,:), '-b')
  % xlim([0 0.1])
  %
  % yyaxis right
  % k = 54;
  % hold on
  % % scatter(t,Idot(k,:))
  % plot(t,Idot(k,:), '-r')
  % % plot(t,I(k,:), '-r')
  % xlim([0 0.1])
  %%
  
  % plot(t,I(end,:))
  
  % k = 8;
  % figure
  % hold on
  % plot(t, Idot(k,:))
  % scatter(t, smooth(Idot(k,:)))
  
  
  struct_to_ws(tok_data_struct);
  mpp = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit.build_inputs.tok_data_struct.mpp;
  
  clear Lp Rp A
  
  for i = 1:neq
    pcurrt = eqs.gdata(i).pcurrt(:);
    mcIp = mpc' * pcurrt / sum(pcurrt);
    mvIp = mpv' * pcurrt / sum(pcurrt);
    Lp(i) = pcurrt' * mpp * pcurrt / sum(pcurrt)^2;
    
    M = circ.Pxx' * [mcIp; mvIp; Lp(i)];
    
    Rp(i) = -1 / ip(i) * M' * Idot(:,i);
    
    A(i) = polyarea(eqs.gdata(i).rbbbs, eqs.gdata(i).zbbbs);
    
    
    %   figure
    %   plot_eq(eqs.gdata(i))
    %   title([num2str(shot) ':' num2str(eqs.time(i))])
    %   contour(eqs.gdata(i).rg, eqs.gdata(i).zg, eqs.gdata(i).psizr, 15, 'color', [1 1 1] * 0.8)
    %   fn = [ROOT '/buildmodel/resistivity_calc/eq/' num2str(shot) '_' num2str(i) '.png'];
    %   set(gcf, 'Position', [668 456 255 414])
    %   saveas(gcf, fn)
  end
  
  minRp = min(Rp(Rp>0)); % enforce nonnegative resistance
  Rp(Rp<0) = minRp;
  eta = Rp.*A;
  
  A(A<=0) = min(A(A>0)); % enforce nonzero area
  
  %%
  % use measured loop voltage in fitting dynamics
  
  tree = 'efit01';
  tag = '.results.aeqdsk.vloopmhd';
  % tag = '.results.aeqdsk.vsurf';
  vloop_sigs = mds_fetch_signal(shot, tree, [], tag, 0);
  vloop = interp1(vloop_sigs.times, vloop_sigs.sigs, t)';
  
  ipdot = Idot(end,:);
  Rp2 = (vloop - Lp.*ipdot)./ip;
  % Rp2 = vloop./ip;
  
  minRp2 = min(Rp2(Rp2>0));  %  enforce nonnegative resistance
  Rp2(Rp2<0) = minRp2;
  eta2 = Rp2.*A;
  
  % figure
  % hold on
  % plot(t, Rp, '--r')
  % plot(t, Rp2, '-b')
  % title('Rp: dynamics fit with Vloop')
  
  % yyaxis right
  % hold on
  % scatter(t,Lp)
  % plot(t,Lp)
  % xlim([0 0.3])
  %%
  % use Te and spitzer resistivity
  
  tree = 'activespec';
  tag = '.mpts.output_data.best:fit_te';
  Te_sigs = mds_fetch_signal(shot, tree, [], tag, 0);
  
  % remove outliers (sensors)
  Te = rmoutliers(Te_sigs.sigs', 'ThresholdFactor', 10)';
  % plot(Te)
  
  % remove outliers (times)
  [Te,i] = rmoutliers(Te, 'ThresholdFactor', 10);
  Te_times = Te_sigs.times(~i);
  
  Te = interp1(Te_times, Te, t, 'nearest', 'extrap');
  
  % figure
  % plot(t,Te)
  
  Te_avg = mean(Te');
  % f = fit(t, Te_avg(:), 'smoothingspline','SmoothingParam',1-1e-6);
  % plot(f, t, Te_avg)
  
  
  % Constants:
  mu0 = 0.4*pi;
  Zeff=1;
  lnlam = 15;
  
  %  Resistivity, Resistance:
  %	From Chen, p. 183:
  %    	    eta = 5.2e-5*Zeff*lnlam/(Te^1.5); %Pl.restvty ohm-m Spitzer
  %	From Hutchinson, p. 19:
  %	    lnlam = 31 - log(sqrt(ne)/Te);  %valid for Te>10 eV
  eta_spitzer = 5.2e-5*Zeff*lnlam./(Te_avg.^1.5);      %Pl.restvty ohm-m Spitzer
  Rp_spitzer = eta_spitzer./A;
  
  
  
  %%
  c = Rp / Rp_spitzer;
  
  figure
  subplot(211)
  hold on
  grid on
  % scatter(t,Rp)
  plot(t,Rp)
  plot(t,Rp2)
  plot(t, c*Rp_spitzer)
  title(['Rp ' num2str(shot)])
  legend('Rp', 'Rp2', 'Rp spitzer')
  
  subplot(212)
  hold on
  grid on
  plot(t, eta)
  plot(t, eta2)
  % scatter(t, eta)
  plot(t, c*eta_spitzer)
  title(['eta ' num2str(shot)])
  legend('eta', 'eta2', 'eta spitzer')
  
  %%
  % c = Rp / Rp_spitzer;
  %
  % figure
  % subplot(211)
  % semilogy(t,Rp)
  % hold on
  % semilogy(t,Rp2)
  % grid on
  % semilogy(t, c*Rp_spitzer)
  % title('Rp')
  % legend('Rp', 'Rp2', 'Rp spitzer')
  %
  % subplot(212)
  % semilogy(t, eta)
  % hold on
  % semilogy(t, eta2)
  % grid on
  % semilogy(t, c*eta_spitzer)
  % title('Eta')
  % legend('eta', 'eta2', 'eta spitzer')
  
  
  if saveplot
    fn = [ROOT 'buildmodel/resistivity_calc/figs/res' num2str(shot) '.png'];
    saveas(gcf, fn)
  end
  
  if saveit
    res = variables2struct(t, A, Lp, Rp, Rp2, Rp_spitzer, eta, eta2, ...
      eta_spitzer, Te_avg, Te);
    
    fn = [ROOT 'buildmodel/resistivity_calc/res/res' num2str(shot) '.mat'];
    save(fn, 'res')
  end
  catch
    warning(['erred on shot ' num2str(shot)])
  end
end


























