clear; clc; close all;

saveit = 0;

mode = 'val';
build_dir = ['/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/raw_data/pertbuilds/geqdsk_' mode '/'];

d = dir(build_dir);
d(1:2) = [];

% geo files for evaluating growth rate
geom_dir = '/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/eval/';
M = load([geom_dir 'M.mat']).M;
Mpp = load([geom_dir 'Mpp.mat']).Mpp;
MpcMpv = load([geom_dir 'MpcMpv.mat']).MpcMpv;
Pxx = load([geom_dir 'Pxx.mat']).Pxx;
Rxx = load([geom_dir 'Rxx.mat']).Rxx;
Mppi = inv(Mpp);

%%
close all

% smoothing + plots
for ishot = 23:23 %1:length(d)
  
  % fn = [build_dir 'val_response_' num2str(shot) '.mat'];
  fn = [d(ishot).folder '/' d(ishot).name];
  
  targs = load(fn).targs;
  t = targs.time;
  nsamples = length(t);
  
  gamma_smooth = smoothdata(targs.gamma, 'movmedian', 15);
  e_gamma = abs(targs.gamma - gamma_smooth);
  igood = targs.gamma < 300 & e_gamma < 20 & targs.ip > 1e5;
  
  shot_is_good = max(t) > 0.3 & sum(igood)/length(igood) > 0.7;
  
  y = targs.dpsidix;
  ysm = smoothdata(y(igood,:,:), 'movmedian', 5);
  ysm = interp1(t(igood), ysm, t, 'linear', 'extrap');
  
  gamma_pred = t*0;
  
  for i = 1:nsamples
    yi = squeeze(ysm(i,:,1:end-1));
    X = Pxx' * MpcMpv' * Mppi * yi;
    amat = -inv(M+X)*Rxx;
    gamma_pred(i) = max(real(eig(amat)));
  end
  
  e_gamma = abs(gamma_smooth - gamma_pred);
  igood2 = e_gamma < 20; 
  igood = igood & igood2;
  
  for icoil = [1 2 40]
    figure
    hold on
    sgtitle(['Shot ' num2str(targs.shot) ' Coil ' num2str(icoil)], 'fontsize', 14)
    set(gcf, 'Position', [1048 106 560 642])
    ax(1) = subplot(311);
    hold on
    grid on
    plot(t, targs.gamma, 'linewidth', 1.5)
    plot(t, gamma_pred, 'linewidth', 1.5, 'linestyle', '--')
    % scatter(t, targs.gamma)
%     ylim([-10 min(max(targs.gamma), 200)])

    ax(2) = subplot(312);
    hold on
    grid on
    plot(t, reshape(y(:,:,icoil), nsamples, []));

    ax(3) = subplot(313);
    hold on
    grid on
    plot(t, reshape(ysm(:,:,icoil), nsamples, []));
    ylimits = ylim;

    linkaxes(ax, 'x')
    linkaxes(ax(2:3), 'y')
    ylim(ylimits)
    drawnow
  end
end

%%
icoil = 2;

ynames = {'ip', 'rcur', 'zcur', 'psibry', 'psimag', 'rmaxis', 'zmaxis'};
ny = length(ynames);
figure
for i = 1:ny
  subplot(ny,1,i)
  plot(targs.(ynames{i}))
  legend(ynames{i}, 'fontsize', 16)
  grid on 
  xlim([105 120])
end

figure
plot(reshape(y(:,:,icoil), nsamples, []));
grid on
xlim([105 120])


for isample = 110:115
  figure
  sgtitle(isample, 'fontweight', 'bold')
  subplot(121);
  contourf(squeeze(targs.psirz(isample,:,:))', 20)
  colorbar
  subplot(122);
  contourf(reshape(y(isample,:,icoil), 65, 65), 20)
  colorbar
  set(gcf, 'Position', [800 330 695 557])
end


load('nstxu_obj_config2016_6565.mat')
zlim = tok_data_struct.limdata(1,:);
rlim = tok_data_struct.limdata(2,:);

figure
hold on
plot(rlim, zlim, 'k', 'linewidth', 2)
for i = 1:20
  sz = 30;
  if ismember(i, 110:115)
    sz = 100;
  end
  plot(targs.rbbbs(i,:), targs.zbbbs(i,:))
  scatter(targs.rmaxis(i), targs.zmaxis(i), sz, 'filled')
  scatter(targs.rcur(i), targs.zcur(i), sz, 'filled')
end

min(rlim)

min(targs.rbbbs(i,:))

islimited = t*0;
for i = 1:nsamples
  if abs(min(targs.rbbbs(i,:)) - min(rlim)) < .01
    islimited(i) = 1;
  end
end

pts = find(islimited);

bar(islimited)










































