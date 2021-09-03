clear; clc; close all;


icoil = 1;

% shot = 204966;
% build_dir = '/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/raw_data/pertbuilds_val/';
% targs = load([build_dir 'val_response_' num2str(shot) '.mat']).targs;


shot = 204966;
build_dir = '/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/raw_data/pertbuilds_val/';
targs = load([build_dir 'val_response_' num2str(shot) '.mat']).targs;



load('nstxu_obj_config2016_6565.mat')
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;

t = targs.actualtime;
nsamples = length(t);


% filter using gamma
thresh = 20;
gamma_smooth = smoothdata(targs.gamma, 'movmedian', 5);       
e = targs.gamma - gamma_smooth;  
ibad = abs(e) > thresh;  
igood = find(~ibad);
ibad = find(ibad);

% smooth response
y = squeeze(targs.dpsidix(:,:,1));
y = interp1(t(igood), y(igood,:), t, 'linear', 'extrap');
dpsidix_smooth = smoothdata(y, 1, 'movmedian', 20);
dpsidix_smooth = reshape(dpsidix_smooth, nsamples, 65, 65);

dpsidix = reshape(targs.dpsidix(:,:,icoil), nsamples, 65, 65);


figure
ax(1) = subplot(211);
hold on
grid on
plot(t, reshape(dpsidix, nsamples, []));
ax(2) = subplot(212);
hold on
grid on
plot(t, reshape(dpsidix_smooth, nsamples, []));
linkaxes(ax)
ylim([-6 1]*1e-5)


%% DPSIDIX

figure(1)
subplot(311)
plot(t, targs.ip)
xl(1) = xline(0);

subplot(312)
plot(t, targs.coil_currents)
xl(2) = xline(0);

subplot(313)
plot(t, targs.gamma)
xl(3) = xline(0);
ylim([0 200])

set(gcf, 'Position', [832 793 534 532])

figure(2)
set(gcf, 'Position', [727 210 823 509])

for i = 1:length(t)
  
  figure(2)
  subplot(131)
  [~, cs] = contourf(rg, zg, squeeze(dpsidix(i,:,:)), 20);  
  title(['t = ' num2str(t(i))], 'fontsize', 14)
  colorbar
   
  subplot(132)
  contourf(rg, zg, squeeze(dpsidix_smooth(i,:,:)), cs.LevelList);  
  title(['t = ' num2str(t(i))], 'fontsize', 14)
  colorbar
      
  subplot(133)
  cla
  hold on
  contour(rg, zg, squeeze(targs.psirz(i,:,:))', 20, 'k')
  scatter(targs.rcur(i), targs.zcur(i), 100, 'r', 'filled')
  
  xl(1).Value = t(i);
  xl(2).Value = t(i);
  xl(3).Value = t(i);
  
  pause(0.5)
  drawnow
end






















































