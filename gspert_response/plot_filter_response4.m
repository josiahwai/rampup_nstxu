clear all; clc; close all;

% shot = 203324;
shot = 203008;
build_dir = '/Users/jwai/Research/rampup_nstxu/gspert_response/old/';
fn = [build_dir 'train_response_' num2str(shot) '.mat'];
targs = load(fn).targs;
load('nstxu_obj_config2016_6565.mat')
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;

t = targs.actualtime;
nsamples = length(t);

if ~isfield(targs, 'dpsidix')
  mpp = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit.build_inputs.tok_data_struct.mpp;
  for i = 1:nsamples
    targs.dpsidix(i,:,:) = mpp * squeeze(targs.dcphidix(i,:,:));
  end
  save(fn, 'targs')
end

icoil = 1;


y = squeeze(targs.dpsidix(:,:,icoil));

figure(1)
plot(t, y)

y = reshape(y, [], 65, 65);
yinit = y;

for iter = 1:3
  ysmooth = smoothdata(y, 1, 'movmedian', 13);
  ysmooth = smoothdata(ysmooth, 1, 'gaussian', 3);
  for i = 1:nsamples
    dy = squeeze(y(i,:,:) - ysmooth(i,:,:));
    e(i) = norm(dy);
  end
  e = e / norm(squeeze(mean(y,1)));
  iuse = find(e < 1);
  y = interp1(t(iuse), y(iuse,:,:), t, 'previous', 'extrap');
  figure
  hold on
  plot(t,e)
  yline(1)
end

figure
y = smoothdata(y, 1, 'gaussian', 5);
plot(t, reshape(y,[],65*65))

%%
figure(1)
ylim([-1 1]*2e-4)
hold on
i = 62;
xline(i);

figure
hold on
contourf(squeeze(yinit(i,:,:)), 20)
colorbar


figure
hold on
dcphi = reshape(targs.dcphidix(i,:,icoil), 65, 65);
contourf(dcphi, 20)
colorbar



figure
plot(t, squeeze(targs.dcphidix(:,:,icoil)))
%%
% % first pass
% y = squeeze(targs.dpsidix(:,:,icoil));
% y = reshape(y, [], 65, 65);
% ysmooth = smoothdata(y, 1, 'movmedian', 13);
% ysmooth = smoothdata(ysmooth, 1, 'gaussian', 3);
% for i = 1:nsamples
%   dy = squeeze(y(i,:,:) - ysmooth(i,:,:));
%   e(i) = norm(dy);
% end
% e = e / norm(squeeze(mean(y,1)));
% iuse = find(e < 1);
% y = interp1(t(iuse), y(iuse,:,:), t, 'previous', 'extrap');
% figure
% hold on
% plot(t,e)
% yline(1)
% 
% % second pass
% ysmooth = smoothdata(y, 1, 'movmedian', 13);
% ysmooth = smoothdata(ysmooth, 1, 'gaussian', 3);
% for i = 1:nsamples
%   dy = squeeze(y(i,:,:) - ysmooth(i,:,:));
%   e(i) = norm(dy);
% end
% e = e / norm(squeeze(mean(y,1)));
% iuse = find(e < 1);
% y = interp1(t(iuse), y(iuse,:,:), t, 'previous', 'extrap');
% figure
% hold on
% plot(t,e)
% yline(1)
% 
% % third pass
% ysmooth = smoothdata(y, 1, 'movmedian', 13);
% ysmooth = smoothdata(ysmooth, 1, 'gaussian', 3);
% for i = 1:nsamples
%   dy = squeeze(y(i,:,:) - ysmooth(i,:,:));
%   e(i) = norm(dy);
% end
% e = e / norm(squeeze(mean(y,1)));
% iuse = find(e < 1);
% y = interp1(t(iuse), y(iuse,:,:), t, 'previous', 'extrap');
% figure
% hold on
% plot(t,e)
% yline(1)
% 
% % third pass
% ysmooth = smoothdata(y, 1, 'movmedian', 13);
% ysmooth = smoothdata(ysmooth, 1, 'gaussian', 3);
% for i = 1:nsamples
%   dy = squeeze(y(i,:,:) - ysmooth(i,:,:));
%   e(i) = norm(dy);
% end
% e = e / norm(squeeze(mean(y,1)));
% iuse = find(e < 1);
% y = interp1(t(iuse), y(iuse,:,:), t, 'previous', 'extrap');
% figure
% hold on
% plot(t,e)
% yline(1)
% 
% % final
% tuse = t(iuse);
% yuse = reshape(y(iuse,:,:), [], 65*65);
% figure
% plot(tuse, yuse)


% figure
% yuse = reshape(y(iuse,:,:), [], 65*65);
% plot(t(iuse), yuse)
% 
% y = interp1(t(iuse), y(iuse,:,:), t, 'previous', 'extrap');
% ysmooth = smoothdata(y, 1, 'gaussian', 5);
% ysmooth = reshape(ysmooth, [], 65*65);
% 
% figure
% plot(t,ysmooth)
% 
% y = interp1(t(iuse), y(iuse,:,:), t, 'previous', 'extrap');




% % filter using gamma
% thresh = 20;
% gamma_smooth = smoothdata(targs.gamma, 'movmedian', 5);
% e = targs.gamma - gamma_smooth;
% ibad = abs(e) > thresh;
% 
% s = vecnorm(squeeze(targs.dpsidix(:,:,icoil)), 2, 2);
% ibad = ibad' | (s > 10*median(s));
% 
% igood = find(~ibad);
% ibad = find(ibad);
% 
% 
% % smooth response
% y = squeeze(targs.dpsidix(:,:,icoil));
% y = interp1(t(igood), y(igood,:), t, 'previous', 'extrap');
% dpsidix_smooth = smoothdata(y, 1, 'movmedian', 13);
% dpsidix_smooth = smoothdata(dpsidix_smooth, 1, 'gaussian', 3);
% targs.dpsidix_smooth(:,:,icoil) = dpsidix_smooth;
% 
% dpsidix = reshape(targs.dpsidix(:,:,icoil), nsamples, 65, 65);
% 
% 
% figure
% hold on
% sgtitle(['Coil ' num2str(icoil)], 'fontsize', 14)
% set(gcf, 'Position', [1000 696 560 642])
% ax(1) = subplot(311);
% hold on
% grid on
% plot(t,targs.gamma)
% ylim([0 200])
% ax(2) = subplot(312);
% hold on
% grid on
% plot(t, reshape(dpsidix, nsamples, []));
% yyaxis right
% % plot(t,e,'linewidth',2)
% ax(3) = subplot(313);
% hold on
% grid on
% plot(t, dpsidix_smooth);
% linkaxes(ax, 'x')
% drawnow
















































