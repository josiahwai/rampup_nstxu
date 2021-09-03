clear; clc; close all;

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

%%
for icoil = 1:54
  icoil
  
  % filter using gamma
  thresh = 20;
  gamma_smooth = smoothdata(targs.gamma, 'movmedian', 5);       
  e = targs.gamma - gamma_smooth;  
  ibad = abs(e) > thresh;  

  s = vecnorm(squeeze(targs.dpsidix(:,:,icoil)), 2, 2);
  ibad = ibad' | (s > 10*median(s));

  igood = find(~ibad);
  ibad = find(ibad);


  % smooth response
  y = squeeze(targs.dpsidix(:,:,icoil));
  y = interp1(t(igood), y(igood,:), t, 'previous', 'extrap');
  dpsidix_smooth = smoothdata(y, 1, 'movmedian', 13);
  dpsidix_smooth = smoothdata(dpsidix_smooth, 1, 'gaussian', 3);
  targs.dpsidix_smooth(:,:,icoil) = dpsidix_smooth;

  dpsidix = reshape(targs.dpsidix(:,:,icoil), nsamples, 65, 65);

  
%   figure
%   hold on
%   sgtitle(['Coil ' num2str(icoil)], 'fontsize', 14)
%   set(gcf, 'Position', [1000 696 560 642])
%   ax(1) = subplot(311);
%   hold on
%   grid on
%   plot(t,targs.gamma)
%   ax(2) = subplot(312);
%   hold on
%   grid on
%   plot(t, reshape(dpsidix, nsamples, []));
%   ax(3) = subplot(313);
%   hold on
%   grid on
%   plot(t, dpsidix_smooth);
%   linkaxes(ax, 'x')
%   drawnow
end
%%
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);
P = circ.Pxx(1:end-1, 1:end-1);

mcc = tok_data_struct.mcc;
mvv = tok_data_struct.mvv;
mcv = tok_data_struct.mcv;
mpc = tok_data_struct.mpc;
mpv = tok_data_struct.mpv;
rxx = P' * diag([tok_data_struct.resc; tok_data_struct.resv]) * P;
mmat = P' * [mcc mcv; mcv' mvv] * P;

mpp = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit.build_inputs.tok_data_struct.mpp;
mppi = inv(mpp);

gamma = [];
for i = 1:nsamples  

  dpsidix = squeeze(targs.dpsidix_smooth(i,:,:));
  dpsidix(:,circ.iipx) = [];
  dcphidix = mppi * dpsidix;
  xmat = P' * [mpc'; mpv'] * dcphidix;
  
  xmat = squeeze(targs.xmat(i,:,:));
  amat = -inv(mmat + xmat)*rxx;
  e = real(esort(eig(amat)));
  gamma(i) = e(1);
end


figure
hold on
plot(t(igood), gamma(igood), '-b')
plot(t(igood), targs.gamma(igood), '--r')
% plot(t, targs.gamma)
ylim([0 500])


















































