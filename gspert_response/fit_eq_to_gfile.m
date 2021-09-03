ccc

tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';

shot = 203324;
itime = 100;

eqs = read_eq(shot, 'all', tree, tokamak, server);
t = eqs.time(itime);
eq = eqs.gdata(itime);
coils = fetch_coilcurrents_nstxu(shot, t);

tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);



%%
% copy circuit information into eq
eq.ecturn = tok_data_struct.ecnturn;
eq.ecid   = ones(size(eq.ecturn));
eq.turnfc = tok_data_struct.fcnturn';
eq.fcturn = circ.fcfrac;
eq.fcid = circ.fccirc';

init = eq;

config = tok_data_struct;
config.max_iterations = 30;
config.constraints = 1;
config.nkn = 10;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

j = init.rbbbs ~= 0 & init.zbbbs ~=0;
spec.targets.rsep = init.rbbbs(j);
spec.targets.zsep = init.zbbbs(j);
spec.weights.sep = ones(length(spec.targets.rsep),1) * 1;
   
spec.targets.psibry = init.psibry;
spec.weights.psibry = 1;

spec.targets.zmaxis = init.zmaxis;
spec.weights.zmaxis = 1;

spec.targets.rmaxis = init.rmaxis;
spec.weights.

gs_configure
gs_initialize
gs_eq_analysis

spec.targets.li = interp1(tspan, traj.li, t);
spec.weights.li = 10;

spec.targets.betap = interp1(tspan, traj.betap, t);
spec.weights.betap = 10;

if t < 0.17 && t > 0.12
  spec.weights.betap = 0.1;
  spec.weights.li = 0.1;
  spec.targets.rbdef = 0.315;
  spec.targets.zbdef = 0.02;
  spec.weights.bdef = 100;
end

config.constraints = 1; % allow for scaling/peaking of profiles
config.pres0 = interp1(tspan, traj.pres, t)';
config.fpol0 = interp1(tspan, traj.fpol, t)';

spec.cccirc = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 ...
    11 12 13];

spec.limits.ic(1:tok_data_struct.nc,1) = [-20 0 0 -8 0 0 -13 -13 -13 ...
    -13 0 0 0 0 0 0 -24 -24 -24 -24 -13 -13 -13 -13 0 0 -8 0 0]'*1000;

spec.limits.ic(1:tok_data_struct.nc,2) = [20 15 0 13.5 15 15 8 8 8 8 ...
    13 13 13 13 13 13 0 0 0 0 8 8 8 8 15 15 13.5 0 15]'*1000;
     
% indices in x for: coil currents, vessel currents, plasma current
iic = 1:8;
iiv = 9:48;
iip = 49;
  
% Lock vessel current
[~,~,~,Pvv] = nstxu2020_vvcirc;
iv = x(iiv);
spec.locks.iv = Pvv*iv;

% Lock plasma current
spec.locks.cpasma = x(iip);

% Lock coil currents
% Some coil currents were turned off for campaign, so lock these to zero

removed_coils = {'PF1BU', 'PF1BL', 'PF1CU', 'PF1CL', 'PF4'};  
coils = {'OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', 'PF5', ...
  'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL'};

for i = 1:length(coils)
  iuse(i) = ~ismember(coils{i}, removed_coils);
end
ic = zeros(size(coils));
ic(iuse) = x(iic);
spec.locks.ic = ic(spec.cccirc);


































  build_inputs.tokamak = 'NSTXU';
  build_inputs.vacuum_objs = tok_data_struct;
  build_inputs.ichooseq = 4;
  build_inputs.irzresp_dynamic = 5;
  build_inputs.irzresp_output = 5;
  build_inputs.iplcirc = 1;
  build_inputs.cccirc = circ.cccirc(:);
  build_inputs.vvcirc = circ.vvcirc(:);
  build_inputs.vvgroup = circ.vvgroup(:);
  
  
  
  [~, t_idx] = min(abs(requesttimes - eqs.time));
  actualtimes = eqs.time(t_idx);
  
  targs = struct;
  isample = 0;
  ibad = 0;
  
  for i = 1:length(t_idx)
    
    requesttime = requesttimes(i);
    actualtime = actualtimes(i);
    
    time = actualtime;
    disp(['time_ms: ' num2str(time*1000)])
    
    try
      eq = eqs.gdata(t_idx(i));
      eq.ecturn = tok_data_struct.ecnturn;
      eq.ecid   = ones(size(eq.ecturn));
      eq.turnfc = tok_data_struct.fcnturn';
      eq.fcturn = circ.fcfrac;
      eq.fcid = circ.fccirc';
      
      build_inputs.equil_data = eq;
      
      sys = build_tokamak_system(build_inputs);
      delete('NSTXU_netlist.dat')
      
      P = circ.Pxx(1:end-1,1:end-1);
      xmatx = sys.xmatx(1:end-1,1:end-1);
      
      xmat = P'*xmatx*P;
      xmat = double(xmat);
      
      e = esort(eig(sys.amat(1:end-1,1:end-1)));
      gamma = real(e(1));
      gamma = double(gamma);
      disp(['gamma: ' num2str(gamma)])
      
      dcphidis = [sys.gspert_data.dcphidis sys.gspert_data.dcphidip(:)];
      dcphidix = dcphidis * circ.Pxx;
      dcphidix = double(dcphidix);
      
      % save to struct
      isample = isample+1;
      
      targs.xmat(isample,:,:) = xmat;
      targs.gamma(isample) = gamma;
      targs.dcphidix(isample,:,:) = dcphidix(:,circ.ikeep);
      targs.requesttime(isample) = requesttime;
      targs.actualtime(isample) = actualtime;
      targs.ip(isample) = eq.cpasma;
      targs.pprime(isample,:) = eq.pprime;
      targs.ffprim(isample,:) = eq.ffprim;
      targs.pres(isample,:) = eq.pres;
      targs.psirz(isample,:,:) = eq.psirz;
      targs.pcurrt(isample,:,:) = eq.pcurrt;
      targs.shot = shot;
      targs.rcur(isample) = sum(tok_data_struct.rgg(:).*eq.pcurrt(:)) / eq.cpasma;
      targs.zcur(isample) = sum(tok_data_struct.zgg(:).*eq.pcurrt(:)) / eq.cpasma;

      
    catch
      warning('Bad shot!')
      ibad = ibad+1;
      bad.shot(ibad) = shot;
      bad.time(ibad) = time;
    end
  end
  
  fn = [save_dir mode '_response_' num2str(shot) '.mat'];
  save(fn, 'targs', '-v7.3')

disp('Done!')




























































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
plot(y)

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
plot(reshape(y,[],65*65))

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
















































