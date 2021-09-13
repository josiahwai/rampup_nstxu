% Use experimental coil and vessel current data to calculate the
% time-dependent plasma resistivity

clear all; clc; close all
ROOT = getenv('RAMPROOT');

% ========
% Settings
% ========
shot = 204960;
cache_dir = [ROOT '/fetch/cache/'];

%%
load('nstxu_obj_config2016_6565.mat');
circ = nstxu2016_circ(tok_data_struct);

eqs = load([cache_dir 'eqs' num2str(shot) '.mat']).eqs;

i = eqs.time(:)' <= 0 | [eqs.gdata(:).cpasma] < sqrt(eps);
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
close all

figure
k = 1;
hold on
grid on
% scatter(t,Idot(k,:))
plot(t,-Idot(k,:), '-b')
% plot(t,I(k,:), '-b')
xlim([0 0.1])

yyaxis right
k = 54;
hold on
% scatter(t,Idot(k,:))
plot(t,Idot(k,:), '-r')
% plot(t,I(k,:), '-r')
xlim([0 0.1])
%%

% plot(t,I(end,:))

% k = 8;
% figure
% hold on
% plot(t, Idot(k,:))
% scatter(t, smooth(Idot(k,:)))


struct_to_ws(tok_data_struct);
mpp = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit.build_inputs.tok_data_struct.mpp;


for i = 1:neq
  pcurrt = eqs.gdata(i).pcurrt(:);
  mcIp = mpc' * pcurrt / sum(pcurrt);
  mvIp = mpv' * pcurrt / sum(pcurrt);
  Lp(i) = pcurrt' * mpp * pcurrt / sum(pcurrt)^2;
  
  M = circ.Pxx' * [mcIp; mvIp; Lp(i)];
  
  Rp(i) = -1 / ip(i) * M' * Idot(:,i);
  
  A(i) = polyarea(eqs.gdata(i).rbbbs, eqs.gdata(i).zbbbs);
  
end

%%
figure
hold on
scatter(t,Rp)
plot(t,Rp)
% xlim([0 0.3])

yyaxis right
hold on
scatter(t,Lp)
plot(t,Lp)
% xlim([0 0.3])

figure
hold on
plot(t, Rp.*A)
scatter(t, Rp.*A)
% xlim([0 0.3])

%%
% compare to Spitzer resistivity



















%%
% %%
% % ================
% % Resistivity calc
% % ================
% 
% t0 = t_snapshots(1);
% tf = t_snapshots(end);
% N = 24;
% tspan = linspace(t0,tf,N+1);
% dt = mean(diff(tspan));
% 
% % load shot trajectory
% traj = make_traj(shot, t_snapshots, tspan);
% 
% 
% xdot = gradient(traj.x', dt)';
% ip = traj.ip;
% lstar_ip = squeeze(traj.lstar(:,end,:));
% 
% iuse = [1:49];
% for i = 1:length(ip)
%   Rp(i) = -1 / ip(i) * lstar_ip(i,iuse) * xdot(i,iuse)';
% end
% 
% % cftool(tspan,Rp)
% 
% ft = fittype( 'smoothingspline' );
% opts = fitoptions( 'Method', 'SmoothingSpline' );
% opts.Normalize = 'on';
% opts.SmoothingParam = 0.975524989202987;
% [fit_Rp, gof] = fit( tspan(:), Rp(:), ft, opts );
% 
% figure
% hold on
% plot(tspan, Rp)
% plot(fit_Rp)
% 
% if saveit
%   fn = '/Users/jwai/Research/rampup_nstxu/buildmodel/resistivity_calc/fit_Rp.mat';
%   save(fn, 'fit_Rp')
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % =========================================
% % Load eq snapshots to make shot trajectory
% % =========================================
% function traj = make_traj(shot, t_snapshots, tspan)
%   dt = mean(diff(tspan));
%   
%   for i = 1:length(t_snapshots)
%     
%     t = t_snapshots(i) * 1000;
%     load([num2str(shot) '_' num2str(t) '_sys.mat']);
%     lstar(i,:,:) = sys.lstar;
%     
%     eq = load(['eq' num2str(shot) '_' num2str(t) '.mat']);
%     eq = eq.eq;   
%     
%     % Load coil currents
%     load(['coils' num2str(shot) '.mat'])
%     
%     ccnames = {'OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', ...
%       'PF5', 'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL'};
%     
%     icoil = [];
%     for k = 1:length(sys.inputs)
%       icoil(k) = find(strcmp(ccnames, sys.inputs{k}));
%     end
%     
%     itime = find( floor(coils.t*1000) == t);
%     ic(i,:) = coils.ic(itime, icoil);
%     iv(i,:) = coils.iv(itime,:);
%     ip(i,:) = eq.cpasma;        
%   end  
%   
%   % finer interpolation
%   method = 'pchip';
%   traj.lstar  = interp1(t_snapshots, lstar,      tspan, 'previous');
%   traj.ic     = interp1(t_snapshots, ic,     tspan, method);
%   traj.iv     = interp1(t_snapshots, iv,     tspan, method);
%   traj.ip     = interp1(t_snapshots, ip,     tspan, method);
%   traj.x      = [traj.ic traj.iv traj.ip'];
% end








































































