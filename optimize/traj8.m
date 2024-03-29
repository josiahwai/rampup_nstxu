clear all; clc; close all
warning('off','MATLAB:polyshape:repairedBySimplify')
warning('off', 'curvefit:fit:nonDoubleYData')
warning('off', 'stats:pca:ColRankDefX')

ROOT = getenv('RAMPROOT');

enforce_stability = 0;

% shot = 204660;  % We will try to recreate this shot
% shot = 203708;
shot = 204069;

% t0 = 0.05;
t0 = 0.07;
tf = 0.9;
N = 50;
fetch_times = linspace(t0, tf, N);

tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';
opts.cache_dir = [ROOT '/fetch/cache/'];
efit01_eqs = fetch_eqs_nstxu(shot, fetch_times, tree, tokamak, server, opts);

% remove some zero-padding, oft causes trouble
for i = 1:N
  k = efit01_eqs.gdata(i).rbbbs==0 & efit01_eqs.gdata(i).zbbbs==0;
  efit01_eqs.gdata(i).rbbbs(k) = [];
  efit01_eqs.gdata(i).zbbbs(k) = [];  
end

init = efit01_eqs.gdata(1);

% load geometry
vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
tok_data_struct.imks = 1;
circ = nstxu2016_circ(tok_data_struct);
nr = tok_data_struct.nr;
nz = tok_data_struct.nz;
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;

% =================
% Optimizer targets
% =================
% targets are: desired boundary, boundary-defining pt, and Ip
targets.time = double(efit01_eqs.time);

% target boundary
gap_opts.plotit = 0;
for i = 1:N
  gaps(i) = get_nstxu_gaps(efit01_eqs.gdata(i), gap_opts);
end
targets.rcp = [gaps(:).r]';
targets.zcp = [gaps(:).z]';
ngaps = size(targets.rcp, 2);


% boundary defining point
bry_opts.plotit = 0;
for i = 1:N
  bry(i) = eq_bdef_analysis(efit01_eqs.gdata(i), tok_data_struct, bry_opts);
end
targets.rbdef = [bry(:).rbdef]';
targets.zbdef = [bry(:).zbdef]';
targets.islimited = [bry(:).islimited]';

% coil and vessel currents (in this case, solving for these so set to 0)
% targets.icx = zeros(N, circ.ncx);
% targets.ivx = zeros(N, circ.nvx);
targets.icx = ones(N,1) * init.icx';
targets.ivx = zeros(N, circ.nvx);


% Ip
targets.ip = [efit01_eqs.gdata(:).cpasma]';


% Wmhd
targets.wmhd = read_wmhd(efit01_eqs, tok_data_struct);

% li
for i = 1:N
  [~,~,~,targets.li(i,1)] = inductance(efit01_eqs.gdata(i), tok_data_struct);
end


% put onto equal-spaced timebase
t = linspace(targets.time(1), targets.time(end), length(targets.time));
ts = mean(diff(t));
N = length(t); 
fns = fieldnames(targets);
for i = 1:length(fns)
  targets.(fns{i})(1:N,:) = interp1(targets.time, targets.(fns{i}), t);
end
  

% rearrange data for easier access later
for i = 1:length(targets.time)  
  for j = 1:length(fns)
    targets_array(i).(fns{j}) = targets.(fns{j})(i,:);
  end
end


% =====================
% Optimizer constraints 
% =====================
% use nan to specify no constraint

constraints.icx = nan(N, length(init.icx)); 
constraints.ivx = nan(N, length(init.ivx));
constraints.ip = nan(N, 1);

constraints.icx(1:N, circ.iremove) = 0;   % these coils turned off
constraints.icx(1,:) = init.icx;

icx_experiment = [efit01_eqs.gdata(:).icx];

constraints.icx(1:N, 10) = icx_experiment(10,:); % PF2L
constraints.icx(1:N, 5) = icx_experiment(5,:);   % PF2U
% constraints.icx(t<0.4, [5 10]) = 0; % PF2U/L constrained to 0 for first part of shot
% constraints.icx(1:N,1) = icx_experiment(1,:); % OH

% =================
% Optimizer weights
% =================
wt.icx = ones(N,circ.ncx) * 1e-6;  % wt.icx(1:N, [2 13]) = 1e-4;
wt.ivx = ones(N,circ.nvx) * 1e-6;
wt.ip = ones(N,circ.np) * 3e-5;
wt.cp = ones(N, ngaps) * 2e8;
wt.bdef = double(~targets.islimited) * 0; 

wt.dicxdt = ones(size(wt.icx)) / ts^2 * 1e-5;
wt.divxdt = ones(size(wt.ivx)) / ts^2 * 0;
wt.dipdt = ones(size(wt.ip))  / ts^2 * 0;
wt.dcpdt = ones(size(wt.cp)) / ts^2 * 0;
wt.dbdefdt = ones(size(wt.bdef)) / ts^2 * 0;

% ====================================
% Estimate plasma current distribution
% ====================================

% for i = 1:N  
%   target = targets_array(i);
%   eq = efit01_eqs.gdata(i);
%   
%   pla_opts.plotit = 1;
%   pla_opts.cold_start = 0;    
%   pla_opts.ref_eq = eq;  
%   pla_opts.debug_use_actual = 0; % t(i) < 0.4;
%   pla_opts.time = t(i);
%   
%   pla(i) = estimate_pla(target, tok_data_struct, eq, pla_opts);
% end


for i = 1:N  
  i
  target = targets_array(i); 
  opts.init = efit01_eqs.gdata(i);
  opts.plotit = 0;
  opts.max_iterations = 3;
  pla(i) = semifreegs(target, tok_data_struct, opts);  
end

% load('pla.mat')

% for i = 1:N
%   eq = efit01_eqs.gdata(i);
%   eq.area = polyarea(eq.rbbbs, eq.zbbbs);
%   eq.psizr_app = reshape(tok_data_struct.mpc*eq.ic + tok_data_struct.mpv*eq.iv, nr, nz);
%   eq.psizr_pla = eq.psizr-eq.psizr_app;
%   pla(i) = eq;
% end



% load('betap')
% for i = 1:N
%   i
%   
%   if i == 1
%     init =  efit01_eqs.gdata(10);
%     opts.max_iterations = 15;
%   else
%     init = pla(i-1);
%     opts.max_iterations = 9;
%   end
%   
%   target = targets_array(i);
%   target.betap = betap(i);
%   pla(i) = estimate_pla_gsdesign(target, tok_data_struct, init, opts);    
%   
%   Lp(i) = inductance(pla(i), tok_data_struct);
%   Lp2(i) = inductance(efit01_eqs.gdata(i), tok_data_struct);
% end
% 
% for i = 1:N
%   pla(i).area = polyarea(pla(i).rbbbs, pla(i).zbbbs);
% end
%
% load('gsdesign_pla.mat')






% ==========================
% Estimate plasma parameters
% ==========================

% Load fitted vacuum model parameters
mxx = vac_sys.sysid_fits.Mxx;
rxx = vac_sys.sysid_fits.Rxx;
mvv = mxx(circ.iivx, circ.iivx);
mcc = mxx(circ.iicx, circ.iicx);
mvc = mxx(circ.iivx, circ.iicx);
Rv = rxx(circ.iivx);
Rc = rxx(circ.iicx);
mpc = tok_data_struct.mpc * circ.Pcc;
mpv = tok_data_struct.mpv * circ.Pvv;
mpp = tok_data_struct.mpp;

% Estimate time-dependent parameters
[eta, eta_t] = load_eta_profile();
params.eta = interp1(eta_t, eta, t)';
for i = 1:N
  pcurrt = pla(i).pcurrt(:);
  ip = sum(pcurrt(:));
  params.mcIp(i,:) = mpc' * pcurrt / ip;
  params.mvIp(i,:) = mpv' * pcurrt / ip;
  params.Lp(i,:) = pcurrt' * mpp * pcurrt / ip^2;
  params.Rp(i,:) = params.eta(i) / pla(i).area;
end
% res1 = load('/Users/jwai/Research/rampup_nstxu/sysid/plasma_resistance/fits/res204660.mat').res;
% res2 = load('/Users/jwai/Research/rampup_nstxu/buildmodel/resistivity_calc/res/res204660.mat').res;
% [rp, rp_t] = load_rp_profile();
% plot(res1.t, res1.Rp, res2.t, res2.Rp, rp_t, rp,'linewidth',3);
% legend('Res1','Res2','rp_profile','fontsize',12)

% res = load('/Users/jwai/Research/rampup_nstxu/buildmodel/resistivity_calc/res/res204660.mat').res;
% params.Rp = double(interp1(res.t, res.Rp2, t)');

[rp, rp_t] = load_rp_profile();
params.Rp = interp1(rp_t, rp, t)';

% load('res203708.mat')
% params.Rp = double(interp1(res.t, res.Rp, t)');

% for i = 1:N
%   [Lp(i), Li(i), Le(i), li(i)] = inductance(efit01_eqs.gdata(i), tok_data_struct);
% end



% Form the time-dependent A,B,C,D matrices
% cdata = build_cmat_data(eq0, circ, tok_data_struct, targ_geo);
% DC = [cdata.x; cdata.dpsicpdix];

for i = 1:N
  % dynamics A, B matrices
  M = [mvv params.mvIp(i,:)'; params.mvIp(i,:) params.Lp(i)];
  R = diag([Rv; params.Rp(i)]);
  Minv = inv(M);
  Ac = -Minv*R;
  Bc = -Minv * [mvc; params.mcIp(i,:)];
  if enforce_stability
    Ac = numerically_stabilize(Ac, 1e3);
  end
  [A{i}, B{i}] = c2d(Ac, Bc, ts);
  
  % output C, D matrices
  r = vacuum_response(pla(i), targets_array(i), tok_data_struct);
  response = [r.disdis; r.dpsicpdis - r.dpsibrydis; r.dpsibrydis_r; r.dpsibrydis_z];
  % response = [r.disdis; r.dpsicpdis - r.dpsibrydis];
  C{i} = response(:, [circ.iivx circ.iipx]);
  D{i} = response(:, circ.iicx);
  
  %   D{i} = DC(:, 1:circ.ncx);
  %   C{i} = DC(:, circ.ncx+1:end);
end


% ================================
% Measure outputs, first iteration
% ================================
% y := [icx ivx ip (psicp-psibry) dpsibrydr dpsibrydz]'
clear ic0hat x0hat y0hat
for i = 1:N
  psizr_app = reshape(mpc*init.icx + mpv*init.ivx, nr, nz);
  psizr_pla = pla(i).psizr_pla;
  psizr = psizr_pla + psizr_app;
  currents = [init.icx; init.ivx; init.cpasma];
  
  ic0hat(:,i) = init.icx;
  x0hat(:,i) = [init.ivx; init.cpasma];
  y0hat(:,i) = measure_y(psizr, currents, targets_array(i), tok_data_struct);
  % y0hat(:,i) = [init.icx; init.ivx; init.cpasma; cp_diff0];
  
end
ic0 = ic0hat(:,1);
x0 = x0hat(:,1);
y0 = y0hat(:,1);

ic0hat = ic0hat(:);
x0hat = x0hat(:);
y0hat = y0hat(:);

xk = x0;
ny = length(y0);

xprev = [x0; zeros((N-1)*(circ.nvx+1),1)];
icprev = [ic0; zeros((N-1)*circ.ncx,1)];
yprev = [y0; zeros((N-1)*ny, 1)];



% ==========================
% form the prediction model
% ==========================
npv = circ.nvx + circ.np;
Apow  = eye(npv);
F  = [];
F_row = zeros(npv, N*circ.ncx);
E  = [];

for i = 1:N
  idx = (circ.ncx*(i-1)+1):(circ.ncx*i);
  F_row = A{i} * F_row;
  F_row(:,idx) = B{i};
  F = [F; F_row];
  Apow = A{i} * Apow;
  E = [E; Apow];
end
F = F / ts;

Chat = blkdiag(C{:});
Dhat = blkdiag(D{:});

% ============
% form costfun
% ============
% y := [icx ivx ip (psicp-psibry) dpsibrydr dpsibrydz]'

% velocity conversion matrices
Svp = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.nvx+circ.np));
Sc = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.ncx));
Sy = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(ny));

% weights and targets
for i = 1:N
  Q{i} = diag([wt.icx(i,:) wt.ivx(i,:) wt.ip(i) wt.cp(i,:) wt.bdef(i,:) wt.bdef(i,:)]);
  Qv{i} = diag([wt.dicxdt(i,:) wt.divxdt(i,:) wt.dipdt(i) wt.dcpdt(i,:) wt.dbdefdt(i,:) wt.dbdefdt(i,:)]);
  
  %   Q{i} = diag([wt.icx(i,:) wt.ivx(i,:) wt.ip(i,:) wt.cp(i,:)]);
  %   Qv{i} = diag([wt.dicxdt(i,:) wt.divxdt(i,:) wt.dipdt(i,:) wt.dcpdt(i,:)]);
  
end
Qhat = blkdiag(Q{:});
Qvhat = blkdiag(Qv{:});
Qbar = Qhat + Sy'*Qvhat*Sy;


%%
for iteration = 1:1

  % update measurements of y
  if iteration > 1
    x0hat = prevsoln.xhat * 0;
    ic0hat = prevsoln.icxhat * 0;
    y0hat = prevsoln.yhat * 0;
    
    for i = 1:N
      icx = prevsoln.icxhat(:,i);
      ivx = prevsoln.ivxhat(:,i);
      ip = prevsoln.iphat(i);
      psizr_app = reshape(mpc*init.icx + mpv*init.ivx, nr, nz);
      psizr_pla = pla(i).psizr_pla;
      psizr = psizr_pla + psizr_app;
      currents = [icx; ivx; ip];
      ic0hat(:,i) = icx;
      x0hat(:,i) = [ivx; ip];
      y0hat(:,i) = measure_y(psizr, currents, targets_array(i), tok_data_struct);  
    end
  end
  
  % targets for y in vector form
  rhat = [targets.icx targets.ivx targets.ip zeros(size(targets.rcp)) zeros(N,2)]';
  % rhat = [targets.icx targets.ivx targets.ip zeros(size(targets.rcp))]';
  % rhat = [targets.icx targets.ivx targets.ip targets.cp_diff]';
  rhat = rhat(:); 

  % cost function
  z = y0hat(:) + Chat * (E*xk - F*icprev - x0hat(:)) - Dhat*ic0hat(:);
  M = Chat*F*Sc + Dhat;
  H = M'*Qbar*M;
  ft = (z'*Qbar - rhat'*Qhat - yprev'*Qvhat*Sy) * M;

  npv = size(H,1);


  % ===============================
  % Constraints: Aeq * ichat = beq
  % ===============================
  clear Aeq beq
  ieq = 1;

  % Equality constraints
  % --------------------

  % coil currents at the specified times
  Aeq{ieq} = eye(N * circ.ncx);
  beq{ieq} = reshape(constraints.icx', [], 1);

  i = isnan(beq{ieq});
  beq{ieq}(i,:) = [];
  Aeq{ieq}(i,:) = [];

  Aeq = cat(1,Aeq{:});
  beq = cat(1,beq{:});

  % % enforce linear independence
  % [u,s,v] = svd_tol(Aeq, 0.9999);
  % Aeq = u'*Aeq;
  % beq = u'*beq;


  % Inequality constraints
  % ----------------------
  Aineq = zeros(0,npv);
  bineq = zeros(0,1);


  % ==============
  % Solve quadprog
  % ==============

  % Solve quadratic program
  opts = mpcActiveSetOptions;
  iA0 = false(length(bineq), 1);

  [icxhat,exitflag,iA,lambda] = mpcActiveSetSolver(H, ft', Aineq, bineq, Aeq, beq, iA0, opts);

  % unpack solution
  yhat = M*icxhat + z;

  yhat = reshape(yhat,[],N);
  icxhat = reshape(icxhat,[],N);
  ivxhat = yhat(circ.iivx,:);
  iphat = yhat(circ.iipx,:);
  xhat = yhat([circ.iivx circ.iipx],:);
 
  prevsoln = variables2struct(xhat, icxhat, ivxhat, iphat, yhat);

  
%   
%   % ====================================
%   % ANALYZE EQULIBRIUM
%   % ====================================
%   for i = 1:N  
%     
%     % analyze
%     icx = icxhat(:,i);
%     ivx = ivxhat(:,i);
%     ip = iphat(i);      
%     psizr_app = reshape(mpc*icx + mpv*ivx, nr, nz);  
%     psizr_pla = pla(i).psizr_pla;
%     psizr = psizr_app + psizr_pla;
% 
%     eq_opts.plotit = 0;
%     eq_opts.robust_tracing = 1;
%     
%     eq = eq_analysis(psizr, pla(i), tok_data_struct, eq_opts);    
%     eq = append2struct(eq, icx, ivx, ip, psizr_app, psizr_pla);
%     eqs(i) = eq;
%     
%     % scatter(targets.rcp(i,:), targets.zcp(i,:), 'k', 'filled')
%   end
%   

  % ====================================
  % Estimate plasma current distribution
  % ====================================
%   for i = 1:N  
%     target = targets_array(i);
%     eq = eqs(i);
%     
%     pla_opts.plotit = 0;
%     pla_opts.cold_start = 0;    
%     pla_opts.ref_eq = efit01_eqs.gdata(i);
%     
%     pla(i) = estimate_pla(target, tok_data_struct, eq, pla_opts);
%   end


end



close all
%%
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% DEBUGGING
% DEBUGGING
% DEBUGGING
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% 
% i = 40;
% 
% icx = icxhat(:,i);
% ivx = ivxhat(:,i);
% ip = iphat(i);      
% psizr_app = reshape(mpc*icx + mpv*ivx, nr, nz);  
% psizr_pla = pla(i).psizr_pla;
% psizr = psizr_app + psizr_pla;
% 
% % eq_opts.plotit = 0;
% % eq_opts.robust_tracing = 1;
% % eq = eq_analysis(psizr, pla(i), tok_data_struct, eq_opts);    
% % eq = append2struct(eq, icx, ivx, ip, psizr_app, psizr_pla);
% % eqs2(i) = eq;
% 
% 
% figure
% hold on
% plot_eq(efit01_eqs.gdata(i));
% contour(rg, zg, eqs(i).psizr, [eqs(i).psibry, eqs(i).psibry], '--b')
% % contour(rg, zg, eqs2(i).psizr, [eqs2(i).psibry, eqs2(i).psibry], '--b')
% target = targets_array(i);
% scatter(target.rcp, target.zcp)

%%

figure
hold on
plot(t, icxhat(2:end,:)/1e3, '-b')
icx_efit = [efit01_eqs.gdata(:).icx];
plot(t, icx_efit(2:end,:)/1e3 , '--r')
for i = circ.iicx_keep(2:end)  
  text(t(end)+.002, icxhat(i,end)/1e3, circ.ccnames{i}, 'fontsize', 14)
end
xlabel('Time [s]', 'fontsize', 16)
ylabel('Coil currents [kA]', 'fontsize', 14)
l = mylegend({'Optimizer', 'Experiment'}, {'-', '--'}, [], {'b', 'r'});
l.FontSize = 14;
title('Coil trajectories', 'fontsize', 18)
set(gcf, 'Position', [669 326 505 314])
xlim([0 1.02])


figure
hold on
plot(t,iphat/1e6,'--b', t,targets.ip/1e6,'-b')
ylabel('Ip [MA]', 'fontsize', 18)
yyaxis right
hold on
plot(t,icxhat(1,:)/1e3, 'r')
plot(t, icx_efit(1,:)/1e3, '--r')
ylabel('OH [kA]', 'fontsize', 18)
title('Ip trajectory', 'fontweight', 'bold', 'fontsize', 16)
xlabel('Time [s]', 'fontsize', 16)
l = mylegend({'Optimizer', 'Experiment'}, {'-', '--'}, [], {[1 1 1]*.02, [1 1 1]*0.2});
l.FontSize = 16;
l.Location = 'best';
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = 'r';
set(gcf, 'Position', [680 723 487 255])



%%


% DEBUGGING: gsdesign comparison
i = 50;

icx = icxhat(:,i);
ivx = ivxhat(:,i);
ip = iphat(i);
% icx = efit01_eqs.gdata(i).icx;
% ivx = efit01_eqs.gdata(i).ivx;
% ip = efit01_eqs.gdata(i).cpasma;

init = efit01_eqs.gdata(i);
opts.pres = efit01_eqs.gdata(i).pres;
opts.fpol = efit01_eqs.gdata(i).fpol;

x = [icx; ivx; ip];
[spec, init, config] = make_gsdesign_inputs2(x, tok_data_struct, init, circ, opts);
eq = gsdesign(spec,init,config);

figure
hold on
eqt = efit01_eqs.gdata(i);
contour(rg, zg, eqt.psizr, [eqt.psibry eqt.psibry], '--b')
contour(rg, zg, eq.psizr, [eq.psibry eq.psibry], '-r')
plot_eq(eq)
contour(rg, zg, eqt.psizr, [eqt.psibry eqt.psibry], '--b')
% contour(rg, zg, init.psizr, [init.psibry init.psibry], '--b')
scatter(targets.rcp(i,:), targets.zcp(i,:), 'k', 'filled')
l = legend('EFIT01', 'optimizer', 'fontsize', 14);
set(gcf, 'Position', [634 449 367 529])
%%
% 
% % DEBUGGING: gsdesign comparison
% 
% for i = 1:5:N  
%   close all
%   target = targets_array(i);
%   eq = efit01_eqs.gdata(i);
%   
%   pla_opts.plotit = 1;
%   pla_opts.cold_start = 0;    
%   pla_opts.ref_eq = eq;
%   pla_opts.debug_use_actual = 0;
%   pla(i) = estimate_pla(target, tok_data_struct, eq, pla_opts);
%   
%   title(t(i))
% end



%%
% %%
% i = 40;
% 
% load('gsdesign_pla.mat')
% 
% eq = efit01_eqs.gdata(i);
% psizr_app = reshape(mpc*eq.icx + mpv*eq.ivx, nr, nz);
% 
% figure
% [~,cs] = contour(eq.rg, eq.zg, eq.psizr-psizr_app, 20, '--r');
% hold on
% contour(eq.rg, eq.zg, pla(i).psizr_pla, cs.LevelList, 'b')
% % contour(eq.rg, eq.zg, x.psizr_pla, cs.LevelList, 'b')
% 
% 
% target = targets_array(i);
% scatter(target.rcp, target.zcp, 'k', 'filled')


%%
% figure
% plot(t,[efit01_eqs.gdata(:).ivx])
% 
% figure
% plot(t,ivxhat)

























