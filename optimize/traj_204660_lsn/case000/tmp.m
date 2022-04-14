% TO DO:
% 
% iterate to get F5/F8 equilibrium (ie can I get w/o first using
% gsdesign?)
%
% make it run fast and on iris
% 
% load inputs from previous shot
%  - equilibria
%  - resistance calc

% clear all; clc; close all
ROOT = getenv('RAMPROOT');

clearvars -except targets targets_array wt constraints pla pla_array ROOT shot init eqs Rp

%%
t = targets.time + 1.5;
N = length(t);
ts = mean(diff(t));

tok_data_struct = load('d3d_obj_mks_struct_6565.mat').tok_data_struct;
load('mpp_full');
tok_data_struct.mpp_full = mpp_full;
struct_to_ws(tok_data_struct);
circ = d3d_circ(tok_data_struct);
iy = circ.iy;

% ==========================
% Estimate plasma parameters
% ==========================

% Load fitted vacuum model parameters
mxx = [mcc mcv; mvv mcv'];
rxx = [resc; resv];
mvc = mcv';
Rv = resv;
Rc = resc;


% Estimate time-dependent parameters
for i = 1:N
  pcurrt = pla_array(i).pcurrt(:);
  ip = sum(pcurrt(:));
  params.mcIp(i,:) = mpc' * pcurrt / ip;
  params.mvIp(i,:) = mpv' * pcurrt / ip;
  params.Lp(i,:) = pcurrt' * mpp_full * pcurrt / ip^2;
end

% [rp, rp_t] = load_rp_profile2(0);
% params.Rp = interp1(rp_t, rp, t, 'linear', 'extrap')';

params.Rp(1:N,:) = Rp; 


% Form the time-dependent A,B,C,D matrices
for i = 1:N
  % dynamics A, B matrices
  M = [mvv params.mvIp(i,:)'; params.mvIp(i,:) params.Lp(i)];
  R = diag([Rv; params.Rp(i)]);
  Minv = inv(M);
  Ac = -Minv*R;
  Bc = -Minv * [mvc; params.mcIp(i,:)];
  % Ac = numerically_stabilize(Ac, 1e3);
  [A{i}, B{i}] = c2d(Ac, Bc, ts);
  
  % output C, D matrices
  r = vacuum_response(pla_array(i), targets_array(i), tok_data_struct);
  response = [r.disdis; r.dpsicpdis - r.dpsibrydis; r.dpsibrydis_r; ...
    r.dpsibrydis_z; r.dpsix2dis_r; r.dpsix2dis_z];
 
  C{i} = response(:, [circ.iiv circ.iip]);
  D{i} = response(:, circ.iic);
end


% ================================
% Measure outputs, first iteration
% ================================
% y := [icx ivx ip (psicp-psibry) dpsibrydr dpsibrydz]'
clear ic0hat x0hat y0hat
for i = 1:N
  psizr_app = reshape(mpc*init.ic + mpv*init.iv, nr, nz);
  psizr_pla = pla_array(i).psizr_pla;
  psizr = psizr_pla + psizr_app;
  currents = [init.ic; init.iv; init.cpasma];
  
  ic0hat(:,i) = init.ic;
  x0hat(:,i) = [init.iv; init.cpasma];
 
  y = measure_y0(psizr, currents, targets_array(i), tok_data_struct);
  
  y0hat(:,i) = [y.currents; y.psicp(:) - y.psibry; y.psibry_r; ...
    y.psibry_z; y.psix2_r; y.psix2_z];
end
ic0 = ic0hat(:,1);
x0 = x0hat(:,1);
y0 = y0hat(:,1);

ic0hat = ic0hat(:);
x0hat = x0hat(:);
y0hat = y0hat(:);

xk = x0;
ny = length(y0);

xprev = [x0; zeros((N-1)*(circ.nv+1),1)];
icprev = [ic0; zeros((N-1)*circ.nc,1)];
yprev = [y0; zeros((N-1)*ny, 1)];



% ==========================
% form the prediction model
% ==========================
npv = circ.nv + circ.np;
Apow  = eye(npv);
F  = [];
F_row = zeros(npv, N*circ.nc);
E  = [];

for i = 1:N
  idx = (circ.nc*(i-1)+1):(circ.nc*i);
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
Svp = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.nv+circ.np));
Sc = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.nc));
Sy = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(ny));

% curvature conversion: d2y = S2y * y
m = tridiag(1, -2, 1, N);
m(N,:) = m(N-1,:);
m(1,:) = m(2,:);
S2y = kron(m, eye(ny));


% weights and targets, Q = weights on values, Qv = weights on derivatives
% (velocities), Q2v = weights on 2nd derivatives 
for i = 1:N
  
  Q{i} = diag([wt.icx(i,:) wt.ivx(i,:) wt.ip(i) wt.cp(i,:) ...
    wt.bdef(i) wt.bdef(i) wt.xp2(i) wt.xp2(i)]);
  
  Qv{i} = diag([wt.dicxdt(i,:) wt.divxdt(i,:) wt.dipdt(i) wt.dcpdt(i,:) ...
    wt.dbdefdt(i) wt.dbdefdt(i)  wt.dxp2dt(i) wt.dxp2dt(i)]); 
  
  Q2v{i} = diag([wt.d2icxdt2(i,:) wt.d2ivxdt2(i,:) wt.d2ipdt2(i) wt.d2cpdt2(i,:) ...
    wt.d2bdefdt2(i) wt.d2bdefdt2(i) wt.d2xp2dt2(i) wt.d2xp2dt2(i)]);  
end

Qhat = blkdiag(Q{:});
Qvhat = blkdiag(Qv{:});
Q2vhat = blkdiag(Q2v{:});
Qbar = Qhat + Sy'*Qvhat*Sy + S2y'*Q2vhat*S2y;


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
      psizr_app = reshape(mpc*init.ic + mpv*init.ivx, nr, nz);
      psizr_pla = pla_array(i).psizr_pla;
      psizr = psizr_pla + psizr_app;
      currents = [icx; ivx; ip];
      ic0hat(:,i) = icx;
      x0hat(:,i) = [ivx; ip];
      y0hat(:,i) = measure_y(psizr, currents, targets_array(i), tok_data_struct);  
    end
  end
  
  % targets for y in vector form
  % rhat = [targets.icx(2:end,:) targets.ivx(2:end,:) targets.ip(2:end) zeros(size(targets.rcp(2:end,:))) zeros(N,2)]';
  rhat = [targets.icx targets.ivx targets.ip zeros(size(targets.rcp)) zeros(N,2) zeros(N,2)]';
  % rhat = rhat(:,2:end);
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
  Aeq{ieq} = eye(N * circ.nc);
  beq{ieq} = reshape(constraints.icx', [], 1);
  i = isnan(beq{ieq});
  beq{ieq}(i,:) = [];
  Aeq{ieq}(i,:) = [];
    
  % constraints of the form A*ic=b, applied to every timestep
  ieq = ieq+1;
  Aeq{ieq} = kron(eye(N), constraints.A1);
  beq{ieq} = ones(N,1) * constraints.b1;
        
  ieq = ieq+1;
  Aeq{ieq} = kron(eye(N), constraints.A2);
  beq{ieq} = ones(N,1) * constraints.b2;
  
  
  Aeq = cat(1,Aeq{:});
  beq = cat(1,beq{:});

  % % enforce linear independence
  [u,s,v] = svd_tol(Aeq, 0.9999);
  Aeq = u'*Aeq;
  beq = u'*beq;

  %   Aeq = zeros(0,npv);
  %   beq = zeros(0,1);
  
  % Inequality constraints
  % ----------------------
  Aineq = kron(eye(N), constraints.Aineq1);
  bineq = repmat(constraints.bineq1, N, 1);
  
  i = isnan(bineq);
  bineq(i) = [];
  Aineq(i,:) = [];
  
%   Aineq = zeros(0,npv);
%   bineq = zeros(0,1);


  % ==============
  % Solve quadprog
  % ==============

  % Solve quadratic program
%   opts = mpcActiveSetOptions;
%   iA0 = false(length(bineq), 1);
% 
%   [icxhat,exitflag,iA,lambda] = mpcActiveSetSolver(H, ft', Aineq, bineq, Aeq, beq, iA0, opts);

  icxhat = quadprog((H+H')/2, ft', Aineq, bineq, Aeq, beq);

  % unpack solution
  yhat = M*icxhat + z;

  yhat = reshape(yhat,[],N);
  icxhat = reshape(icxhat,[],N);
  ivxhat = yhat(circ.iiv,:);
  iphat = yhat(circ.iip,:);
  xhat = yhat([circ.iiv circ.iip],:);
 
  prevsoln = variables2struct(xhat, icxhat, ivxhat, iphat, yhat);

end

close all

%%
ic_true = targets.icx_true';
ic_pred = icxhat;
icoil = 3:20;

icoil = max(abs(ic_pred' - ic_true')) > 200;
icoil(1:2) = 0;
icoil([iy.F3A iy.F4A iy.F5A iy.F8A iy.F9A]) = 1;
icoil = find(icoil);


figure
hold on
plot(t, ic_pred(icoil,:), '-b')
plot(t, ic_true(icoil,:), '--r')
for i = icoil
  text(t(end)+.002, ic_pred(i,end), tok_data_struct.ccnames(i,:), 'fontsize', 14)
end
xlabel('Time [s]', 'fontsize', 16)
ylabel('Coil currents [A]', 'fontsize', 14)
l = mylegend({'Optimizer', 'Experiment'}, {'-', '--'}, [], {'b', 'r'});
l.FontSize = 14;
title('Coil trajectories', 'fontsize', 18)




%%
ip_pred = iphat;

figure
subplot(211)
hold on
plot(t,ip_pred/1e6,'-b', t,targets.ip/1e6,'--r')
ylabel('Ip [MA]', 'fontsize', 18)

subplot(212)
hold on
plot(t,ic_pred(1:2,:)/1e3, '-b')
plot(t, ic_true(1:2,:)/1e3, '--r')
ylabel('OH [kA]', 'fontsize', 18)
xlabel('Time [s]', 'fontsize', 16)




%%
% DEBUGGING: gsdesign comparison
if 0
  i = 5;
  
  %   icx = icxhat(:,i);
  
  icx = init.ic;
  icoil = [iy.F3A iy.F4A iy.F4B iy.F5A iy.F8A iy.F9A];
  % icoil = [iy.F4A iy.F5A];
  icx(icoil) = icxhat(icoil,i);
  
  ivx = ivxhat(:,i);
  ip = iphat(i);
  
  init = eqs{2};
  opts.pres = init.pres;
  opts.fpol = init.fpol;
  
  x = [icx; ivx; ip];
  [spec, init, config] = make_gsdesign_inputs2(x, tok_data_struct, init, circ, opts);
  eq = gsdesign(spec,init,config);
  
  figure
  hold on
  contour(rg, zg, init.psizr, [init.psibry init.psibry], '--b')
  contour(rg, zg, eq.psizr, [eq.psibry eq.psibry], '-r')
  plot_eq(eq)
  contour(rg, zg, init.psizr, [init.psibry init.psibry], '--b')
  scatter(targets.rcp(i,:), targets.zcp(i,:), 'k', 'filled')
  
  
  [rx2, zx2] = isoflux_xpFinder(init.psizr, 1.15, 1.15, rg, zg);
  scatter(rx2, zx2, 60, 'b', 'filled')
  scatter(eq.nulls.r, eq.nulls.z, 60, 'r', 'filled')
  % scatter(targets.rx2(i), targets.zx2(i), '*', 'b')
  
  l = legend('Init', 'Final', 'fontsize', 14);
  set(gcf, 'Position', [634 449 367 529])
end


%% write to file

t = t + 5-t(1);

ccnames = cellstr(tok_data_struct.ccnames);
coils.t = t;
for i = 1:nc
  coils.(ccnames{i}) = ic_pred(i,:);
end

if 0
  ic2file(coils)
end


%%

psizr_app = reshape(mpc*icxhat + mpv*ivxhat, nz, nr, []);
psizr_pla = pla_array(1).psizr_pla;
psizr = psizr_app + psizr_pla;

for i = 1:N
  psi = squeeze(psizr(:,:,i));
  [rx2(i), zx2(i)] = isoflux_xpFinder(psi, 1.18, 1.15, rg, zg);
end

figure
plot_eq(init)
scatter(rx2, zx2, 'filled')

% figure
% subplot(211)
% plot(t, rx2)
% subplot(212)
% plot(t, zx2)


































