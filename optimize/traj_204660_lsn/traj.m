
targets = struct_fields_to_double(targets);

opts.cache_dir = [ROOT '/fetch/cache/'];
efit01_eqs = fetch_eqs_nstxu(shot, t, 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts);


% rearrange data for easier access later
fns = fieldnames(targets);
for i = 1:length(targets.time)  
  for j = 1:length(fns)
    targets_array(i).(fns{j}) = targets.(fns{j})(i,:);
  end
end

init = efit01_eqs.gdata(1);
ngaps = size(targets.rcp,2);
ts = mean(diff(t));

params.Rp = interp1(rp_sig.times, rp_sig.sigs, t);
targets.icx = zeros(N, circ.ncx);
targets.ivx = zeros(N, circ.nvx);

% ==================================
% Optimizer constraints 
% (use nan to specify no constraint)
% ==================================

% Equality constraints 
% ----------------------
constraints.icx = nan(N, length(init.icx)); 
constraints.ivx = nan(N, length(init.ivx));
constraints.ip = nan(N, 1);

constraints.icx(1:N, circ.iremove) = 0;   % these coils turned off
constraints.icx(1,:) = init.icx;

% icx_experiment = [efit01_eqs.gdata(:).icx];
% constraints.icx(1:N, 10) = icx_experiment(10,:); % PF2L
% constraints.icx(1:N, 5) = icx_experiment(5,:);   % PF2U
constraints.icx(t<0.4, [5 10]) = 0;   % PF2U/L constrained to 0 for first part of shot


% Inequality constraints 
% ----------------------
constraints_min.icx = nan(N, circ.ncx);  
constraints_min.icx(:, circ.ii_unipolar) = 0;
  

% =================
% Optimizer weights
% =================
wt.icx = ones(N,circ.ncx) * 1e-5; 
% wt.icx(iexclude,:) = 1e-10;
wt.icx(:,:) = 0;

wt.ivx = ones(N,circ.nvx) * 1e-6;
% wt.ivx(iexclude,:) = 1e-7;
wt.ivx(:,:) = 0;

wt.ip = ones(N,circ.np) * 3e-5;
% wt.ip(iexclude,:) = 3e-5;

wt.cp = ones(N, ngaps) * 2e8;
% wt.cp(iexclude,:) = 2e6;
% wt.cp(:,:) = 2e6;

wt.bdef = double(~targets.islimited) * 2e8; 
% wt.bdef(iexclude,:) = wt.bdef(iexclude,:) * 1e-2;
% wt.bdef(:,:) = 0;


% first derivative (velocity) weights
wt.dicxdt = ones(size(wt.icx)) / ts^2 * 1e-7;  % wt.dicxdt(:,1) = 1e-8;
wt.divxdt = ones(size(wt.ivx)) / ts^2 * 0;
wt.dipdt = ones(size(wt.ip))  / ts^2 * 0;
wt.dcpdt = ones(size(wt.cp)) / ts^2 * 0;
wt.dbdefdt = ones(size(wt.bdef)) / ts^2 * 0;


% second derivative (curvature) weights
wt.d2icxdt2 = ones(size(wt.icx))   / ts^4 * 2e-10;  
wt.d2ivxdt2 = ones(size(wt.ivx))   / ts^4 * 0;
wt.d2ipdt2 = ones(size(wt.ip))     / ts^4 * 0;
wt.d2cpdt2 = ones(size(wt.cp))     / ts^4 * 0;
wt.d2bdefdt2 = ones(size(wt.bdef)) / ts^4 * 0;










% ====================================
% Estimate plasma current distribution
% ====================================

if 0
  for i = 1:N
    i
    close all
    target = targets_array(i);
    opts.init = efit01_eqs.gdata(i);
    opts.plotit = 0;
    opts.max_iterations = 10;
    pla(i) = semifreegs(target, tok_data_struct, opts);

    pla(i).li
    target.li
    % plot( efit01_eqs.gdata(i).rbbbs, efit01_eqs.gdata(i).zbbbs, 'g')
  end
else
  load('pla204660.mat')
end


if 1
  psizr_pla = [pla(:).psizr_pla];
  psizr_pla = reshape(psizr_pla, nz, nr, []);
  x = timeseries(psizr_pla, t);
  y = smooth_ts(x);
  psizr_pla = y.Data;
  for i = 1:N
    pla(i).psizr_pla = squeeze(psizr_pla(:,:,i));
  end
end




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
% [eta, eta_t] = load_eta_profile();
% params.eta = interp1(eta_t, eta, t)';
for i = 1:N
  pcurrt = pla(i).pcurrt(:);
  ip = sum(pcurrt(:));
  params.mcIp(i,:) = mpc' * pcurrt / ip;
  params.mvIp(i,:) = mpv' * pcurrt / ip;
  params.Lp(i,:) = pcurrt' * mpp * pcurrt / ip^2;  
end




% Form the time-dependent A,B,C,D matrices

for i = 1:N
  % dynamics A, B matrices
  M = [mvv params.mvIp(i,:)'; params.mvIp(i,:) params.Lp(i)];
  R = diag([Rv; params.Rp(i)]);
  Minv = inv(M);
  Ac = -Minv*R;
  Bc = -Minv * [mvc; params.mcIp(i,:)];
  %   if enforce_stability
  %     Ac = numerically_stabilize(Ac, 1e3);
  %   end
  [A{i}, B{i}] = c2d(Ac, Bc, ts);
  
  % output C, D matrices
  r = vacuum_response(pla(i), targets_array(i), tok_data_struct);
  response = [r.disdis; r.dpsicpdis - r.dpsibrydis; r.dpsibrydis_r; r.dpsibrydis_z];

  C{i} = response(:, [circ.iivx circ.iipx]);
  D{i} = response(:, circ.iicx);   
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

% velocity conversion matrices: dy = Sy*y - yprev
Svp = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.nvx+circ.np));
Sc = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(circ.ncx));
Sy = kron(diag(ones(N,1)) + diag(-1*ones(N-1,1), -1), eye(ny));

% curvature conversion: d2y = S2y * y
m = tridiag(1, -2, 1, N);
m(N,:) = m(N-1,:);
m(1,:) = m(2,:);
S2y = kron(m, eye(ny));


% weights and targets
for i = 1:N
  Q{i} = diag([wt.icx(i,:) wt.ivx(i,:) wt.ip(i) wt.cp(i,:) wt.bdef(i,:) wt.bdef(i,:)]);
  Qv{i} = diag([wt.dicxdt(i,:) wt.divxdt(i,:) wt.dipdt(i) wt.dcpdt(i,:) wt.dbdefdt(i,:) wt.dbdefdt(i,:)]);
  Q2v{i} = diag([wt.d2icxdt2(i,:) wt.d2ivxdt2(i,:) wt.d2ipdt2(i) wt.d2cpdt2(i,:) wt.d2bdefdt2(i,:) wt.d2bdefdt2(i,:)]);   
end
Qhat = blkdiag(Q{:});
Qvhat = blkdiag(Qv{:});
Q2vhat = blkdiag(Q2v{:});
Qbar = Qhat + Sy'*Qvhat*Sy + S2y'*Q2vhat*S2y;


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

  % enforce linear independence
  [u,s,v] = svd_tol(Aeq, 0.9999);
  Aeq = u'*Aeq;
  beq = u'*beq;


  % Inequality constraints
  % ----------------------
  
  % Use this if no inequality constraints:
  %     Aineq = zeros(0,npv);
  %     bineq = zeros(0,1);
  
  Aineq = -eye(N*circ.ncx);
  bineq = -reshape(constraints_min.icx', [], 1);
   
  i = isnan(bineq);
  bineq(i,:) = [];
  Aineq(i,:) = [];
  
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
  cphat = yhat(circ.iipx+1:circ.iipx+ngaps, :);
  
  prevsoln = variables2struct(xhat, icxhat, ivxhat, iphat, yhat);

  end


close all

figure
hold on
% plot(t, [pla(:).icx]/1e3, 'g')
plot(t, icxhat(2:end,:)/1e3, '-b')
icx_efit = [efit01_eqs.gdata(:).icx];
ivx_efit = [efit01_eqs.gdata(:).ivx];
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

if 0
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
end
%%

if 0
  close all

  psizr_app = mpc*icxhat + mpv*ivxhat;
  % psizr_app = mpc*[efit01_eqs.gdata(:).icx] + mpv*[efit01_eqs.gdata(:).ivx];


  psizr_app = reshape(psizr_app, nz, nr, []);
  psizr_pla = [pla(:).psizr_pla];
  psizr_pla = reshape(psizr_pla, nz, nr, []);
  psizr = psizr_app + psizr_pla;

  f = figure;
  f.Position =  [731 485 437 692];
  for i = 1:N
    clf
    psi = squeeze(psizr(:,:,i));
    eq = eq_params(psi, tok_data_struct, 1);

    target = targets_array(i);
    scatter(target.rcp, target.zcp, 'k', 'filled');
    scatter(target.rbdef, target.zbdef, 120, 'ob')
    title(i)

    drawnow
  end

end













%%
if 0

  % DEBUGGING: gsdesign comparison
  
  close all
  
  i = 35;
  
  icx = icxhat(:,i);         
  ivx = ivxhat(:,i);
  ip = iphat(i);
%   icx = efit01_eqs.gdata(i).icx;
%   ivx = efit01_eqs.gdata(i).ivx;
%   ip = efit01_eqs.gdata(i).cpasma;
  
  init = efit01_eqs.gdata(i);
  % init = copyfields(pla(i), efit01_eqs.gdata(i), [], false);
  
  opts.pres = efit01_eqs.gdata(i).pres;
  opts.fpol = efit01_eqs.gdata(i).fpol;
  
  % opts.pres = pla(i).pres;
  % opts.fpol = pla(i).fpol;
  
  x = [icx; ivx; ip];
  [spec, init, config] = make_gsdesign_inputs2(x, tok_data_struct, init, circ, opts);
  
  
  % % DEBUG
  spec.weights.pres = ones(size(opts.pres)) * 1e-10;
  spec.weights.fpol = ones(size(opts.fpol)) * 1e-10;
  spec.targets.li = targets.li(i);
  spec.weights.li = 1;
  % spec.targets.betap = 0.84;
  % spec.weights.betap = 0;
  % read_wmhd(eq, tok_data_struct)
  % read_wmhd(pla(i), tok_data_struct)
  % read_wmhd(init, tok_data_struct)
  spec.targets.Wth = targets.wmhd(i);
  spec.weights.Wth = 1;
  
  
  eq = gsdesign(spec,init,config);
  
  spec.targets.Wth
  eq.Wth
  
  
  
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
end















































