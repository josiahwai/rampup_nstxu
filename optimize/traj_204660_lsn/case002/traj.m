opts.cache_dir = [ROOT '/fetch/cache/'];
efit01_eqs = fetch_eqs_nstxu(shot, t, 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts);


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
  pcurrt = pla_array(i).pcurrt(:);
  ip = sum(pcurrt(:));
  params.mcIp(i,:) = mpc' * pcurrt / ip;
  params.mvIp(i,:) = mpv' * pcurrt / ip;
  params.Lp(i,:) = pcurrt' * mpp * pcurrt / ip^2;  
end

params.Rp = interp1(rp_sig.times, rp_sig.sigs, t);


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
  r = vacuum_response(targets_array(i), tok_data_struct);


  response = [r.disdis; 
              r.dpsicpdis-r.dpsitouchdis; 
              r.dpsicpdis-r.dpsixlodis;
              r.dpsicpdis-r.dpsixupdis;
              r.dpsixlodis_r;
              r.dpsixlodis_z;
              r.dpsixupdis_r;
              r.dpsixupdis_z];

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
  psizr_pla = pla_array(i).psizr_pla;
  psizr = psizr_pla + psizr_app;
  currents = [init.icx; init.ivx; init.cpasma];
  
  ic0hat(:,i) = init.icx;
  x0hat(:,i) = [init.ivx; init.cpasma];
  
  y = measure_y2(psizr, init, targets_array(i), tok_data_struct);
  y0hat(:,i) = struct2vec(y, cv.names);
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
  Q{i} = diag( struct2vec( wts_array(i), cv.names)); 
  Qv{i} = diag( struct2vec( dwts_array(i), cv.names)); 
  Q2v{i} = diag( struct2vec( d2wts_array(i), cv.names));    
end
Qhat = blkdiag(Q{:});
Qvhat = blkdiag(Qv{:});
Q2vhat = blkdiag(Q2v{:});
Qbar = Qhat + Sy'*Qvhat*Sy + S2y'*Q2vhat*S2y;


for iteration = 1:1

  % targets for y in vector form
  rhat = [];
  for i = 1:N    
    rhat = [rhat; struct2vec(targets_array(i), cv.names)];
  end
  % rhat = struct2vec(targets, cv.names);  
  % rhat = [targets.icx targets.ivx targets.ip zeros(size(targets.rcp)) zeros(N,2)]';
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
  psizr_pla = [pla_array(:).psizr_pla];
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
if 1

  % DEBUGGING: gsdesign comparison
  
  % close all   
  
  i = 15;
  
  icx = icxhat(:,i);    
  ivx = ivxhat(:,i);
  ip = iphat(i);
%   icx = efit01_eqs.gdata(i).icx;
%   ivx = efit01_eqs.gdata(i).ivx;
%   ip = efit01_eqs.gdata(i).cpasma;
  
% icx0 = icx;
% icx(iy.PF1AL) = icx0(iy.PF1AU);
% icx(iy.PF1AU) = icx0(iy.PF1AL);
% icx(iy.PF3L) = icx0(iy.PF3U);
% icx(iy.PF3U) = icx0(iy.PF3L);

  

  init = efit01_eqs.gdata(i);  
  opts.pres = efit01_eqs.gdata(i).pres;
  opts.fpol = efit01_eqs.gdata(i).fpol;

  % init = copyfields(pla(i), efit01_eqs.gdata(i), [], false);  
  % opts.pres = pla(i).pres;
  % opts.fpol = pla(i).fpol;
  
  x = [icx; ivx; ip];
  [spec, init, config] = make_gsdesign_inputs2(x, tok_data_struct, init, circ, opts);
  
  
  % % DEBUG
  spec.weights.pres = ones(size(opts.pres)) * 1e-10;
  spec.weights.fpol = ones(size(opts.fpol)) * 1e-10;
%   [~,~,~,li] = inductance(efit01_eqs.gdata(i), tok_data_struct);
%   wmhd = read_wmhd(efit01_eqs.gdata(i), tok_data_struct);
  li = pla_array(i).li;
  wmhd = read_wmhd(pla_array(i), tok_data_struct); 
  spec.targets.li = li;
  spec.weights.li = 1;
  spec.targets.Wth = wmhd;
  spec.weights.Wth = 1;

%   spec.targets.rbdef = targets.rbdef(i);
%   spec.targets.zbdef = targets.zbdef(i);
%   spec.weights.bdef = 1;
  
  spec.weights.sep(1:end) = 0.1;

%   spec.targets.zcur = -0.01;
%   spec.weights.zcur = 1;
%   ipf3u = 7:10;
%   ipf3l = 21:24;
%   spec.targets.ic(ipf3u) = spec.locks.ic(ipf3u);
%   spec.targets.ic(ipf3l) = spec.locks.ic(ipf3l);  
%   spec.weights.ic(ipf3l) = 1e-8;
%   spec.weights.ic(ipf3u) = 1e-8;  
%   spec.locks.ic(ipf3u) = nan;
%   spec.locks.ic(ipf3l) = nan;

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
  plot(targets.rbdef(i), targets.zbdef(i), 'xb', 'markersize', 20, 'linewidth', 3)
  l = legend('EFIT01', 'optimizer', 'fontsize', 14);  
  set(gcf, 'Position', [634 449 367 529])
end

%%
% 
%
% close all
% 
% eq = efit01_eqs.gdata(i);
% psi = eq.psizr;
% psi_app = reshape(mpc*eq.icx, nz, nr);
% psi_pla = psi - psi_app;
% 
% icx = eq.icx;
% icx(iy.PF1AL) = eq.icx(iy.PF1AL) + 1e3;
% psi_app = reshape(mpc*icx, nz, nr);
% psi = psi_pla + psi_app;
% [rx, zx, psix] = isoflux_xpFinder(psi, 0.6, -1.1, rg, zg);
% 
% plot_eq(eq)
% contour(rg, zg, psi, [psix psix], 'b')
% 
% %%
% close all
% 
% 
% i = 12;
% eq = efit01_eqs.gdata(i);
% 
% icx = efit01_eqs.gdata(i).icx;
% ivx = efit01_eqs.gdata(i).ivx;
% ip = efit01_eqs.gdata(i).cpasma;
% 
% icx1 = eq.icx;
% 
% x = [icx1; ivx; ip];
% [spec, init, config] = make_gsdesign_inputs2(x, tok_data_struct, init, circ, opts);
% spec.weights.pres = ones(size(opts.pres)) * 1e-10;
% spec.weights.fpol = ones(size(opts.fpol)) * 1e-10;
% [~,~,~,li] = inductance(efit01_eqs.gdata(i), tok_data_struct);
% wmhd = read_wmhd(efit01_eqs.gdata(i), tok_data_struct);
% spec.targets.li = li;
% spec.weights.li = 1;
% spec.targets.Wth = wmhd;
% spec.weights.Wth = 1;
% eq1 = gsdesign(spec,init,config);
% 
% %%
% icx2 = eq.icx;
% icx2(iy.PF1AU) = icx2(iy.PF1AU) + 1e3;
% 
% x = [icx2; ivx; ip];
% [spec, init, config] = make_gsdesign_inputs2(x, tok_data_struct, init, circ, opts);
% spec.weights.pres = ones(size(opts.pres)) * 1e-10;
% spec.weights.fpol = ones(size(opts.fpol)) * 1e-10;
% [~,~,~,li] = inductance(efit01_eqs.gdata(i), tok_data_struct);
% wmhd = read_wmhd(efit01_eqs.gdata(i), tok_data_struct);
% spec.targets.li = li;
% spec.weights.li = 1;
% spec.targets.Wth = wmhd;
% spec.weights.Wth = 1;
% eq2 = gsdesign(spec,init,config);
% 
% figure
% plot_eq(eq1)
% contour(rg, zg, eq2.psizr, [eq2.psibry eq2.psibry], 'b', 'linewidth', 2)
% 
% 
% %%
% eq = efit01_eqs.gdata(i);
% r = gspert(eq, tok_data_struct);
% 
% dpsizrdis = mpp*r.dcphidis;
% 
% r.dpsizrdix = dpsizrdis * blkdiag(circ.Pcc, circ.Pvv);
% 
% m1 = reshape(mpc(:,iy.PF1AL), nz, nr);
% m2 = reshape(r.dpsizrdix(:,iy.PF1AL), nz, nr);
%   
% figure; plot_nstxu_geo(tok_data_struct); contour(rg, zg, m1); colorbar
% figure; plot_nstxu_geo(tok_data_struct); contour(rg, zg, m2); colorbar
% figure; plot_nstxu_geo(tok_data_struct); contour(rg, zg, m1+m2); colorbar




















