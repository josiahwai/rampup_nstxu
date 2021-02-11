% EXAMPLE:
% shot = 204660;
% tspan0 = [30:10:960] / 1e3;
% tspanf = tspan0;
% modeldir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/std/';
% eqdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk_import/';
% load('nstxu_obj_config2016_6565.mat')
% circ = nstxu2016_circ(tok_data_struct);
% [traj, eqs] = load_interp_trajectory(shot, circ, tspan0, tspanf, eqdir, modeldir)
% plot(traj.tspan, traj.ic)

% =========================================
% Load eq snapshots to make shot trajectory
% =========================================
function [traj, eqs] = load_interp_trajectory4(shot, circ, tspan, eqdir, modeldir)

%   ccc
%   load('matlab.mat')
%   tspan = (30:10:960)'/1e3;
   
  
  dt = mean(diff(tspan)); 
  N = length(tspan);
  
  for i = 1:N
    t = tspan(i) * 1000;
    
    load([modeldir num2str(shot) '_' num2str(t) '_sys.mat']);
    eq = load([eqdir 'eq' num2str(shot) '_' num2str(t) '.mat']);
    
    eq = eq.eq;    
    eqs{i} = eq;
    rxx(i,:) = sys.rxx;
    A(i,:,:) = sys.A;
    B(i,:,:) = sys.B;    
    psibry(i) = eq.psibry;
    li(i) = eq.li;
    betap(i) = eq.betap;
    pres(i,:) = eq.pres;
    fpol(i,:) = eq.fpol;
    zcur(i) = eq.zcur;
    rcur(i) = eq.rcur;
    ic(i,:) = circ.Pcc \ eq.ic;  % eq stores as ungrouped, convert to circuit
    iv(i,:) = circ.Pvv \ eq.iv;
    ip(i,:) = eq.cpasma;
  end
  x = [ic iv ip];
  
  %%
  % Use smoothing splines to smooth the A matrix
  Asm = smoothdata(A, 1, 'movmedian', 10);  % removes outliers    
  y = squeeze(reshape(Asm, N, 1, [])); 
  for k = 1:numel(sys.A)
    f = fit(tspan(:), y(:,k), 'smoothingspline', 'SmoothingParam', 0.9999);
    y(:,k) = f(tspan);
  end  
  Asm = reshape(y, N, circ.nx, circ.nx);
  
  % Calculate the smoothed lstar and smoothed B matrix
  for i = 1:N
    Adum = squeeze(Asm(i,:,:));    
    lstarinv = -Adum * diag(1./ rxx(i,:));    
    Bsm(i,:,:) = lstarinv(:, circ.iicx);
    lstari(i,:,:) = lstarinv;
  end
  
  % stabilize and discretize A matrix
  for i = 1:N
    Astab = numerically_stabilize( squeeze(Asm(i,:,:)), 1e5);
    [Ad(i,:,:), Bd(i,:,:)] = c2d(Astab, squeeze(Bsm(i,:,:)), dt);
    e(i,:) = eig(squeeze(Ad(i,:,:)));
  end
  
%   PLOTS:   
%   ad = squeeze(reshape(Ad, N, 1, [])); 
%   figure
%   plot(tspan,ad)
%   
%   bd = squeeze(reshape(Bd, N, 1, [])); 
%   figure
%   plot(tspan,bd)
%   
%   figure
%   bar(tspan,  max(abs(e')) - 1)
  
  
  traj = variables2struct(tspan, Ad, Bd, lstari, psibry, li, betap, pres, fpol, zcur, ...
    rcur, ic, iv, ip, x); 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

