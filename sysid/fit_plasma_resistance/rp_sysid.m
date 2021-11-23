% Plasma Resistance Fitting Algorithm
ccc
ROOT = getenv('RAMPROOT');

% shotlist = [204660];
try
  s = load('/p/nstxusr/nstx-users/jwai/nstxu-nns/data/matlab/train_test_val_split.mat').split;
catch
  s = load('/Users/jwai/Research/nstxu-nns/data/matlab/train_test_val_split.mat').split;
end
shotlist = [s.trainshots; s.valshots; s.testshots];


savedir = [ROOT 'sysid/fit_plasma_resistance/fits_all/'];
saveit = 1;
plotit = 0;

load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);

for ishot = 1:length(shotlist)
  
  try
    
    shot = shotlist(ishot);
    disp(shot)
    
    % =========
    % LOAD DATA
    % =========
    tree = 'EFIT01';
    tokamak = 'nstxu';
    server = 'skylark.pppl.gov:8501';
    opts.cache_it = 1;
    opts.cache_dir = [ROOT '/fetch/cache/'];
    eqs = fetch_eqs_nstxu(shot, 'all', tree, tokamak, server, opts);
    
    i = eqs.time(:)' <= 0.02 | gradient([eqs.gdata(:).cpasma]) < -1e5 | [eqs.gdata(:).cpasma] < 1e3;
    eqs.adata(i) = [];
    eqs.gdata(i) = [];
    eqs.tms(i) = [];
    eqs.time(i) = [];
    eqs = struct_fields_to_double(eqs);
    
    tstart = eqs.time(1);
    tend = eqs.time(end);
    N = length(eqs.time);
    tsample = linspace(tstart, tend, N);
    t = linspace(tstart, tend, N);
    dt = mean(diff(tsample));
    
    icx = [eqs.gdata(:).icx];
    ivx = [eqs.gdata(:).ivx];
    ip = [eqs.gdata(:).cpasma];
    
    icts = timeseries(icx', eqs.time);
    ivts = timeseries(ivx', eqs.time);
    ipts = timeseries(ip', eqs.time);
    
    icts = resample(icts,tsample);
    ivts = resample(ivts,tsample);
    ipts = resample(ipts,tsample);
    
    % obtain derivatives
    Tsmooth = 10;  % [ms]
    nsmooth = ceil(Tsmooth/1000/dt);
    
    ic = smoothdata(icts.Data, 1, 'movmean', nsmooth);
    iv = smoothdata(ivts.Data, 1, 'movmean', nsmooth);
    ip = smoothdata(ipts.Data, 1, 'movmean', nsmooth);
    
    icdot = gradient(ic', dt)';
    ivdot = gradient(iv', dt)';
    ipdot = gradient(ip', dt)';
    
    y = double(ipts.Data);
    u = double([icdot ivdot]);
    shotdata = iddata(y, u, dt);
    
    x0 = y(1);
    
    % =======================================
    % estimate inductances, plasma resistance
    % =======================================
    
    % load geometry
    circ = nstxu2016_circ(tok_data_struct);
    nr = tok_data_struct.nr;
    nz = tok_data_struct.nz;
    rg = tok_data_struct.rg;
    zg = tok_data_struct.zg;
    
    % Load fitted vacuum model parameters
    mpc = tok_data_struct.mpc * circ.Pcc;
    mpv = tok_data_struct.mpv * circ.Pvv;
    mpp = tok_data_struct.mpp_full;
    
    % Estimate time-dependent parameters
    clear traj
    traj.t = double(eqs.time);
    for i = 1:length(eqs.time)
      pcurrt = eqs.gdata(i).pcurrt(:);
      ip = sum(pcurrt(:));
      traj.mcIp(i,:) = mpc' * pcurrt / ip;
      traj.mvIp(i,:) = mpv' * pcurrt / ip;
      traj.Lp(i,1) = pcurrt' * mpp * pcurrt / ip^2;
      traj.A(i,:) = polyarea(eqs.gdata(i).rbbbs, eqs.gdata(i).zbbbs);
      traj.vloop(i,1) = -traj.mcIp(i,:)*icdot(i,:)' - traj.mvIp(i,:)*ivdot(i,:)';
      traj.Rp(i,1) = (traj.vloop(i) - traj.Lp(i)*ipdot(i)) / ip;
    end
    
    
    file_args = {traj, t, dt};
    
    % simulate
    x = x0;
    y = zeros(N,1);
    for i = 1:N
      y(i) = x;
      [x,~] = rp_dynamics(t(i), x, u(i,:), traj.Rp, traj.Lp, file_args);
    end
    
    if plotit
      figure
      hold on
      title(shot)
      xlabel('Time [s]')
      ylabel('Ip [A]')
      plot(t,y)
      plot(ipts)
      legend('Sim', 'EFIT')
      
      figure
      hold on
      title(shot)
      xlabel('Time [s]')
      ylabel('Resistance [Ohm]')
      plot(t,traj.Rp)
    
      drawnow
    end        
    
    if saveit
      res = traj;
      fn = [savedir '/res' num2str(shot) '.mat'];
      save(fn, 'res')
    end
    
  catch
    warning(['Error on shot ' num2str(shot)])
  end
end



















































