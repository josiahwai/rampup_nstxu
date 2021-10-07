% Plasma Resistance Fitting Algorithm
ccc
ROOT = getenv('RAMPROOT');


% shotlist = [204660];
shotlist = load('/p/nstxusr/nstx-users/jwai/nstxu-nns/data/matlab/train_test_val_split.mat').split.valshots;
shotlist = shotlist(4:end);
saveit = 1;



load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);

for ishot = 1:length(shotlist)

  shot = shotlist(ishot);

  % =========
  % LOAD DATA
  % =========
  tree = 'EFIT01';
  tokamak = 'nstxu';
  server = 'skylark.pppl.gov:8501';
  opts.cache_it = 1;
  opts.cache_dir = [ROOT '/fetch/cache/'];
  eqs = fetch_eqs_nstxu(shot, 'all', tree, tokamak, server, opts);

  i = eqs.time(:)' <= 0 | gradient([eqs.gdata(:).cpasma]) < -1e5 | [eqs.gdata(:).cpasma] < 1e3;
  eqs.adata(i) = [];
  eqs.gdata(i) = [];
  eqs.tms(i) = [];
  eqs.time(i) = [];
  eqs = struct_fields_to_double(eqs);


  tstart = eqs.time(2);
  tend = eqs.time(end-1);
  N = length(eqs.time)-2;
  tsample = linspace(tstart, tend, N);
  parameter_times = linspace(tstart, tend, min(N,100));
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
  vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
  tok_data_struct = vac_sys.build_inputs.tok_data_struct;
  tok_data_struct.imks = 1;
  circ = nstxu2016_circ(tok_data_struct);
  nr = tok_data_struct.nr;
  nz = tok_data_struct.nz;
  rg = tok_data_struct.rg;
  zg = tok_data_struct.zg;

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
  clear traj
  traj.t = double(eqs.time);
  for i = 1:length(eqs.time)
    pcurrt = eqs.gdata(i).pcurrt(:);
    ip = sum(pcurrt(:));
    traj.mcIp(i,:) = mpc' * pcurrt / ip;
    traj.mvIp(i,:) = mpv' * pcurrt / ip;
    traj.Lp(i,1) = pcurrt' * mpp * pcurrt / ip^2;  
    traj.A(i,:) = polyarea(eqs.gdata(i).rbbbs, eqs.gdata(i).zbbbs);    
  end

  Rp0 = max(1.5e-5-5e-5*parameter_times, 3e-6); % starting guess, [Ohms]
  Lp0 = interp1(traj.t, traj.Lp, parameter_times);

  % fit_Lp = load('fit_Lp.mat').fit_Lp;
  % Lp0 = fit_Lp(parameter_times);

  odefun = 'rp_dynamics';
  order = [1 (circ.ncx+circ.nvx) 1];
  file_args = {traj, parameter_times, dt};
  parameters = {Rp0, Lp0};

  sys = idnlgrey(odefun,order,parameters,x0,dt,'FileArgument',file_args);


  sys.Parameters(1).Fixed = false;  % Rp
  sys.Parameters(2).Fixed = true; % Lp

  sys.Parameters(1).Minimum = 0;  % Rp
  sys.Parameters(2).Minimum = 0;  % Lp



  % ============
  % Fit to data
  % ============
  opt = nlgreyestOptions('Display', 'on', 'SearchMethod', 'auto');

  opt.SearchOptions.MaxIterations = 7;
  opt.SearchOptions.StepTolerance = 1e-8;
  opt.SearchOptions.FunctionTolerance = 1e-7;

  %%
  sys_est = nlgreyest(shotdata, sys, opt);


  %%
  Rp = sys_est.Parameters(1).Value;
  figure
  plot(parameter_times, Rp)
  drawnow

  traj.Rp = interp1(parameter_times, Rp, traj.t);
  traj.eta = traj.Rp .* traj.A;

  res = traj;

  if saveit 
    fn = [ROOT 'sysid/plasma_resistance/fits/res' num2str(shot) '.mat'];
    save(fn, 'res')
  end

end























































