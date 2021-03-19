% Vessel Resistance Fitting Algorithm: fits the vessel elements
ccc
RAMPROOT = getenv('RAMPROOT');

% ========
% Settings
% ========
load_data_from_mdsplus = 0;
enforce_stability = 1;

load([RAMPROOT '/sysid/vessel_sysid2/train20_shotdata.mat']);

saveit = 0;

% shotlist = [204660];
shotlist = [204660 203012 204330 204146 204659 204186 ...
  204963 204328 204092 204324 204306 204074 ...
  203018 204651 203321 203502 203849 204062 ...
  204307 204114];

nshots = length(shotlist);
starttimes = zeros(nshots,1);
endtimes = [0.85 0.55 0.6 1.7 0.9 0.3 1.3 0.5 1.2 0.75 1 0.25 ...
  0.5 1.6 0.9 0.4 1.2 1.8 0.9 0.35]';

% =========
% LOAD DATA
% =========

if load_data_from_mdsplus % load data afresh
  for ishot = 1:nshots
    
    shot = shotlist(ishot);
    disp(['Fetching shot ' num2str(shot)])
    
    Ts = 1e-3;
    tstart = starttimes(ishot);
    tend = endtimes(ishot);
    tsample = tstart:Ts:tend;
    
    % Coil Currents
    include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
      'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
    
    icsignals = get_icsignals(shot, [], [], include_coils);
    ivsignals = get_ivsignals(shot);
    ipsignals = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:IPMEAS');
    
    icts = timeseries(icsignals.sigs,icsignals.times);
    ivts = timeseries(ivsignals.sigs, ivsignals.times);
    ipts = timeseries(ipsignals.sigs,ipsignals.times);
    
    %     figure
    %     plot(ipts)
    %     title(num2str(shot))
    
    icts = resample(icts,tsample);
    ivts = resample(ivts,tsample);
    ipts = resample(ipts,tsample);
    
    % obtain derivatives
    Tsmooth = 10;  % [ms]
    nsmooth = floor(Tsmooth/1000/Ts);
    
    ic = smoothdata(icts.Data, 1, 'movmean', nsmooth);
    iv = smoothdata(ivts.Data, 1, 'movmean', nsmooth);
    ip = smoothdata(ipts.Data, 1, 'movmean', nsmooth);  
    
    icdot = gradient(ic', Ts)';
    ivdot = gradient(iv', Ts)';
    ipdot = gradient(ip', Ts)';
    
    icdot = smoothdata(icdot,1,'movmean',nsmooth);
    ivdot = smoothdata(ivdot,1,'movmean',nsmooth);
    ipdot = smoothdata(ipdot,1,'movmean',nsmooth);
    
    % DO NOT use the filtered values for y here    
    y = double(ivts.Data);
    tshift = 0; % ms
    nshift = floor(tshift/1000/Ts);
    y(end-nshift+1:end,:) = [];    
    yall{ishot} = y;
    tsample(1:nshift) = [];
    
    % DO NOT filter the voltages used in u here
    u = double([icdot ipdot]);
    u(1:nshift,:) = [];
    uall{ishot} = u;
    
    shotdata{ishot} = iddata(y, u, Ts);
    x0{ishot} = y(1,:)';
  end
  
  merge_shotdata = merge(shotdata{:});
  
  all_shotdata = variables2struct(yall, uall, Ts);
  fn = [RAMPROOT '/sysid/vessel_sysid2/all_shotdata.mat'];
  save(fn, 'all_shotdata')
  
else % load already saved data
  
  struct_to_ws(all_shotdata);
  
  for i = 1:length(yall)    
    shotdata{i} = iddata(yall{i}, uall{i}, Ts);
  end
  merge_shotdata = merge(shotdata{:});
end



%%
% ================
% Initialize Model
% ================
vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);

Mxx = vac_sys.sysid_fits.Mxx;
Rxx = vac_sys.sysid_fits.Rxx;

Mvc = Mxx(circ.iivx, circ.iicx);
Mvv = Mxx(circ.iivx, circ.iivx);
Rv = Rxx(circ.iivx);

dum = load([RAMPROOT '/buildmodel/built_models/mcc/204660_100_sys.mat']);
Mvp = dum.sys.lstar(circ.iivx,end) * 0;

file_args = {enforce_stability};
parameters = {'Rv' Rv; 'Mvc' Mvc; 'Mvv' Mvv; 'Mvp' Mvp};
odefun = 'vessel_dynamics2';
sys = idgrey(odefun, parameters, 'd', file_args, Ts, 'InputDelay', 3);

sys.Structure.Parameters(1).Free = true;   % Rv
sys.Structure.Parameters(2).Free = true;   % Mvc
sys.Structure.Parameters(3).Free = false;  % Mvv
sys.Structure.Parameters(4).Free = true;  % Mvp

f = 100;

sys.Structure.Parameters(1).Minimum = Rv / f;    % Rv
sys.Structure.Parameters(2).Minimum = -1e-4;     % Mvc
sys.Structure.Parameters(3).Minimum = -1e-4;     % Mvv
sys.Structure.Parameters(4).Minimum = -1e-4;     % Mvp

sys.Structure.Parameters(1).Maximum = Rv * f;   % Rv
sys.Structure.Parameters(2).Maximum = 5e-4;     % Mvc
sys.Structure.Parameters(3).Maximum = 5e-4;     % Mvv
sys.Structure.Parameters(4).Maximum = 5e-4;     % Mvp

search_options.Algorithm = 'sqp';
search_options.FunctionTolerance = 5e-6;
search_options.StepTolerance = 5e-6;
search_options.MaxIterations = 30;
search_options.Advanced.TolFun = search_options.FunctionTolerance;
search_options.Advanced.TolX = search_options.StepTolerance;
search_options.Advanced.MaxIter = search_options.MaxIterations;
search_options.Advanced.Algorithm = search_options.Algorithm;

wt = eye(circ.nvx);

opt = greyestOptions('Display', 'on', 'InitialState', 'backcast', ...
  'DisturbanceModel', 'none', 'Focus', 'simulation', ...
  'SearchMethod', 'auto', 'EnforceStability', true, 'OutputWeight', wt);

opt.SearchOptions.MaxIterations = 10;
opt.SearchOptions.Tolerance = 1e-5;
% opt.SearchOptions.StepTolerance = 1e-12;
% opt.SearchOptions.FunctionTolerance = 1e-12;

%%
% ==============
% ESTIMATE MODEL
% ==============
sys_est = greyest(merge_shotdata, sys, opt);

%%

Rv = sys_est.Structure.Parameters(1).Value;
Mvc = sys_est.Structure.Parameters(2).Value;
Mvv = sys_est.Structure.Parameters(3).Value;
Mvp = sys_est.Structure.Parameters(4).Value;

vess_fits = variables2struct(Rv, Mvc, Mvv, Mvp);

if saveit
  % save([RAMPROOT '/sysid/vessel_sysid2/fits/vess_fits_BEST.mat'], 'vess_fits')
  % save([RAMPROOT '/sysid/vessel_sysid2/fits/sys_est_BEST.mat'], 'sys_est')
end



for k = 1:nshots
  compare_opt = compareOptions;
  ypred = compare(shotdata{k}, sys_est, compare_opt).OutputData;
  yinit = compare(shotdata{k}, sys, compare_opt).OutputData;
  ytrue = yall{k};
  
  tstart = starttimes(k);
  tend = endtimes(k);
  tsample = tstart:Ts:tend;
  
  figure
  hold on
  plot(tsample, ypred, 'b')
  plot(tsample, ytrue, '--r')
  % plot(yinit, 'k')
end



k = 3;
tstart = starttimes(k);
tend = endtimes(k);
tsample = tstart:Ts:tend;
ypred = compare(shotdata{k}, sys_est, compare_opt).OutputData;
ytrue = yall{k};

for i = 1:circ.nvx
  figure
  hold on
  plot(tsample, ytrue, 'Color', [1 1 1]*0.8)
  plot(tsample, ypred(:,i), 'b', 'linewidth', 2)
  plot(tsample, ytrue(:,i), '--r', 'linewidth', 2)
  title(circ.vvnames{i}, 'fontsize', 18)
  xlim([0 0.85])
end
  






























































