% Vessel Fitting Algorithm: 
% simultaneously fits vessel resistances, mutual inductances with coils & 
% other vessel elements. This can take several hours to run! Suggest only
% fitting to 2 or 3 shots as a test, before fitting to more. 

ccc
ROOT = getenv('RAMPROOT');

savedir = [ROOT '/sysid/fit_coils_vessels/'];
savefn = 'coil_vessel_fit.mat';
load_data_from_mdsplus = 0;
use_previous_fit = 1;  % don't run a new fit, just make some plots with an old fit

saveit = 0;

Ts = 1e-2;

% SETTINGS FOR FITTING TO 10 SHOTS
% shotlist = [204077 202930 203834 204660 204146 203868 203322 203324 203716 203890];
% nshots = length(shotlist);
% starttimes = 0.07 * ones(nshots,1);
% endtimes = [1 0.9 1.2 0.9 1.7 1.6 1 1.2 0.8 1.2];


% SETTINGS FOR FITTING TO 4 SHOTS
shotlist = [204660 204146 204077 203890];
nshots = length(shotlist);
starttimes =[.07 .07 .07 0.07];
endtimes = [0.85 1.7 1 1.2];


tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);

%%
% =========
% LOAD DATA
% =========

if load_data_from_mdsplus
  for ishot = 1:nshots
    
    shot = shotlist(ishot);
    disp(['Fetching shot ' num2str(shot)])
        
    tstart = starttimes(ishot);
    tend = endtimes(ishot);
    tsample = tstart:Ts:tend;
    N = length(tsample);
    
    % Coil Currents
    include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
      'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
    
    icsignals = get_icsignals(shot, [], [], include_coils);
    ivsignals = get_ivsignals(shot);
    ipsignals = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:IPMEAS');
    vsignals = get_vobjcsignals(shot, [], [], include_coils);
    
    icts = timeseries(icsignals.sigs,icsignals.times);
    ivts = timeseries(ivsignals.sigs, ivsignals.times);
    ipts = timeseries(ipsignals.sigs,ipsignals.times);
    vts = timeseries(vsignals.sigs,vsignals.times);
    
    icts = resample(icts,tsample);
    ivts = resample(ivts,tsample);
    ipts = resample(ipts,tsample);
    vts = resample(vts,tsample);
    
    % obtain derivatives
    Tsmooth = 10;  % [ms]
    nsmooth = floor(Tsmooth/1000/Ts);
    
    ic = smoothdata(icts.Data, 1, 'movmean', nsmooth);
    iv = smoothdata(ivts.Data, 1, 'movmean', nsmooth);
    ip = smoothdata(ipts.Data, 1, 'movmean', nsmooth);
    v = smoothdata(vts.Data, 1, 'movmean', nsmooth);
    
    icdot = gradient(ic', Ts)';
    ivdot = gradient(iv', Ts)';
    ipdot = gradient(ip', Ts)';
    
    icdot = smoothdata(icdot,1,'movmean',nsmooth);
    ivdot = smoothdata(ivdot,1,'movmean',nsmooth);
    ipdot = smoothdata(ipdot,1,'movmean',nsmooth);
        
    tree = 'EFIT01';
    tokamak = 'nstxu';
    server = 'skylark.pppl.gov:8501';
    opts.cache_dir = [ROOT '/fetch/cache/'];
    opts.cache_it = 1;
    efit01_eqs = fetch_eqs_nstxu(shot, tsample, tree, tokamak, server, opts);
    
    % calculate inductances
    Lp=[]; mcIp=[]; mvIp=[];
    for i = 1:length(efit01_eqs.time)
      [Lp(i,1), ~, ~, ~, mcIp(i,:), mvIp(i,:)] = inductance(efit01_eqs.gdata(i), tok_data_struct);
    end
    Lp = interp1(efit01_eqs.time, Lp, tsample);
    mcIp = interp1(efit01_eqs.time, mcIp, tsample);
    mvIp = interp1(efit01_eqs.time, mvIp, tsample);
            
    shotdata.(['s' num2str(shot)]) = variables2struct(tsample,ic,iv,ip,v,icdot,ivdot,ipdot,Lp,mcIp,mvIp);
    
  end
  
  save('shotdata.mat', 'shotdata')
else
  load('shotdata.mat')
end

%%
% ==================
% SETUP SYSTEM MODEL
% ==================

% setup the input/output data
for ishot=1:nshots
  shot = shotlist(ishot);

  s = shotdata.(['s' num2str(shot)]);
  s = struct_fields_to_double(s);

  % remove nans
  idxnan = any(isnan(s.ic')) | any(isnan(s.icdot')) | any(isnan(s.mcIp'));
  fns = fieldnames(s);
  for j = 1:length(fns)
    x = s.(fns{j});
    if size(x,1) == 1, x = x'; end
    s.(fns{j}) = x(~idxnan,:);
  end
  
  y = [s.ic s.iv];
  u = [s.v -diag(s.ipdot)*[s.mcIp s.mvIp]];
     
  data{ishot} = iddata(y, u, Ts);  
  
  x0(:,ishot) = y(1,:);
  
end

iodata = merge(data{:});


% initialize model
Pxx = circ.Pxx(1:end-1,1:end-1); % remove Ip circuit from transition map
Rxx = Pxx' * diag([tok_data_struct.resc; tok_data_struct.resv]) * Pxx;
Mxx = Pxx' * [tok_data_struct.mcc tok_data_struct.mcv; tok_data_struct.mcv' tok_data_struct.mvv] * Pxx;

% add the fitted external inductances and resistances from power supplies
ext_fit = load([ROOT '/sysid/fit_Rext_Lext/ext_fit.mat']).ext_fit;
Mxx_ext = diag([ext_fit.Lext; zeros(circ.nvx,1)]);
Rxx_ext = diag([ext_fit.Rext; zeros(circ.nvx,1)]);

Mxx = Mxx + Mxx_ext;
Rxx = diag(Rxx + Rxx_ext);

% separate out the inductances into their constituent parts
mcc = Mxx(circ.iicx, circ.iicx);
lc = diag(mcc);
mcc = mcc - diag(lc);
mcv = Mxx(circ.iicx, circ.iivx);
mvv = Mxx(circ.iivx, circ.iivx);
lv = diag(mvv);
mvv = mvv - diag(lv);
rc = Rxx(circ.iicx);
rv = Rxx(circ.iivx);
mcc_triu = triu(mcc); % since mcc=mcc', we only have to estimate half the parameters
mvv_triu = triu(mvv); 

parameters_init = variables2struct(rc,rv,lc,lv,mcc,mvv,mcv,Mxx,Rxx);

% set up the greybox fit
args.circ = circ;
args.ts = Ts;
file_args = {args};

parameters = {rc, rv, lc, lv, mcc_triu, mvv_triu, mcv};
odefun = 'coil_vessel_dynamics';

order.ny = circ.ncx + circ.nvx;
order.nx = order.ny;
order.nu = size(iodata.InputData{1},2);

sys = idnlgrey(odefun,order,parameters,x0,Ts,'FileArgument',file_args);

sys.Parameters(1).Fixed = false;  % rc
sys.Parameters(2).Fixed = false;  % rv
sys.Parameters(3).Fixed = false;  % lc
sys.Parameters(4).Fixed = false;  % lv
sys.Parameters(5).Fixed = true;  % mcc_triu
sys.Parameters(6).Fixed = true;  % mvv_triu
sys.Parameters(7).Fixed = true;  % mcv

f = 1000;

sys.Parameters(1).Minimum = rc / f;  % rc
sys.Parameters(2).Minimum = rv / f;  % rv
sys.Parameters(3).Minimum = lc / f;  % lc
sys.Parameters(4).Minimum = lv / f;  % lv
sys.Parameters(5).Minimum = -1e-4;  % mcc_triu
sys.Parameters(6).Minimum = -1e-4;  % mvv_triu
sys.Parameters(7).Minimum = -1e-4;  % mcv

sys.Parameters(1).Maximum = rc * f;  % rc
sys.Parameters(2).Maximum = rv * f;  % rv
sys.Parameters(3).Maximum = lc * f;  % lc
sys.Parameters(4).Maximum = lv * f;  % lv
sys.Parameters(5).Maximum = max(mcc(:))*2;  % mcc_triu
sys.Parameters(6).Maximum = max(mvv(:))*2;  % mvv_triu
sys.Parameters(7).Maximum = max(mcv(:))*2;  % mcv

opt = nlgreyestOptions('Display', 'on', 'SearchMethod', 'auto');
opt.SearchOptions.MaxIterations = 20;
opt.SearchOptions.StepTolerance = 1e-8;
opt.SearchOptions.FunctionTolerance = 1e-7;
opt.Regularization.Nominal = 'model';
wt.ic = 5;
wt.iv = 1;
opt.OutputWeight = diag([wt.ic*ones(circ.ncx,1); wt.iv*ones(circ.nvx,1)]);


  
%%
% ==============
% ESTIMATE MODEL
% ==============
if ~use_previous_fit
  
  sys_est = nlgreyest(iodata, sys, opt);
  
  rc = sys_est.Parameters(1).Value;
  rv = sys_est.Parameters(2).Value;
  lc = sys_est.Parameters(3).Value;
  lv = sys_est.Parameters(4).Value;
  mcc_triu = sys_est.Parameters(5).Value;
  mvv_triu = sys_est.Parameters(6).Value;
  mcv = sys_est.Parameters(7).Value;

  mcc = mcc_triu + mcc_triu' + diag(lc);
  mvv = mvv_triu + mvv_triu' + diag(lv);
  Mxx = [mcc mcv; mcv' mvv];
  Rxx = [rc; rv];

  coil_vessel_fit = variables2struct(rc,rv,lc,lv,mcc,mvv,mcv,Mxx,Rxx);
  
  if saveit
    save([savedir savefn], 'coil_vessel_fit')
  end
      
else
  coil_vessel_fit = load([savedir savefn]).coil_vessel_fit;
  struct_to_ws(coil_vessel_fit);
  mcc_triu = triu(mcc - diag(diag(mcc)));
  mvv_triu = triu(mvv - diag(diag(mvv)));
end
 

%% 
% ==============================================
% SIMULATE SOME SHOTS WITH THE FITTED PARAMETERS
% ==============================================

for ishot = 1:nshots

  shot = shotlist(ishot);

  x = x0(:,ishot);
  uall = iodata.InputData{ishot};
  yall = iodata.OutputData{ishot};
  N = size(uall,1);
  t = starttimes(ishot) + Ts*(0:N-1);

  xall = [];
  for i = 1:N
    u = uall(i,:);
    [x,y] = coil_vessel_dynamics(t(i), x, u, rc, rv, lc, lv, mcc_triu, mvv_triu, mcv, file_args);
    xall(i,:) = x;
  end

  figure
  title(['Coil currents: ' num2str(shot)])
  hold on
  plot(t, yall(:,circ.iicx), '--r')
  plot(t, xall(:,circ.iicx), '--b')

  figure
  title(['Vessel currents: ' num2str(shot)])
  hold on
  plot(t, yall(:,circ.iivx), '--r')
  plot(t, xall(:,circ.iivx), '--b')

end


%% 
% ==========================
% PLOT THE FITTED PARAMETERS
% ==========================
ccnames = categorical(circ.ccnames);
vvnames = categorical(circ.vvnames);
p0 = parameters_init;
p1 = coil_vessel_fit;

figure
hold on
title('Coil Resistances')
bar(ccnames, [p0.rc p1.rc])
legend('Init', 'Fit')
ylim([0 5*median(p0.rc)])

figure
hold on
title('Coil Inductances')
bar(ccnames, [p0.lc p1.lc])
legend('Init', 'Fit')

figure
hold on
title('Vessel Resistances')
bar(vvnames, [p0.rv p1.rv])
legend('Init', 'Fit')
ylim([0 5*median(p0.rv)])

figure
hold on
title('Vessel Inductances')
bar(vvnames, [p0.lv p1.lv])
legend('Init', 'Fit')
ylim([0 5*median(p0.lv)])
















































