% This does voltage control at power supplies and simulates response

ccc
ROOT = getenv('RAMPROOT');

load_data_from_mdsplus = 0;
saveit = 0;
Ts = 1e-2;
shot = 204660;
tstart = 0.07;
tend = 0.85;
t = (tstart:Ts:tend)';
N = length(t);

tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);

rp_fit = load([ROOT '/sysid/fit_plasma_resistance/fits/res' num2str(shot) '.mat']).res;
sys = load([ROOT 'sysid/fit_coils_vessels/coil_vessel_fit.mat']).coil_vessel_fit;

% =========
% LOAD DATA
% =========
if load_data_from_mdsplus
  
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
  
  icts = resample(icts,t);
  ivts = resample(ivts,t);
  ipts = resample(ipts,t);
  vts = resample(vts,t);
  
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
  efit01_eqs = fetch_eqs_nstxu(shot, t, tree, tokamak, server, opts);
  
  % calculate inductances
  Lp=[]; mcIp=[]; mvIp=[];
  for i = 1:length(efit01_eqs.time)
    [Lp(i,1), ~, ~, ~, mcIp(i,:), mvIp(i,:)] = inductance(efit01_eqs.gdata(i), tok_data_struct);
  end
  Lp = interp1(efit01_eqs.time, Lp, t);
  mcIp = interp1(efit01_eqs.time, mcIp, t);
  mvIp = interp1(efit01_eqs.time, mvIp, t);
  
  shotdata.(['s' num2str(shot)]) = variables2struct(t,ic,iv,ip,v,icdot,ivdot,ipdot,Lp,mcIp,mvIp);      
  save('shotdata.mat', 'shotdata')
else
  load('shotdata.mat')
end

sdata = shotdata.(['s' num2str(shot)]);

% Interpolate fits/experimental data onto simulation timebase
fns = fieldnames(rp_fit);
t0 = rp_fit.t(:);
for i = 1:length(fns)
  x = rp_fit.(fns{i});
  if size(x,1) == 1, x = x'; end
  rp_fit.(fns{i}) = interp1(t0, x, t);
end

fns = fieldnames(sdata);
t0 = sdata.tsample(:);
for i = 1:length(fns)
  x = sdata.(fns{i});
  if size(x,1) == 1, x = x'; end
  sdata.(fns{i}) = interp1(t0, x, t);
end

Rp = rp_fit.Rp;
mcIp = rp_fit.mcIp;
mvIp = rp_fit.mvIp;
Lp = rp_fit.Lp;
uall = sdata.v; % power supply voltages
x0 = [sdata.ic(1,:) sdata.iv(1,:) sdata.ip(1)]'; 
x = x0;

for i = 1:N
  
  M = [sys.Mxx [mcIp(i,:) mvIp(i,:)]'];
  M = [M; [mcIp(i,:) mvIp(i,:) Lp(i)]];
  R = diag([sys.rc; sys.rv; Rp(i)]);
  
  Minv = inv(M);
  A = -Minv*R;
  B = Minv(:,circ.iicx);  
  A = numerically_stabilize(A,100);
  
  [Ad,Bd] = c2d(A,B,Ts);
  
  u = uall(i,:)'; 
  
  x = Ad*x + Bd*u;
  
  idx = x(circ.ii_unipolar) < 0;
  x(circ.ii_unipolar(idx)) = 0;
  
  xall(i,:) = x;
end


figure
hold on
plot(t, [sdata.ic], '--r')
plot(t, xall(:,circ.iicx), '--b')


figure
hold on
plot(t, [sdata.iv], '--r')
plot(t, xall(:,circ.iivx), '--b')


figure
hold on
plot(t, [sdata.ip], '--r')
plot(t, xall(:,circ.iipx), '--b')














































