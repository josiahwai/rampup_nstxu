clear all; clc; close all
warning('off','MATLAB:polyshape:repairedBySimplify')
warning('off', 'curvefit:fit:nonDoubleYData')
warning('off', 'stats:pca:ColRankDefX')
ROOT = getenv('RAMPROOT');


% settings
shot = 204660;  
tref = [.07 .22 .4 .9];
Nref = length(tref);
N = 50;
t = linspace(.07, .9, N);
ts = mean(diff(t));

% load geometry
vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
tok_data_struct.imks = 1;
circ = nstxu2016_circ(tok_data_struct);
iy = circ.iy;
struct_to_ws(tok_data_struct);


% =========================
% Parameters: Ip, W, li, Rp
% =========================
ip_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:IPMHD');
wmhd_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:WMHD');
li_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:LI');
[rp_sig.sigs, rp_sig.times] = load_rp_profile(0);
params = variables2struct(ip_sig, wmhd_sig, li_sig, rp_sig);


% ===============================
% Load reference data for targets
% ===============================
refs.time = tref';

opts.cache_dir = [ROOT '/fetch/cache/'];
ref_eqs = fetch_eqs_nstxu(shot, tref, 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts);
init = ref_eqs.gdata(1);

% boundary definition
bry_opts.plotit = 0;
for i = 1:Nref
  bry(i) = eq_bdef_analysis(ref_eqs.gdata(i), tok_data_struct, bry_opts);
end

fds = {'islimited', 'rx_lo', 'zx_lo', 'rx_up', 'zx_up', 'rtouch', 'ztouch', 'rbdef', 'zbdef'};
for i = 1:length(fds)
  y = [bry(:).(fds{i})];
  refs.(fds{i}) = fillmissing(y', 'nearest');
end

% control points
gap_opts.plotit = 0;
gap_opts.use_out_up_lo = 0;
for i = 1:Nref
  gaps(i) = get_nstxu_gaps(ref_eqs.gdata(i), gap_opts);
end
refs.rcp = [gaps(:).r]';
refs.zcp = [gaps(:).z]';

refs.li = interp1(li_sig.times, li_sig.sigs, refs.time);
refs.wmhd = interp1(wmhd_sig.times, wmhd_sig.sigs, refs.time);
refs.ip = interp1(ip_sig.times, ip_sig.sigs, refs.time);


% ========
% Targets
% ========
targets.time = t'; 

fds = {'islimited', 'rx_lo', 'zx_lo', 'rx_up', 'zx_up', 'rtouch', ...
  'ztouch', 'rbdef', 'zbdef', 'rcp', 'zcp'};

for i = 1:length(fds)  
  fd = fds{i};
  targets.(fd) = interp1(tref, refs.(fd), t(:));
end

k = (targets.islimited > 0.98);
targets.islimited = k; 

rx_lo = mds_fetch_signal(shot, 'efit01', t, '.RESULTS.AEQDSK:RXPT1', 0);
targets.rx_lo(~k) = smooth(rx_lo.sigs(~k));

zx_lo = mds_fetch_signal(shot, 'efit01', t, '.RESULTS.AEQDSK:ZXPT1', 0);
targets.zx_lo(~k) = smooth(zx_lo.sigs(~k));

targets.rbdef = [targets.rtouch(k); targets.rx_lo(~k)];
targets.zbdef = [targets.ztouch(k); targets.zx_lo(~k)];

targets.rcp = smoothdata(targets.rcp);
targets.zcp = smoothdata(targets.zcp);
targets.ip = interp1(ip_sig.times, ip_sig.sigs, t(:));

targets = struct_fields_to_double(targets);

targets.icx = zeros(N, circ.ncx);
targets.ivx = zeros(N, circ.nvx);


% =================
% Plasma parameters
% =================


% rearrange data for easier access later
fns = fieldnames(targets);
for i = 1:length(targets.time)  
  for j = 1:length(fns)
    targets_array(i).(fns{j}) = targets.(fns{j})(i,:);
  end
end
fns = fieldnames(refs);
for i = 1:length(refs.time)  
  for j = 1:length(fns)
    refs_array(i).(fns{j}) = refs.(fns{j})(i,:);
  end
end

% solve free-bry GS to get equilibria
% this version fits to: li, Wmhd, Ip, and boundary
for i = 1:Nref
  i
  close all
  ref = refs_array(i);
  opts.init = ref_eqs.gdata(i);
  opts.plotit = 0;
  opts.max_iterations = 10;
  pla_array(i) = semifreegs(ref, tok_data_struct, opts);
  % pla_array(i).wmhd
  % ref.wmhd

end


% interpolate pla estimates at tref onto the finer timebase t
% looks complicated but mostly just wrangling with data types
fns = fieldnames(pla_array);
for i = 1:length(fns)
  for j = 1:Nref
    pla.(fns{i})(j,:,:) = pla_array(j).(fns{i});    
  end
end

pla.islimited = double(pla.islimited);
for i = 1:length(fns)
  pla.(fns{i})(1:N,:,:) = interp1(tref, pla.(fns{i}), t);
end
pla.islimited = pla.islimited > 0.98;

for i = 1:N
  for j = 1:length(fns)
    pla_array(i).(fns{j}) = squeeze(pla.(fns{j})(i,:,:));
  end
end

if 1
  psizr_pla = permute(pla.psizr_pla, [2 3 1]);   
  x = timeseries(psizr_pla, t);
  y = smooth_ts(x);
  psizr_pla = y.Data;  
  for i = 1:N
    pla_array(i).psizr_pla = squeeze(psizr_pla(:,:,i));
  end
  pla.psizr_pla = permute(psizr_pla, [3 1 2]);
end


% ==================================
% Optimizer constraints 
% (use nan to specify no constraint)
% ==================================
init = ref_eqs.gdata(1);
ngaps = size(targets.rcp,2);

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
% wt.bdef(11) = 2e6;
% wt.bdef(12) = 2e7;
% wt.bdef(13) = 4e7;
% wt.bdef(14) = 1e8;


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





















