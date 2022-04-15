clear all; clc; close all
warning('off','MATLAB:polyshape:repairedBySimplify')
warning('off', 'curvefit:fit:nonDoubleYData')
warning('off', 'stats:pca:ColRankDefX')
ROOT = getenv('RAMPROOT');


% settings
shot = 204660;  
tref = [.07 0.22 0.4 .9];
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
tfetch = [t(1) t(end) t(end)]; 
ip_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:IPMHD', 0);
wmhd_sig = mds_fetch_signal(shot, 'efit01', tfetch, '.RESULTS.AEQDSK:WMHD', 0);
li_sig = mds_fetch_signal(shot, 'efit01', tfetch, '.RESULTS.AEQDSK:LI', 0);
[rp_sig.sigs, rp_sig.times] = load_rp_profile(0);

tflattop = 0.4;
wmhd_sig.times = [t(1) tflattop t(end)];
li_sig.times = [t(1) tflattop t(end)];


% ===============================
% Load reference data for targets
% ===============================
refs.time = tref';

opts.cache_dir = [ROOT '/fetch/cache/'];
eq0 = fetch_eqs_nstxu(shot, tref(1), 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts).gdata;
eqf = fetch_eqs_nstxu(shot, tref(end), 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts).gdata;
init = eq0;

% boundary definition
bry_opts.plotit = 0;
bry0 = eq_bdef_analysis(eq0, tok_data_struct, bry_opts);
bryf = eq_bdef_analysis(eqf, tok_data_struct, bry_opts);

refs.li = interp1(li_sig.times, li_sig.sigs, refs.time);
refs.wmhd = interp1(wmhd_sig.times, wmhd_sig.sigs, refs.time);
refs.ip = interp1(ip_sig.times, ip_sig.sigs, refs.time);

[refs.rcp(1,:), refs.zcp(1,:)] = interparc(bry0.rbbbs, bry0.zbbbs, 20, 0, 0);
[refs.rcp(4,:), refs.zcp(4,:)] = interparc(bryf.rbbbs, bryf.zbbbs, 20, 0, 0);

ds = min(bryf.rbbbs) - min(tok_data_struct.limdata(2,:));  % distance from plasma to inner wall

refs.rcp(2,:) = refs.rcp(4,:) - ds;
refs.zcp(2,:) = refs.zcp(4,:);

refs.rcp(3,:) = refs.rcp(4,:) - ds*0.3;
refs.zcp(3,:) = refs.zcp(4,:);


refs.rtouch = ones(Nref,1)*bry0.rtouch;
refs.ztouch = ones(Nref,1)*bry0.ztouch;
refs.rx_lo = bryf.rx_lo - [ds ds ds*0.3 0]';
refs.zx_lo = ones(Nref,1)*bryf.rx_lo;
refs.rx_up = bryf.rx_up - [ds ds ds*0.3 0]';
refs.zx_up = ones(Nref,1)*bryf.rx_up;

refs.rbdef = [refs.rtouch(1:2); refs.rx_lo(3:4)];
refs.zbdef = [refs.ztouch(1:2); refs.zx_lo(3:4)];
refs.islimited = [1 1 0 0]';



%%
% ========
% Targets
% ========
init = eq0;
ngaps = size(refs.rcp,2);

targets.time = t';
targets.icx = zeros(N, circ.ncx);
targets.ivx = zeros(N, circ.nvx);
targets.ip = interp1(ip_sig.times, ip_sig.sigs, t(:));
targets.diff_psicp_psitouch = zeros(N, ngaps);
targets.diff_psicp_psixlo = zeros(N, ngaps);
targets.diff_psicp_psixup = zeros(N, ngaps);
targets.psixlo_r = zeros(N,1);
targets.psixlo_z = zeros(N,1);
targets.psixup_r = zeros(N,1);
targets.psixup_z = zeros(N,1);



fds = {'rcp', 'zcp', 'rtouch', 'ztouch', 'rx_lo', 'zx_lo', 'rx_up', ...
  'zx_up', 'rbdef', 'zbdef', 'islimited'};

for i = 1:length(fds)  
  fd = fds{i};
  targets.(fd) = interp1(tref, refs.(fd), t(:));
end

k = (targets.islimited > 0.98);
targets.islimited = k; 

% rx_lo = mds_fetch_signal(shot, 'efit01', t, '.RESULTS.AEQDSK:RXPT1', 0);
% targets.rx_lo(~k) = smooth(rx_lo.sigs(~k));
targets.rx_lo = smooth(interp1([0.05 0.22 0.32 0.53 0.9], [0.56 0.57 0.65 0.66 0.69], t));
% targets.rx_lo = smooth(targets.rx_lo);

% zx_lo = mds_fetch_signal(shot, 'efit01', t, '.RESULTS.AEQDSK:ZXPT1', 0);
% targets.zx_lo(~k) = smooth(zx_lo.sigs(~k));
targets.zx_lo = smooth(interp1([0.06 0.22 0.34 0.5 0.9], -[1.14 1.14 1.07 1.02 1.02], t));
% targets.zx_lo = smooth(targets.zx_lo);

targets.rbdef = [targets.rtouch(k); targets.rx_lo(~k)];
targets.zbdef = [targets.ztouch(k); targets.zx_lo(~k)];

targets.rcp = smoothdata(targets.rcp);
targets.zcp = smoothdata(targets.zcp);

targets = struct_fields_to_double(targets);


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

% opts.cache_dir = [ROOT '/fetch/cache/'];
% efit01_eqs = fetch_eqs_nstxu(shot, t, 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts);
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

% indices used for the controlled variables (cv):
clear cv
cv.icx = circ.iicx;
cv.ivx = circ.iivx;
cv.ip = circ.iipx;
cv.diff_psicp_psitouch = cv.ip + (1:ngaps);
cv.diff_psicp_psixlo = cv.diff_psicp_psitouch(end) + (1:ngaps);
cv.diff_psicp_psixup = cv.diff_psicp_psixlo(end) + (1:ngaps);
cv.psixlo_r = cv.diff_psicp_psixup(end) + 1;
cv.psixlo_z = cv.psixlo_r + 1;
cv.psixup_r = cv.psixlo_z + 1;
cv.psixup_z = cv.psixup_r + 1;
cv.names = fieldnames(cv);

% initialize weights to zero
for i = 1:length(cv.names)
  fd = cv.names{i};
  wts.(fd) = zeros(N, length(cv.(fd)));   % will hold weights on the value of that variable
  dwts.(fd) = zeros(N, length(cv.(fd)));  % weights on the derivative/velocity
  d2wts.(fd) = zeros(N, length(cv.(fd))); % weights on the 2nd derivative
end

% Populate the weights
% ....................

wts.icx = ones(N,circ.ncx) * 3e-5; 
% wts.ivx = ones(N,circ.nvx) * 1e-6;
% activation = sigmoid(t, 30, 0.25, 1);
% activation = interp1([0 0.14 0.25 0.9], [0 0 1 1], t(:), 'linear');
activation = double(~targets.islimited);

wts.ip = ones(N,circ.np) * 3e-5;

wts.diff_psicp_psitouch = (1-activation) * ones(1,ngaps) * 2e8;
wts.diff_psicp_psixlo = activation * ones(1,ngaps) * 2e8;
wts.diff_psicp_psixup = ones(N, ngaps) * 0;

wts.psixlo_r = activation * 1e8;
wts.psixlo_z = activation * 1e8;
wts.psixup_r = activation * 0;
wts.psixup_z = activation * 0;

dwts.icx = ones(N,circ.ncx) / ts^2 * 1e-7;
d2wts.icx = ones(N, circ.ncx) / ts^4 * 2e-10;


% rearrange data for easier access later
targets_array = struct2structarray(targets);
refs_array = struct2structarray(refs);
wts_array = struct2structarray(wts);
dwts_array = struct2structarray(dwts);
d2wts_array = struct2structarray(d2wts);


%%
% =================
% Plasma parameters
% =================

% solve free-bry GS to get equilibria
% this version fits to: li, Wmhd, Ip, and boundary
for i = 1:Nref
  i  
  ref = refs_array(i);
  opts.init = eq0; % ref_eqs.gdata(i);
  opts.plotit = 0;
  opts.max_iterations = 10;
  pla_array(i) = semifreegs(ref, tok_data_struct, opts);
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
  pla.psizr_pla = smoothdata(pla.psizr_pla, 1, 'sgolay', 5);
end
























































