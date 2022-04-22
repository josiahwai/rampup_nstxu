clear all; clc; close all
warning('off','MATLAB:polyshape:repairedBySimplify')
warning('off', 'curvefit:fit:nonDoubleYData')
warning('off', 'stats:pca:ColRankDefX')
ROOT = getenv('RAMPROOT');


% settings
shot = 204660;  
tref = [.07 .22 .4 .9]';
Nref = length(tref);
N = 50;
t = linspace(.07, .9, N)';
ts = mean(diff(t));

opts.cache_dir = [ROOT '/fetch/cache/'];
efit01_eqs = fetch_eqs_nstxu(shot, t, 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts);
icx_efit = [efit01_eqs.gdata(:).icx];
ivx_efit = [efit01_eqs.gdata(:).ivx];

% load geometry
vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
tok_data_struct.mpp_full = tok_data_struct.mpp;
tok_data_struct.mpp = load('nstxu_obj_config2016_6565.mat').tok_data_struct.mpp;
tok_data_struct.imks = 1;
circ = nstxu2016_circ(tok_data_struct);
iy = circ.iy;
struct_to_ws(tok_data_struct);


% =========================
% Parameters: Ip, W, li, Rp
% =========================
ip_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:IPMHD');
wmhd_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:WMHD', 0);
li_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:LI', 0);
[rp_sig.sigs, rp_sig.times] = load_rp_profile(0);
params = variables2struct(ip_sig, wmhd_sig, li_sig, rp_sig);

% ===============================
% Load reference data for targets
% ===============================
refs.time = tref;

opts.cache_dir = [ROOT '/fetch/cache/'];
ref_eqs = fetch_eqs_nstxu(shot, tref, 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts);
eq0 = ref_eqs.gdata(1);

refs.li = interp1(li_sig.times, li_sig.sigs, refs.time);
refs.wmhd = interp1(wmhd_sig.times, wmhd_sig.sigs, refs.time);
refs.ip = interp1(ip_sig.times, ip_sig.sigs, refs.time);


refs.icx = nan(Nref, circ.ncx);
refs.icx(1,:) = ref_eqs.gdata(1).icx;

refs.icx(tref<=0.4, iy.PF2U) = 0;
refs.icx(tref<=0.4, iy.PF2L) = 0;

refs.islimited = tref < 0.23;

refs.ivx = nan(Nref, circ.nvx);
refs.ivx(1,:) = ref_eqs.gdata(1).ivx;


% shape targets
refs.rbbbs = {ref_eqs.gdata(:).rbbbs}';
refs.zbbbs = {ref_eqs.gdata(:).zbbbs}';

refs.rbbbs{2} = refs.rbbbs{3} - (min(refs.rbbbs{3}) - min(limdata(2,:)));
refs.zbbbs{2} = refs.zbbbs{3}; 

refs.rbbbs{4} = refs.rbbbs{3} - (min(refs.rbbbs{3}) - min(refs.rbbbs{4}));
refs.zbbbs{4} = refs.zbbbs{3}; 


for i = 1:Nref
  gaps(i) = get_nstxu_gaps_from_bbbs(refs.rbbbs{i}, refs.zbbbs{i});
end
refs.rcp = [gaps(:).r]';
refs.zcp = [gaps(:).z]';

refs_array = struct2structarray(refs);

clear eqs
for i = 1:Nref
  ref = refs_array(i);
  eqs(i) = gsdesign_fit(ref, tok_data_struct);
end

% =================
% Plasma parameters
% =================

fds = {'pcurrt', 'psizr_pla'};

for i = 1:length(fds)
  fd = fds{i};
  x = [eqs(:).(fd)];
  x = reshape(x, nz, nr, []);
  x = permute(x, [3 1 2]);
  x = interp1(tref, x, t);
  x = smoothdata(x, 1, 'sgolay', 15, 'degree', 2);
  pla.(fd) = x;
  for j = 1:N
    pla_array(j).(fd) = squeeze(x(j,:,:));
  end
end


% ========
% Targets
% ========
init = ref_eqs.gdata(1);
ngaps = size(refs.rcp,2);

targets.time = t(:);
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


targets.rcp = smoothdata(interp1(tref, refs.rcp, t(:)));
targets.zcp = smoothdata(interp1(tref, refs.zcp, t(:)));
targets.rtouch(1:N,1) = min(limdata(2,:));
targets.ztouch(1:N,1) = 0;
targets.islimited= t < 0.23;

targets.rx_lo = smooth(interp1([0.05 0.1 0.23 0.4 0.9], [0.56 0.59 0.595 0.64 0.66], t));
targets.zx_lo = smooth(interp1([0.05 0.1 0.23 0.4 0.9], -[1.14 1.09 1.055 1.055 1.055], t));

targets.rx_up(1:N,1) = 0.5;
targets.zx_up(1:N,1) = 0;

k = targets.islimited;
targets.rbdef = [targets.rtouch(k); targets.rx_lo(~k)];
targets.zbdef = [targets.ztouch(k); targets.zx_lo(~k)];

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

%%
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
icx = [eqs(:).icx];
icx = interp1(tref, icx', t);
x = timeseries(icx', t);
opts.plotit = 0;
y = smooth_ts(x, opts);
icx = squeeze(y.Data);

close all
plot(t, icx)

pcurrts = [eqs(:).pcurrt];
pcurrts = reshape(pcurrts, nz, nr, []);
pcurrts = permute(pcurrts, [3 1 2]);
pcurrts = interp1(tref, pcurrts, t);

% contour(squeeze(pcurrts(30,:,:)))



































