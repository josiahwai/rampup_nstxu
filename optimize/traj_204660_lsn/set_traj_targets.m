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


% load geometry
vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
tok_data_struct.imks = 1;
circ = nstxu2016_circ(tok_data_struct);
struct_to_ws(tok_data_struct);

% =========================
% Parameters: Ip, W, li, Rp
% =========================
ip_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:IPMEAS');
wmhd_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:WMHD');
li_sig = mds_fetch_signal(shot, 'efit01', [], '.RESULTS.AEQDSK:LI');
[rp_sig.sigs, rp_sig.times] = load_rp_profile(0);
params = variables2struct(ip_sig, wmhd_sig, li_sig, rp_sig);

% ===============================
% Load reference data for targets
% ===============================
ref.t = tref;

opts.cache_dir = [ROOT '/fetch/cache/'];
efit01_eqs = fetch_eqs_nstxu(shot, tref, 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts);

% boundary definition
bry_opts.plotit = 0;
for i = 1:Nref
  bry(i) = eq_bdef_analysis(efit01_eqs.gdata(i), tok_data_struct, bry_opts);
end

fds = {'islimited', 'rx_lo', 'zx_lo', 'rx_up', 'zx_up', 'rtouch', 'ztouch', 'rbdef', 'zbdef'};
for i = 1:length(fds)
  y = [bry(:).(fds{i})];
  ref.(fds{i}) = fillmissing(y, 'nearest');
end

% control points
gap_opts.plotit = 0;
gap_opts.use_out_up_lo = 0;
for i = 1:Nref
  gaps(i) = get_nstxu_gaps(efit01_eqs.gdata(i), gap_opts);
end
ref.rcp = [gaps(:).r]';
ref.zcp = [gaps(:).z]';


%%
% ========
% Targets
% ========
targets.time = t'; 

fds = {'islimited', 'rx_lo', 'zx_lo', 'rx_up', 'zx_up', 'rtouch', ...
  'ztouch', 'rbdef', 'zbdef', 'rcp', 'zcp'};

for i = 1:length(fds)  
  fd = fds{i};
  targets.(fd) = interp1(tref, ref.(fd), t(:));
end

k = (targets.islimited==1);
targets.islimited = k; 
targets.rbdef = [targets.rtouch(k); targets.rx_lo(~k)];
targets.zbdef = [targets.ztouch(k); targets.zx_lo(~k)];

targets.rcp = smoothdata(targets.rcp);
targets.zcp = smoothdata(targets.zcp);
targets.ip = smooth(interp1(ip_sig.times, ip_sig.sigs, t));






























