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

targets.ip = smooth(interp1(ip_sig.times, ip_sig.sigs, t));

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

refs.rbbbs = {ref_eqs.gdata(:).rbbbs}';
refs.zbbbs = {ref_eqs.gdata(:).zbbbs}';

refs.icx = nan(Nref, circ.ncx);
refs.icx(1,:) = ref_eqs.gdata(1).icx;
% refs.icx(4,:) = ref_eqs.gdata(4).icx;

refs.icx(tref<=0.4, iy.PF2U) = 0;
refs.icx(tref<=0.4, iy.PF2L) = 0;

% refs.icx(:,1) = 0;

refs.islimited = tref < 0.23;

refs.ivx = nan(Nref, circ.nvx);
refs.ivx(1,:) = ref_eqs.gdata(1).ivx;

refs_array = struct2structarray(refs);

clear eqs
for i = 1:Nref
  ref = refs_array(i);
  eqs(i) = gsdesign_fit(ref, tok_data_struct);
end


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


%%

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
mpp = tok_data_struct.mpp_full;


% Estimate time-dependent parameters
for i = 1:N
  pcurrt = squeeze(pcurrts(i,:,:));
  pcurrt = pcurrt(:);
  ip = sum(pcurrt(:));
  params.mcIp(i,:) = mpc' * pcurrt / ip;
  params.mvIp(i,:) = mpv' * pcurrt / ip;
  params.Lp(i,:) = pcurrt' * mpp * pcurrt / ip^2;  
end

params.Rp = interp1(rp_sig.times, rp_sig.sigs, t);

%%
ipfx = icx(2:end,:);
ipfxdot = gradient(ipfx, ts);
params.mpfIp = params.mcIp(:,2:end);
psidot_pf = diag(params.mpfIp * ipfxdot);

Mip_oh = params.mcIp(:,1);

ip = targets.ip;
ipdot = gradient(targets.ip, ts);

iohdot = -(params.Rp.*ip + psidot_pf + params.Lp.*ipdot) ./ Mip_oh;
ioh = eq0.icx(1) + cumtrapz(t, iohdot);

figure
hold on
plot(t, ioh)
plot(t, icx_efit(1,:))
%%

pcurrts2 = reshape(pcurrts, nz*nr, N);
u1 = mpv' * mpp*pcurrts2;

icx(1,:) = ioh;
icxdot = gradient(icx, ts);
u2 = mvc * icxdot;

u = u1 + u2;

A = -inv(mvv) * diag(Rv);
B = -inv(mvv);
C = eye(circ.nvx);
sys = ss(A,B,C,0);

ivx0 = eq0.ivx;

ivx = lsim(sys, u, t, ivx0);

%%

% update refs with ioh and ivx
refs.icx(:,1) = interp1(t, ioh, tref);
refs.ivx = interp1(t, ivx, tref);

refs_array = struct2structarray(refs);


for i = 1:Nref
  ref = refs_array(i);
  eqs(i) = gsdesign_fit(ref, tok_data_struct);
end

%%
icx = [eqs(:).icx];
icx = interp1(tref, icx', t);
x = timeseries(icx', t);
opts.plotit = 0;
y = smooth_ts(x, opts);
icx = squeeze(y.Data);

icxhat = icx;
ivxhat = ivx';
iphat = [efit01_eqs.gdata(:).cpasma];

%%
figure
hold on
plot(t, icxhat, 'b')
plot(t, icx_efit, '--r')


figure
hold on
plot(t, ivxhat, 'b')
plot(t, ivx_efit, '-r')



%%


if 0

  % GSDESIGN COMPARISON
  
  close all   
  
  i = 40;
  
  icx = icxhat(:,i);    
  ivx = ivxhat(:,i);
  ip = iphat(i);
  %   icx = efit01_eqs.gdata(i).icx;
  %   ivx = efit01_eqs.gdata(i).ivx;
  %   ip = efit01_eqs.gdata(i).cpasma;   

  init = efit01_eqs.gdata(i);  
  opts.pres = efit01_eqs.gdata(i).pres;
  opts.fpol = efit01_eqs.gdata(i).fpol;

  % init = copyfields(pla(i), efit01_eqs.gdata(i), [], false);  
  % opts.pres = pla(i).pres;
  % opts.fpol = pla(i).fpol;
  
  x = [icx; ivx; ip];
  [spec, init, config] = make_gsdesign_inputs2(x, tok_data_struct, init, circ, opts);
  
  
  % Modify weights
  spec.weights.pres = ones(size(opts.pres)) * 1e-10;
  spec.weights.fpol = ones(size(opts.fpol)) * 1e-10;
    [~,~,~,li] = inductance(efit01_eqs.gdata(i), tok_data_struct);
    wmhd = read_wmhd(efit01_eqs.gdata(i), tok_data_struct);
%   li = pla_array(i).li;
%   wmhd = read_wmhd(pla_array(i), tok_data_struct); 
  spec.targets.li = li;
  spec.weights.li = 1;
  spec.targets.Wth = wmhd;
  spec.weights.Wth = 1;
  
  spec.weights.sep(1:end) = 0.1;


  spec.targets.zcur = -0.01;
  spec.weights.zcur = 1;
  ipf3u = 7:10;
  ipf3l = 21:24;
  spec.targets.ic(ipf3u) = spec.locks.ic(ipf3u);
  spec.targets.ic(ipf3l) = spec.locks.ic(ipf3l);
  spec.weights.ic(ipf3l) = 1e-8;
  spec.weights.ic(ipf3u) = 1e-8;
  spec.locks.ic(ipf3u) = nan;
  spec.locks.ic(ipf3l) = nan;


  eq = gsdesign(spec,init,config);
  
  spec.targets.Wth
  eq.Wth
   
  eq.ic(ipf3u(1))
  spec.targets.ic(ipf3u(1))
  
  eq.ic(ipf3l(1))
  spec.targets.ic(ipf3l(1))
  

%   figure
%   hold on
%   eqt = efit01_eqs.gdata(i);
%   contour(rg, zg, eqt.psizr, [eqt.psibry eqt.psibry], '--b')
%   contour(rg, zg, eq.psizr, [eq.psibry eq.psibry], '-r')
%   plot_eq(eq)
%   contour(rg, zg, eqt.psizr, [eqt.psibry eqt.psibry], '--b')
%   % contour(rg, zg, init.psizr, [init.psibry init.psibry], '--b')
%   scatter(targets.rcp(i,:), targets.zcp(i,:), 'k', 'filled')
%   plot(targets.rbdef(i), targets.zbdef(i), 'xb', 'markersize', 20, 'linewidth', 3)
%   l = legend('EFIT01', 'optimizer', 'fontsize', 14);  
%   set(gcf, 'Position', [634 449 367 529])
end



































