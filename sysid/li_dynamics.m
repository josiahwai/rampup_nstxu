clear; clc; close all;

% shot = 203898;
shot = 204660;

MDSPLUS_DIR=getenv('MDSPLUS_DIR');
addpath(genpath(MDSPLUS_DIR))

ROOT = getenv('RAMPROOT');

res_fn = [ROOT 'sysid/fit_plasma_resistance/fits_all/res' num2str(shot) '.mat'];
res = load(res_fn).res;
times = res.t;

res.Rp(1:2) = res.Rp(3:4);
% cftool(res.t, res.Rp)

f = fit(times, res.Rp,'smoothingspline','SmoothingParam', 0.999830445339403);
Rp = f(times);

figure
hold on
plot(times, Rp)
plot(times, res.Rp)

res.Rp = Rp;

%%
load('nstxu_obj_config2016_6565.mat')

tree = 'efit01';
tag = '.RESULTS.AEQDSK:VLOOPMHD'; 
plotit = 1;

% signal = mds_fetch_signal(shot, tree, times, tag, plotit);


ROOT = getenv('RAMPROOT');
tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';
opts.cache_dir = [ROOT '/fetch/cache/'];
opts.plotit = 0;
eqs = fetch_eqs_nstxu(shot, times, tree, tokamak, server, opts);

t = eqs.time;
N = length(t);
dt = mean(diff(t));
psibry = [eqs.gdata(:).psibry];
psibrydot = gradient(psibry, dt);

hold on
plot(t, -psibrydot)
xlim([0.1 0.8])


Lp = nan(N,1);
Li = nan(N,1);
Le = nan(N,1);
li = nan(N,1);
psic = nan(N,1);
for i = 1:N
  try
    [Lp(i), Li(i), Le(i), li(i), ~, ~, psic(i)] = inductance(eqs.gdata(i), tok_data_struct);
  catch
  end
end

Ip = [eqs.gdata(:).cpasma]';
ic = [eqs.gdata(:).ic];


% close all
% figure
% hold on
% plot(t,Lp)
% plot(t, res.Rp)


Vb = -smooth(gradient(psibry,dt));
Vc = -smooth(gradient(psic,dt));
Vr = res.Rp .* Ip;

% figure
% hold on
% plot(t,smooth(Vr))
% plot(t,smooth(Vc))
%%
Vb = smooth(Vb);
Vr = smooth(Vr);
Vc = smooth(Vc);
%%
L = Li;
L(9:N) = 0;

for i = 2:N
  dldt = 0.9 * (Vr(i) - Vc(i)) / Ip(i);
  L(i) = L(i-1) + dt*dldt;
end

figure
hold on
plot(t,L)
plot(t,Li)
plot(t,Ip / median(Ip) * median(Li))
xlim([0 1])

%%
close all
I = Ip;
I(9:N) = 0;
for i = 2:N
  didt = (Vb(i) + Vc(i) - 1.7*Vr(i)) / Li(i);
  % didt = (Vb(i) - Vr(i)) / Lp(i);
  I(i) = I(i-1) + dt*didt;
end

figure
hold on
plot(t,I)
plot(t,Ip)

% figure
% hold 
% plot(t, Vr, t, Vc, t, Vb)

%%

struct_to_ws(tok_data_struct);

ip = [eqs.gdata(:).cpasma]';
ic = [eqs.gdata(:).ic];
ipdot = gradient(ip, dt);

psi_app = mpc*ic;
r = [0.4528    0.7560    1.3974    0.7327];
z = [0.0466    0.9913   -0.0816   -1.2012];
for i = 1:N
  tmp = reshape(psi_app(:,i), nz, nr);
  psib_app(i) = mean(bicubicHermite(rg, zg, tmp, r, z));
end

Le_approx = mean(Le);
psib_pla = Le_approx * ip;

psib = psib_app + psib_pla';

%%
L = Li;
L(9:N) = 0;

for i = 8:N-1
  dldt = 2 / ip(i) * ( -0.5*ipdot(i)*L(i) + 1*(Vb(i) - Vr(i)));
  L(i+1) = L(i) + dt*dldt;
end

figure
hold on
plot(t,L)
plot(t,Li)
plot(t,Ip / median(Ip) * median(Li))
xlim([0 1])








































