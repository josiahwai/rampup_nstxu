clear all; clc; close all
ROOT = getenv('RAMPROOT');

shot = 204660;  

load('nstxu_obj_config2016_6565.mat')
struct_to_ws(tok_data_struct);
circ = nstxu2016_circ(tok_data_struct);
mcc = circ.Pcc' * mcc * circ.Pcc;

fetch_times = [];
tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';
opts.cache_dir = [ROOT '/fetch/cache/'];
efit01_eqs = fetch_eqs_nstxu(shot, fetch_times, tree, tokamak, server, opts);

eq  = efit01_eqs.gdata(20);
t = double(efit01_eqs.time);

%% weights for ip tracking


tmp = [efit01_eqs.gdata(:).ic];
ioh = tmp(1,:);
tefit = double(efit01_eqs.time);

lc = mcc(1,1);
[lp, ~, ~, ~, mcIp] = inductance(eq, tok_data_struct);
mcp = mcIp(1);
M = [lc mcp; mcp lp];
invM = inv(M);

rp = 1e-6;
rc = resc(1);
R = diag([rc rp]);

A = -invM*R;
B = invM(:,1);

t_ref = [0.0018 0.0682 0.2194 0.4000 0.9604];
ip_ref = [0.0328 2.5292 5.4460 7.7058 7.6270] * 1e5;

x0 = [efit01_eqs.gdata(1).ic(1) efit01_eqs.gdata(1).cpasma]';
x = x0;
N = 100;
t = linspace(0,1,N);
dt = mean(diff(t));

% v0 = rc * x(1);
v0 = 0;
dv = 0;
v = v0 + dv;
int_e = 0;

iptarg = interp1(t_ref, ip_ref, t, 'linear', 'extrap');

kp = 20e-3;
ki = 2e-1; %8e-2;

for i = 1:N
  xdot = A*x + B*v;
  x = x + xdot*dt;
 
  ic = x(1);
  ip = x(2);
  
  e = iptarg(i) - ip;
  int_e = int_e + dt*e;
  
  dv = -kp*e - ki*int_e;
%   v0 = rc*ic;
  v0 = 0;
  v = v0+dv;
  
  xall(i,:) = x;
  vall(i,:) = v;
  
end

ic = xall(:,1);
ip = xall(:,2);

figure
subplot(311)
plot(t, iptarg, '--', t, ip)
subplot(312)
plot(tefit, ioh, '--', t, ic)
subplot(313)
plot(t, vall)


%% weights for general current tracking
close all

lc = diag(mcc);
rc = circ.Pcc \ (resc .* circ.ccfrac');
ic = [efit01_eqs.gdata(:).icx];
t = double(efit01_eqs.time);
N = length(t);
dt = mean(diff(t));

kp(1:13) = 5e-2;
ki(1:13) = 8e-1;
kp(8) = 3e-1;
ki(8) = 8e-1;
kp = kp(:);
ki = ki(:);


x0 = ic(:,1);
x = x0;

v0 = rc.*x;
dv = 0;
v = v0+dv;

e_int = 0;
xall =  []; vall = []; v0all = [];

for i = 1:N
  
  xdot = -rc ./ lc .* x + 1./lc .* v;
  x = x + xdot*dt;
  
  e = ic(:,i) - x;
  e_int = e_int + e*dt;
  
  dv = kp.*e + ki.*e_int;
  v0 = rc.*x;
  v = v0 + dv;
  
  xall(:,i) = x;
  vall(:,i) = v;
  v0all(:,i) = v0;
end

icoil = 2:13;

figure
subplot(211)
plot(t, ic(icoil,:), '--r', t, xall(icoil,:), '-b')
subplot(212)
hold on
plot(t, vall(icoil,:), 'b')
plot(t, v0all(icoil,:), 'g')
  


%%
% plot_eq(eq)

close all
mpcx = mpc * circ.Pcc;


% % PF5
icoil = 8;
P.r = 1.45;
P.z = 0;
scale = bicubicHermite(rg, zg, reshape(mpcx(:,icoil), nz, nr), P.r, P.z);
kp = 3e-1 / scale;
ki = 0;


% PF1AU/L
% icoil = 2;
% P.z = 0.6;
% P.z = 1.15;
% scale = bicubicHermite(rg, zg, reshape(mpcx(:,icoil), nz, nr), P.r, P.z);
% kp = 2e-2 / scale;
% ki = 0;

% PF2U/L
% icoil = 5;
% P.r = 0.9;
% P.z = 1;
% scale = bicubicHermite(rg, zg, reshape(mpcx(:,icoil), nz, nr), P.r, P.z);
% kp = 2e-2 / scale;
% ki = 0;

% PF3U/L
% icoil = 6;
% P.r = 1.03;
% P.z = 0.9;
% scale = bicubicHermite(rg, zg, reshape(mpcx(:,icoil), nz, nr), P.r, P.z);
% kp = 4e-2 / scale;
% ki = 0;


[P.psi, P.psi_r, P.psi_z] = bicubicHermite(rg, zg, eq.psizr, P.r, P.z);
P.dpsi = P.psi - eq.psibry;


lc = diag(mcc);
rc = circ.Pcc \ (resc .* circ.ccfrac');
rc = rc(icoil);
lc = lc(icoil);

eq = efit01_eqs.gdata(100);
eq.psizr_pla = eq.psizr - reshape(mpc*eq.ic, nz, nr);

x = eq.ic(icoil);
N = 100;
t = linspace(0,1,N);
dt = mean(diff(t));

dv = 0;
v0 = rc*x;
v = v0 + dv;

e_int = 0;

xall=[]; vall=[]; eall=[];

for i = 1:N
  
  xdot = -rc / lc * x + 1/lc * v;
  x = x + xdot*dt;
  
  icx = eq.icx;
  icx(icoil) = x;
  psizr_app = reshape(mpcx*icx, nz, nr);
  psizr = psizr_app + eq.psizr_pla;
  
  psi = bicubicHermite(rg, zg, psizr , P.r, P.z);
  y = psi - eq.psibry;
  r = 0;
  
  e = r-y;
  e_int = e_int + e*dt;
  
  dv = kp.*e + ki.*e_int;
  v0 = rc.*x;
  v = v0 + dv;
  
  xall(:,i) = x;
  vall(:,i) = v;
  eall(:,i) = e;
  psiall(:,i) = psi;
end

figure
subplot(311)
plot(t, eall)
subplot(312)
plot(t,xall)
subplot(313)
hold on
plot(t, vall)


% figure
% ic = [efit01_eqs.gdata(:).icx];
% t = efit01_eqs.time;
% plot(t,ic(8,:))




















