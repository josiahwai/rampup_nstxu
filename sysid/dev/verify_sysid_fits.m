clear all; clc; close all

ROOT = getenv('RAMPROOT');
enforce_stability = 0;

shot = 204660;  % We will try to recreate this shot

% t0 = 0.07;
% tf = 0.9;
% N = 50;
% t = linspace(t0, tf, N);
% ts = mean(diff(t));
% 
% % sysid fit for plasma resistance
% res = load(['res' num2str(shot) '.mat']).res;
% Rp = interp1(res.t, res.Rp, t);
% 
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);
% 
% tree = 'EFIT01';
% tokamak = 'nstxu';
% server = 'skylark.pppl.gov:8501';
% opts.cache_dir = [ROOT '/fetch/cache/'];
% efit01_eqs = fetch_eqs_nstxu(shot, t, tree, tokamak, server, opts);
% 
% % remove some zero-padding, oft causes trouble
% for i = 1:N 
%   k = efit01_eqs.gdata(i).rbbbs==0 & efit01_eqs.gdata(i).zbbbs==0;
%   efit01_eqs.gdata(i).rbbbs(k) = [];
%   efit01_eqs.gdata(i).zbbbs(k) = [];  
% end
% 
% 
% % calculate inductances
% for i = 1:N
%   [Lp{i}, Li{i}, Le{i}, li{i}, mcIp{i}, mvIp{i}] = inductance(efit01_eqs.gdata(i), tok_data_struct);
% end
% 
% include_coils = {'OH', 'PF1aU', 'PF1bU', 'PF1cU', 'PF2U', 'PF3U', 'PF4', ...
%   'PF5', 'PF3L', 'PF2L', 'PF1cL', 'PF1bL', 'PF1aL'};
% 
% vobjcsignals = get_vobjcsignals(shot, [], [], include_coils);
% v = smoothdata(vobjcsignals.sigs, 'movmean', 50);
% v = interp1(vobjcsignals.times, v, t);



%%
% Pxx = circ.Pxx(1:end-1,1:end-1); % remove Ip circuit from transition map
% Rxx = Pxx' * diag([tok_data_struct.resc; tok_data_struct.resv]) * Pxx;
% Mxx = Pxx' * [tok_data_struct.mcc tok_data_struct.mcv; tok_data_struct.mcv' tok_data_struct.mvv] * Pxx;
% 
% ext_fit = load([ROOT '/sysid/fit_Rext_Lext/ext_fit.mat']).ext_fit;
% Mxx_ext = diag([ext_fit.Lext; zeros(circ.nvx,1)]);
% Rxx_ext = diag([ext_fit.Rext; zeros(circ.nvx,1)]);
% 
% Mxx = Mxx + Mxx_ext;
% Rxx = Rxx + Rxx_ext;

s = load('shotdata.mat').shotdata.s204660;
s = struct_fields_to_double(s);

% remove nans
idxnan = any(isnan(s.ic')) | any(isnan(s.icdot')) | any(isnan(s.mcIp'));
fns = fieldnames(s);
for j = 1:length(fns)
  x = s.(fns{j});
  if size(x,1) == 1, x = x'; end
  s.(fns{j}) = x(~idxnan,:);
end

ts = mean(diff(s.tsample));
t = s.tsample;
N = length(s.tsample);

load('sysid_fits.mat')
Mxx = sysid_fits.Mxx;
Rxx = sysid_fits.Rxx;
rc = Rxx(circ.iicx);
mcv = Mxx(circ.iicx, circ.iivx);
mcc = Mxx(circ.iicx, circ.iicx);
mvv = Mxx(circ.iivx, circ.iivx);
Mxxi = inv(Mxx);

x0 = [s.ic(1,:) s.iv(1,:)]';
x = x0;

for i = 1:N
  
  A = -Mxxi * diag(Rxx);
  B = Mxxi(:,circ.iicx);
  B(:,end+1) = -[s.mcIp(i,:) s.mvIp(i,:)]';
  
  u = [s.v(i,:) s.ipdot(i)]';
  
  [Ad,Bd] = c2d(A,B,ts);
  
  x = Ad*x + Bd*u;
  
  % enforce nonnegative current in unipolar coils
  k = x(circ.ii_unipolar) < 0;
  x(circ.ii_unipolar(k)) = 0;
    
  xall(i,:) = x;
      
end

figure
hold on
plot(t, s.ic, '--r')
plot(t, xall(:,circ.iicx), '--b')


figure
hold on
plot(t, s.iv, '--r')
plot(t, xall(:,circ.iivx), '--b')


% figure
% hold on
% plot(t, s.ip, '--r')
% plot(t, xall(:,circ.iipx), '--b')












































