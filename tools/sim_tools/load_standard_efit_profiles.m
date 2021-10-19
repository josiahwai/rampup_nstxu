
function [pprime_profile, ffprim_profile] = load_standard_efit_profiles(plotit, trange)

load('/Users/jwai/Research/nstxu-nns/data/rawdata/data_by_var/train/pprime.mat')
load('/Users/jwai/Research/nstxu-nns/data/rawdata/data_by_var/train/ffprim.mat')
load('/Users/jwai/Research/nstxu-nns/data/rawdata/data_by_var/train/time.mat')

if ~exist('plotit', 'var'), plotit = 0; end
if ~exist('trange', 'var'), trange = [min(time) max(time)]; end

i = time > trange(1) & time < trange(2);
ffprim = ffprim(i,:);
pprime = pprime(i,:);


[~,j] = rmoutliers(ffprim(:,1), 'ThresholdFactor', 5);
[~,k] = rmoutliers(ffprim(:,end), 'ThresholdFactor', 5);
i = ~(j|k);
ffprim = ffprim(i,:);


% the first component is v(:,1), the other terms here are just to get
% the sign right and a good average scaling

[u,s,v] = svd(pprime, 'econ');
pprime_profile = mode(sign(u(:,1))) * mean(abs(u(:,1))) * s(1) * v(:,1);  


[u,s,v] = svd(ffprim, 'econ');
ffprim_profile = mode(sign(u(:,1))) * mean(abs(u(:,1))) * s(1) * v(:,1);

if plotit
  psin = linspace(0,1,length(ffprim_profile));
  
  figure
  plot(psin, pprime_profile)
  yyaxis right
  plot(psin, ffprim_profile)
  legend('PPRIME', 'FFPRIME', 'fontsize', 18)
  
%   figure
%   plot(psin, ffprim)
end











