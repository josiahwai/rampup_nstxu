
function [Rp,t] = load_rp_profile(plotit)

if ~exist('plotit', 'var'), plotit = 0; end

ROOT = getenv('RAMPROOT');

% res_dir = [ROOT 'buildmodel/resistivity_calc/res/'];
% res_dir = [ROOT 'buildmodel/resistivity_calc/res2/'];
res_dir = [ROOT '/sysid/old/old_sysids/fit_plasma_resistance/fits/'];   % USED TO BE: res_dir = [ROOT 'sysid/plasma_resistance/fits/'];
% res_dir = [ROOT 'sysid/fit_plasma_resistance/fits/'];

d = dir([res_dir 'res*']);  

t_sample = .025:.01:0.8;

for i = 1:length(d) 
  res = load([d(i).folder '/' d(i).name]).res;  
  [~,isoutl] = rmoutliers(res.Rp, 'movmean', 15);       % remove shot outliers 
%   isoutl = false(size(isoutl));
  Rp_sample(i,:) = interp1(res.t(~isoutl), res.Rp(~isoutl), t_sample);  
end

for i = 1:length(t_sample)  
  [~, isoutl] = rmoutliers(Rp_sample(:,i));   % remove time outliers
%   isoutl = false(size(isoutl));
  Rp_sample(isoutl,i) = nan;
end

Rp = nanmedian(Rp_sample);
% Rp = nanmean(Rp_sample);

% extend endpoints for more robust fitting
i = ~isnan(Rp);
t_extend = linspace(t_sample(1) - 0.02, t_sample(end) + 0.02);
Rp = interp1(t_sample(i), Rp(i), t_extend, 'nearest', 'extrap');

% smoothing fit
% f = fit(t_extend(:), Rp(:),'smoothingspline','SmoothingParam', 1-1e-4);
f = fit(t_extend(:), Rp(:),'smoothingspline','SmoothingParam', 1-1e-5);

% extend to future times
t = linspace(0,2,200);
Rp = interp1(t_sample, f(t_sample), t, 'nearest', 'extrap');

if plotit 
  figure
  ax  = gca;
  co1 = ax.ColorOrder(1,:);
  co2 = ax.ColorOrder(2,:);
  hold on
  plot(t_sample, Rp_sample(1,:), 'color', co1)
  plot(t, Rp, 'color', co2, 'linewidth', 3)    
  plot(t_sample, Rp_sample, 'color', co1)
  plot(t, Rp, 'color', co2, 'linewidth', 3)  
  legend('NSTXU Experimental Fits', 'Average', 'fontsize', 16)  
  ylabel('Resistance [Ohm]', 'fontsize', 14)
  xlabel('Time [s]', 'fontsize', 14)  
  xlim([0.05 0.8])
  ylim([-0.2 2.8]*1e-5)
end
















