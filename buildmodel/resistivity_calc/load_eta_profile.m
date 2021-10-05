
function [eta,t] = load_eta_profile(plotit)

if ~exist('plotit', 'var'), plotit = 0; end

ROOT = getenv('RAMPROOT');
res_dir = [ROOT 'buildmodel/resistivity_calc/res/'];
d = dir([res_dir 'res*']);  % load all the resistances that were fit with fit_rp.m

t_sample = .025:.01:1;

for i = 1:length(d) 
  res = load([d(i).folder '/' d(i).name]).res;  
  [~,j] = rmoutliers(res.eta, 'movmean', 15);       % remove shot outliers 
  eta_sample(i,:) = interp1(res.t(~j), res.eta(~j), t_sample);  
end

for i = 1:length(t_sample)  
  [~, isoutl] = rmoutliers(eta_sample(:,i));   % remove time outliers
  eta_sample(isoutl,i) = nan;
end

eta = nanmedian(eta_sample);


% extend endpoints for more robust fitting
i = ~isnan(eta);
t_extend = linspace(t_sample(1) - 0.02, t_sample(end) + 0.02);
eta = interp1(t_sample(i), eta(i), t_extend, 'nearest', 'extrap');

% smoothing fit
f = fit(t_extend(:), eta(:),'smoothingspline','SmoothingParam', 1-1e-4);

% extend to future times
t = linspace(0,2,200);
eta = interp1(t_sample, f(t_sample), t, 'nearest', 'extrap');

if plotit 
  figure
  hold on
  plot(t_sample, eta_sample(1,:), 'color', 'b')
  plot(t, eta, 'r', 'linewidth', 3)    
  plot(t_sample, eta_sample, 'color', 'b')
  plot(t, eta, 'r', 'linewidth', 3)  
  legend('2015-2016 all shots', 'Average', 'fontsize', 16)  
  ylabel('Resistivity [Ohm/m^2]', 'fontsize', 14)
  xlabel('Time [s]', 'fontsize', 14)
  
end
















