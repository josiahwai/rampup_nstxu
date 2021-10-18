function p = profile_test2(np, nf, plotit, trange)

load('/Users/jwai/Research/nstxu-nns/data/rawdata/data_by_var/train/pprime.mat')
load('/Users/jwai/Research/nstxu-nns/data/rawdata/data_by_var/train/ffprim.mat')
load('/Users/jwai/Research/nstxu-nns/data/rawdata/data_by_var/train/time.mat')

if ~exist('plotit', 'var'), plotit = 0; end
if ~exist('trange', 'var'), trange = [min(time) max(time)]; end

% restrict to trange
i = time > trange(1) & time < trange(2);
ffprim = ffprim(i,:);
pprime = pprime(i,:);

% remove outliers
[~,j] = rmoutliers(ffprim(:,1), 'ThresholdFactor', 5);
[~,k] = rmoutliers(ffprim(:,end), 'ThresholdFactor', 5);
i = ~(j|k);
ffprim = ffprim(i,:);


[score, coeff, ~, ~, explained, mu] = pca(ffprim, 'NumComponents', nf, 'Centered', 'on');
p.ffprim = score * diag(mean(abs(coeff)) .* mode(sign(coeff)));
p.ffprim_mean = mu;

% explained(1:5)

[score, coeff] = pca(pprime, 'NumComponents', np, 'Centered', 'off');
p.pprime = score * diag(mean(abs(coeff)) .* mode(sign(coeff)));


if plotit
  psin = linspace(0,1,length(p.ffprim_mean));
  
  figure
  hold on
  plot(psin, p.pprime)
  yyaxis right
  hold on  
  plot(psin, p.ffprim)
  plot(psin, p.ffprim_mean, '--')
  legend('PPRIME', 'FFPRIME', 'fontsize', 18)
end











