clear; clc; close all;

saveit = 1;

builds_dir = '/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/raw_data/pertbuilds_train';

d = dir(builds_dir);
d(1:2) = [];

icoil = 2;


for i = 1:length(d)
  
  disp(i / length(d))
  
  fn = [d(i).folder '/' d(i).name];
  targs = load(fn).targs;
  
  % filter using gamma
  t = targs.actualtime;
  thresh = 20;
  gamma_smooth = smoothdata(targs.gamma, 'movmedian', 5);       
  e = targs.gamma - gamma_smooth;  
  ibad = abs(e) > thresh;  
  s = vecnorm(squeeze(targs.dpsidix(:,:,icoil)), 2, 2);
  ibad = ibad' | (s > 10*median(s));
  igood = find(~ibad);
  ibad = find(ibad);

  % smooth response
  y = squeeze(targs.dpsidix(:,:,icoil));
  y = interp1(t(igood), y(igood,:), t, 'previous', 'extrap');
  dpsidix_smooth = smoothdata(y, 1, 'movmedian', 13);
  dpsidix_smooth = smoothdata(dpsidix_smooth, 1, 'gaussian', 3);
  
%   figure
%   plot(dpsidix_smooth)
  
  tag = ['dpsidix' num2str(icoil-1) '_smooth'];
  targs.(tag) = dpsidix_smooth;
  
  if isfield(targs, 'dpsidix0_smooth_win_13')    
    targs.dpsidix0_smooth = targs.dpsidix0_smooth_win13;
    targs = rmfield(targs, 'dpsidix0_smooth_win13');
  end
  
  if saveit
    save(fn, 'targs')
  end
end












































