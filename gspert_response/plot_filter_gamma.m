clear; clc; close all

thresh = 20;
saveit = 0;
plotit = 1;
builds_dir = '/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/raw_data/pertbuilds_train';

d = dir(builds_dir);
d(1:2) = [];
nbadshots = 0;

for i = 1:length(d)
  
  disp(i / length(d))
  
  targs = load([d(i).folder '/' d(i).name]).targs;
    
  gamma_smooth = smoothdata(targs.gamma, 'movmedian', 20);       
  e = targs.gamma - gamma_smooth;  
  ibad = abs(e) > thresh;  
  igood = find(~ibad);
  ibad = find(ibad);
  
  %   if length(targs.ibad) / length(targs.gamma) > 0.3
  %     disp(targs.shot)
  %     nbadshots = nbadshots+1;
  %   end
  
  if plotit
    figure
    subplot(211)
    hold on
    title(num2str(targs.shot))
    scatter(igood, targs.gamma(igood))
    plot(e)
    scatter(ibad, targs.gamma(ibad), 'filled', 'r')
    ylim([0 200])
    subplot(212)
    scatter(igood, targs.gamma(igood))
  end    
  
  if saveit
    targs.igood = igood;
    targs.ibad = ibad;
    save([d(i).folder '/' d(i).name], 'targs')
  end
  
end
  
 
