clc; close all; clear all; clear gsdesign spec init config

shot = 204660;
saveit = 1;
plotit = 1;
times = 60:10:300;
shotdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk/';
savedir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk_import/';


for i = 1:length(times)
  time_ms = times(i);
  disp(time_ms)
  try
    eq = import_geqdsk(shot, time_ms, shotdir, savedir, saveit, plotit);
  catch
    warning(['Failed to import ' num2str(time_ms)])
  end
end



