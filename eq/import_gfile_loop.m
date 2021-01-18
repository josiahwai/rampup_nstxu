clc; close all; clear all; clear gsdesign spec init config


shot = 204660;
saveit = 0;
plotit = 1;
% times = [60:10:120 140:10:300];
times = 60;

for i = 1:length(times)
  time_ms = times(i);
  disp(time_ms)
  try
    eq = import_gfile(shot, time_ms, saveit, plotit);
  catch
    warning(['Failed to import ' num2str(time_ms)])
  end
end



