function [Lp, times] = read_Lp()

shot = 204660;
modeldir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/mcc/';
times = [30:10:960] / 1e3;
Lp = zeros(length(times), 1);
for i = 1:length(times)
  time_ms = times(i) * 1e3;
  sys = load([modeldir num2str(shot) '_' num2str(time_ms) '_sys.mat']).sys;      
  Lp(i) = sys.Lp;
end