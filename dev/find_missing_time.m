clear; clc; close all


% A = load('/Users/jwai/Research/rampup_nstxu/gspert_response/target_xmats/val_responses_A.mat').response_data;
% B = load('/Users/jwai/Research/rampup_nstxu/gspert_response/target_xmats/val_responses_B.mat').response_data;
% shots_times = load('/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/shots_times/actual/val_shots_times.mat').shots_times;


A = load('/Users/jwai/Research/rampup_nstxu/gspert_response/target_xmats/train_responses_A.mat').response_data;
B = load('/Users/jwai/Research/rampup_nstxu/gspert_response/target_xmats/train_responses_B.mat').response_data;
shots_times = load('/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/shots_times/actual/train_shots_times.mat').shots_times;
shots = shots_times.shots;
times = shots_times.times;

shotsAB = [A.validshots B.validshots];
timesAB = [A.validtimes B.validtimes];

uniq_shots = unique(shots);

for shot = uniq_shots
  
  n = sum(shots==shot);
  m = sum(shotsAB==shot);
  
  if m ~= n
    disp(shot)
  end  
end

shot = 204082;
t = times(shots==shot);
tab = timesAB(shotsAB==shot);
shot = 204082;
time = 1.4;
time = 1.45;





















