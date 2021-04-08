clear all; clc; close all;

addpath(genpath('/Users/jwai/Research/rampup_nstxu'));

shots_times =load('/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/shots_times/actual/train_shots_times.mat').shots_times;
xmats = load('/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/raw_data/train_xmat.mat').xmat;

vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);

mcc = tok_data_struct.mcc;
mcv = tok_data_struct.mcv;
mvv = tok_data_struct.mvv;
M = [mcc mcv; mcv' mvv];

mxx = circ.Pxx(1:end-1,1:end-1)' * M * circ.Pxx(1:end-1,1:end-1);

shots = shots_times.shots;
uniq_shots = unique(shots);

%%
close all
figure
hold on

for ishot = 1:30 %length(uniq_shots)
  
  shot = uniq_shots(ishot);
  
  idx = find(shots==shot);
  t{ishot} = shots_times.times(idx);

  gamma{ishot} = [];
  for i = idx
    xmat = squeeze(xmats(i,:,1:end-1));
    A = -inv(xmat + mxx)*vac_sys.Rxx;
    eigs = esort(eig(A));
    e = max(real(eigs(1)), 0);
    gamma{ishot} = [gamma{ishot}; e];
  end

  plot(t{ishot}, gamma{ishot})

end




























