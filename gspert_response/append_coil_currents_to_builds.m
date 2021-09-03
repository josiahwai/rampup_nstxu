% append coil currents

clear; clc; close all
saveit = 1;

builds_dir = '/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/raw_data/pertbuilds_val';
d = dir(builds_dir);
d(1:2) = [];

vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);


for i = 1:length(d)
  
  disp(i / length(d))
  
  targs = load([d(i).folder '/' d(i).name]).targs;
  
  coils = fetch_coilcurrents_nstxu(targs.shot, targs.actualtime);
  
  targs.coil_currents = coils.icx(circ.ikeep,:)';
  targs.vessel_currents = coils.ivx';
  
  if saveit
    save([d(i).folder '/' d(i).name], 'targs')
  end
end




















































