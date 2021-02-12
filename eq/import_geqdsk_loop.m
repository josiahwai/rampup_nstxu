clc; close all; clear all; clear gsdesign spec init config

shot = 204660;
saveit = 1;
plotit = 1;
times = 20:10:1000;
shotdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk/';
savedir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk_import/';

coils = mds_fetch_current_voltage(shot, 0);

load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);
iused = circ.ikeep;

pause('off')

for i = 1:length(times)
  time_ms = times(i);
  disp(time_ms)
  try    
    [~,k] = min(abs(coils.t - time_ms/1000));    
    ic(iused) = coils.ic(:,k);    
    iv = coils.iv(:,k); 
    import_geqdsk(shot, ic(:), iv(:), time_ms, shotdir, savedir, saveit, plotit);
  catch
    warning(['Failed to import ' num2str(time_ms)])
  end
end



