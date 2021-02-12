% sim settings
shot = 204660;
times = [20:10:1000];
saveit = 0;
pause('off')
savedir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/lstar/';
eqdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk_import/';
load('nstxu_obj_config2016_6565.mat')

% Build Loop
times_failed = [];
for i = 1:length(times)
  t = times(i);  
  try
    load([eqdir 'eq' num2str(shot) '_' num2str(t) '.mat']);    
    build_inputs.tokamak = 'NSTXU';
    build_inputs.vacuum_objs = tok_data_struct;
    build_inputs.ichooseq = 4;
    build_inputs.equil_data = eq;
    build_inputs.irzresp_dynamic = 5;
    build_inputs.irzresp_output = 5;
    build_inputs.iplcirc = 1;
    circ = nstxu2016_circ(tok_data_struct);
    build_inputs.cccirc = circ.cccirc;
    build_inputs.vvcirc = circ.vvcirc;
    
    nstxu_sys = build_tokamak_system(build_inputs);
    delete('NSTXU_netlist.dat')
    lstar = nstxu_sys.lstar;
    
    if saveit
      fn = [savedir num2str(shot) '_' num2str(time_ms) '_lstar.mat'];
      save(fn, 'lstar');
    end
    
  catch
    times_failed = [times_failed; t];
  end
end

disp('Did not import for times: ')
disp(times_failed)

% Failed for: 20 970 980 990 1000






