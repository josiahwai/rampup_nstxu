ccc

save_dir = '/Users/jwai/Research/rampup_nstxu/gspert_response/builds_indiv/';
mode = 'val';

tree = 'EFIT01';
tokamak = 'nstxu';
server = 'skylark.pppl.gov:8501';

shots_times = load('/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/shots_times/requested/val_shots_times.mat').shots_times;

shotlist = shots_times.shots;
timelist = shots_times.times;
targ_shots = shotlist;

for ishot = 1:length(targ_shots)
  
  shot = targ_shots(ishot);
  
  fprintf('\n\n')
  disp(['Fetching shot ' num2str(ishot) ' of ' num2str(length(targ_shots))])
  fprintf('\n\n')
  
  requesttimes = timelist{shotlist==shot};
  
  tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
  
  circ = nstxu2016_circ(tok_data_struct);
  
  build_inputs.tokamak = 'NSTXU';
  build_inputs.vacuum_objs = tok_data_struct;
  build_inputs.ichooseq = 4;
  build_inputs.irzresp_dynamic = 5;
  build_inputs.irzresp_output = 5;
  build_inputs.iplcirc = 1;
  build_inputs.cccirc = circ.cccirc(:);
  build_inputs.vvcirc = circ.vvcirc(:);
  build_inputs.vvgroup = circ.vvgroup(:);
  
  eqs = read_eq(shot, 'all', tree, tokamak, server);
  
  [~, t_idx] = min(abs(requesttimes - eqs.time));
  actualtimes = eqs.time(t_idx);
  
  targs = struct;
  isample = 0;
  ibad = 0;
  
  for i = 1:length(t_idx)
    
    requesttime = requesttimes(i);
    actualtime = actualtimes(i);
    
    time = actualtime;
    disp(['time_ms: ' num2str(time*1000)])
    
    try
      eq = eqs.gdata(t_idx(i));
      eq.ecturn = tok_data_struct.ecnturn;
      eq.ecid   = ones(size(eq.ecturn));
      eq.turnfc = tok_data_struct.fcnturn';
      eq.fcturn = circ.fcfrac;
      eq.fcid = circ.fccirc';
      
      build_inputs.equil_data = eq;
      
      sys = build_tokamak_system(build_inputs);
      delete('NSTXU_netlist.dat')
      
      P = circ.Pxx(1:end-1,1:end-1);
      xmatx = sys.xmatx(1:end-1,1:end-1);
      
      xmat = P'*xmatx*P;
      xmat = double(xmat);
      
      e = esort(eig(sys.amat(1:end-1,1:end-1)));
      gamma = real(e(1));
      gamma = double(gamma);
      disp(['gamma: ' num2str(gamma)])
      
      dcphidis = [sys.gspert_data.dcphidis sys.gspert_data.dcphidip(:)];
      dcphidix = dcphidis * circ.Pxx;
      dcphidix = double(dcphidix);
      
      % save to struct
      isample = isample+1;
      
      targs.xmat(isample,:,:) = xmat;
      targs.gamma(isample) = gamma;
      targs.dcphidix(isample,:,:) = dcphidix(:,circ.ikeep);
      targs.requesttime(isample) = requesttime;
      targs.actualtime(isample) = actualtime;
      targs.ip(isample) = eq.cpasma;
      targs.pprime(isample,:) = eq.pprime;
      targs.ffprim(isample,:) = eq.ffprim;
      targs.pres(isample,:) = eq.pres;
      targs.psirz(isample,:,:) = eq.psirz;
      targs.pcurrt(isample,:,:) = eq.pcurrt;
      targs.shot = shot;
      targs.rcur(isample) = sum(tok_data_struct.rgg(:).*eq.pcurrt(:)) / eq.cpasma;
      targs.zcur(isample) = sum(tok_data_struct.zgg(:).*eq.pcurrt(:)) / eq.cpasma;

      
    catch
      warning('Bad shot!')
      ibad = ibad+1;
      bad.shot(ibad) = shot;
      bad.time(ibad) = time;
    end
  end
  
  fn = [save_dir mode '_response_' num2str(shot) '.mat'];
  save(fn, 'targs', '-v7.3')
end
disp('Done!')



















































