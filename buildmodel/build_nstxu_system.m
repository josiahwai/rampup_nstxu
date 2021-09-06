clc; close all; clear

% ========
% SETTINGS
% ========
shot = 204660;
times = 30:10:960;
saveit = 0;
savedir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/mcc_coneqt/';
eqdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk_import/';
use_coneqt_tok_data_struct = 1;

% ===========
% Build loop
% ===========
if use_coneqt_tok_data_struct
  new = load('coneqt_nstxu_obj_config2016_6565.mat');
  old = load('nstxu_obj_config2016_6565.mat');
  tok_data_struct = copyfields(new.tok_data_struct, old.tok_data_struct, {}, 0);
  Rvv_fit = load('Rvv_fit.mat').Rvv_fit;
else
  tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
end
circ = nstxu2016_circ(tok_data_struct); % circuit connections
fit_Rp = load('fit_Rp5.mat').fit_Rp;  % custom plasma resistance

% Build it
for itime = 1:length(times)
  
  time_ms = times(itime)
  
  % load equilibrium       
  load([eqdir 'eq' num2str(shot) '_' num2str(time_ms) '.mat']);
  
  % define build_inputs
  build_inputs.tokamak = 'NSTXU';
  build_inputs.vacuum_objs = tok_data_struct;
  build_inputs.ichooseq = 4;
  build_inputs.equil_data = eq;
  build_inputs.irzresp_dynamic = 5;
  build_inputs.irzresp_output = 5;
  build_inputs.iplcirc = 1;
  build_inputs.cccirc = circ.cccirc(:);
  build_inputs.vvcirc = circ.vvcirc(:);
  build_inputs.vvgroup = circ.vvgroup(:);
  build_inputs.Rp = fit_Rp(time_ms/1000);    
  
  nstxu_sys = build_tokamak_system(build_inputs); 
  delete('NSTXU_netlist.dat')
    
  if use_coneqt_tok_data_struct
    Pxx = nstxu_sys.Pxx;
    Pcc = nstxu_sys.Pcc;
    
    nstxu_sys.lstar = Pxx' * nstxu_sys.mxx * Pxx; 
    
    lstari = inv(nstxu_sys.lstar);
    rxx = [diag(Pcc' * diag(tok_data_struct.resc) * Pcc); Rvv_fit; nstxu_sys.Resp];
    A = -lstari * diag(rxx);
    B = lstari(:,1:circ.ncx);
    
    
  else
    Pxx = nstxu_sys.Pxx;
    rxx = Pxx' * nstxu_sys.rxx * Pxx;
    nstxu_sys.lstar = Pxx' * nstxu_sys.mxx * Pxx;
    lstari = inv(nstxu_sys.lstar);
    A = -lstari * rxx;
    B = lstari(:,1:circ.ncx);
    
%     A = nstxu_sys.amat;
%     B = nstxu_sys.bmat;
  end
  
  % Some coils turned off in campaign, remove from circuit equations  
  A(circ.iremove,:) = 0;
  A(:,circ.iremove) = 0;
  B(circ.iremove,:) = 0;
  B(:,circ.iremove) = 0;
      
  if saveit
    sys.Lp = nstxu_sys.Lp;
    sys.A = A;
    sys.As = removeFirstEigval(A);                
    sys.B = B;
    sys.lstar = nstxu_sys.lstar;
    sys.rxx = diag(nstxu_sys.Pxx' * nstxu_sys.rxx * nstxu_sys.Pxx);
    sys.inputs = circ.ccnames;
    sys.states = [circ.ccnames circ.vvnames 'Ip'];
    fn = [savedir num2str(shot) '_' num2str(time_ms) '_sys.mat'];
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    save(fn, 'sys')
  end
end















