clear all; clc; close all

% ========
% SETTINGS
% ========
shot = 204660;
times = 60:10:300;
shotdir = '/u/jwai/rampup_nstxu/eq/geqdsk/';
load('nstxu_obj_config2016_6565.mat')
savedir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/modify_resv/';
saveit = 1;

% DEBGUG TEST: Modify resistance in vessel circuit 3 and 11
% [vvgroup, vvcirc, vvnames] = nstxu2020_vvcirc;
% i = find(vvgroup==3);
% j = find(vvgroup==11);
% tok_data_struct.resv(i) = tok_data_struct.resv(i) * 2;
% tok_data_struct.resv(j) = tok_data_struct.resv(j) / 3;


% coils were turned off for previous campaign, remove from circuit eqns
remove_coils = {'PF1BU', 'PF1BL', 'PF1CU', 'PF1CL', 'PF4'};  

% ===========
% Build loop
% ===========
for itime = 1:length(times)
  
  time_ms = times(itime)
  
  load(['eq' num2str(shot) '_' num2str(time_ms) '.mat'])

  % define build_inputs
  build_inputs.tokamak = 'NSTXU';
  build_inputs.vacuum_objs = tok_data_struct;
  build_inputs.ichooseq = 4;
  build_inputs.equil_data = eq;
  build_inputs.irzresp_dynamic = 5;
  build_inputs.irzresp_output = 5;
  build_inputs.iplcirc = 1;
  build_inputs.cccirc = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 ...
      11 12 13];
  [vvgroup, vvcirc, vvnames] = nstxu2020_vvcirc;
  build_inputs.vvcirc = vvcirc;
  build_inputs.vvgroup = vvgroup';
  
  % plasma resistance is much larger when limited
  islimited = time_ms < 400;
  if islimited
    if isfield(build_inputs, 'Te_res')
      build_inputs = rmfield(build_inputs, 'Te_res');
    end
    build_inputs.Rp = 2.44e-8 * 130 * interp1([60 70 80 100 150 190 200 210 220 300], [2.8 2.2 1.3 0.7 0.55 0.55 1.5 1.7 0.15 0.1], time_ms);
  else
    if isfield(build_inputs, 'Rp')
      build_inputs = rmfield(build_inputs, 'Rp');
    end
    build_inputs.Te_res = min( time_ms/0.3, 4000);
  end
%   load('fit_Rp3.mat')
%   build_inputs.Rp = fit_Rp(time_ms/1000) * 1.2;


  nstxu_sys = build_tokamak_system(build_inputs); 
  delete('NSTXU_netlist.dat')

  ccnames = {'OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', ...
    'PF5', 'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL'};
  
  % Remove unconnected coils from the circuit equations
  if exist('remove_coils', 'var') && ~isempty('remove_coils')
    iremove = zeros(size(ccnames));
    for i = 1:length(remove_coils)
      iremove = iremove | strcmp(ccnames, remove_coils{i});
    end
    nstxu_sys.amat(iremove,:) = [];
    nstxu_sys.amat(:,iremove) = [];
    nstxu_sys.bmat(:,iremove) = [];
    nstxu_sys.bmat(iremove,:) = [];
    nstxu_sys.lstar(iremove,:) = [];
    nstxu_sys.lstar(:,iremove) = [];        
    ccnames = {ccnames{~iremove}};
  end        
  
  % [A,B] = c2d(nstxu_sys.amat, nstxu_sys.bmat, .01);
  % A(end,end)
  
  if saveit
    sys.A = nstxu_sys.amat;
    sys.As = removeFirstEigval(nstxu_sys.amat);
    sys.B = nstxu_sys.bmat;
    sys.inputs = ccnames;
    sys.states = [ccnames cellstr(vvnames)' {'IP'}];  
    sys.lstar = nstxu_sys.lstar;
    fn = [savedir num2str(shot) '_' num2str(time_ms) '_sys.mat'];
    save(fn, 'sys')
  end
end















