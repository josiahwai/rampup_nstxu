clear; clc; close all

% ========
% SETTINGS
% ========
shot = 204660;
time_ms = 240;
shotdir = '/u/jwai/rampup_nstxu/eq/geqdsk/';
load('nstxu_obj_config2016_6565.mat')

savedir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/';
saveAB = 1;

% ====================================
% Define spec,init,config for gsdesign
% ====================================

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
[vvgroup, vvcirc] = nstxu2020_vvcirc;
build_inputs.vvcirc = vvcirc;
build_inputs.vvgroup = vvgroup';
  
nstxu_sys = build_nstxu_system(build_inputs); 
delete('NSTXU_netlist.dat')

ccnames = {'OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', ...
  'PF5', 'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL'};

vvnames = {};
for i = 1:max(vvgroup)
  vvnames{i} = ['vv' num2str(i)];
end
 
if saveAB
  sys.A = nstxu_sys.amat;
  sys.As = removeFirstEigval(nstxu_sys.amat);
  sys.B = nstxu_sys.bmat;
  sys.inputs = ccnames;
  sys.states = [ccnames vvnames {'IP'}];  
  fn = [savedir num2str(shot) '_' num2str(time_ms) '_sys.mat'];
  save(fn, 'sys')
end

  


















