clear; clc; close all

% ========
% SETTINGS
% ========
shot = 204660;
time_ms = 100;
shotdir = '/u/jwai/rampup_nstxu/eq/geqdsk/';
load('eq204660_100.mat')
load('nstxu_obj_config2016_6565.mat')


% ====================================
% Define spec,init,config for gsdesign
% ====================================

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





















