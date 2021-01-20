clear all; clc; close all

% ========
% SETTINGS
% ========
shot = 204660;
% times = [60:10:120 140:10:300];
times = 200;
shotdir = '/u/jwai/rampup_nstxu/eq/geqdsk/';
load('nstxu_obj_config2016_6565.mat')
savedir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/geom1/';
eqdir = '/Users/jwai/Research/rampup_nstxu/eq/eq1/';
saveit = 0;

% some coils turned off for previous campaign, remove from circuit eqns
remove_coils = {'PF1BU', 'PF1BL', 'PF1CU', 'PF1CL', 'PF4'};  

% ===========
% Build loop
% ===========
for itime = 1:length(times)
  
  time_ms = times(itime)
  
  load([eqdir 'eq' num2str(shot) '_' num2str(time_ms) '.mat'])
 
  % Turn troubles
  fcnturn = tok_data_struct.fcnturn';
  Ifrac = zeros(size(fcnturn));
  cccirc = [2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 11 12 13] - 1;
  for icirc = 1:max(cccirc)
    icoils = find(cccirc == icirc);
    Ifrac(icoils) = fcnturn(icoils) / sum(fcnturn(icoils));
  end
  eq.turnfc = tok_data_struct.fcnturn';
  eq.fcturn = Ifrac;
  eq.fcid = cccirc;
  eq.ecid = [1 1 1 1 1 1 1 1];
  eq.ecturn = [112 110 109.5 108.5 108.5 109.5 110 112];
  
  
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
  
  % custom plasma resistance
  res.Rp = 3.17e-6 * [2.8 2.2 1.3 0.7 0.55 0.55 1.5 1.7 0.15 0.1];
  res.time = [60 70 80 100 150 190 200 210 220 300];
  build_inputs.Rp = interp1(res.time, res.Rp, time_ms);

  nstxu_sys = build_tokamak_system(build_inputs); 
  delete('NSTXU_netlist.dat')

  % Now we do some modification from the GA codes in order to include r and
  % z of current centroid as state variables in the model equations. 
  % state:=[Icx; Ivx; Ip; rcur; zcur];
%   [rzrig_data] = rzrig(eq,'NSTXU',tok_data_struct, 0, 0, [], 0, 0);
  
  
%   lstar = sys.lstar;
%   rxx = sys.Pxx' * sys.rxx * sys.Pxx;
%   drdix = sys.Pxx \ [sys.gspert_data.drdis sys.gspert_data.drdip]';
%   dzdix = sys.Pxx \ [sys.gspert_data.dzdis sys.gspert_data.dzdip]';
%   
%   Ip0 = eq.cpasma;
%   lstar(:,end+1) = 1./drdix;
%   lstar(:,end+1) = 1./dzdix;
%   lstar(end+1,:) = [drdix' -1 0];
%   lstar(end+1,:) = [dzdix' 0 -1];
%   rxx = blkdiag(rxx, zeros(2));
%   
%   A = -inv(lstar) * rxx;
%   B = inv(lstar);
%   B(:,sys.ncx+1:end) = [];
  
  
  % Remove unconnected coils from the circuit equations
  ccnames = {'OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', ...
    'PF5', 'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL'};
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
  
  if saveit
    A = nstxu_sys.amat;
    B = nstxu_sys.bmat;
    sys.A = A;
    sys.As = removeFirstEigval(A);
    sys.B = B;
    sys.inputs = ccnames;
    sys.states = [ccnames cellstr(vvnames)' {'IP'}];
    sys.lstar = nstxu_sys.lstar;
    fn = [savedir num2str(shot) '_' num2str(time_ms) '_sys.mat'];
    save(fn, 'sys')
  end
end















