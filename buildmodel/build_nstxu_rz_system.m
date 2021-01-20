clc; close all

% ========
% SETTINGS
% ========
shot = 204660;
times = 60:10:200;
saveit = 1;

savedir = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/rz/';
eqdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk_import/';
load('nstxu_obj_config2016_6565.mat')

% some coils turned off for previous campaign, remove from circuit eqns
remove_coils = {'PF1BU', 'PF1BL', 'PF1CU', 'PF1CL', 'PF4'};  

% ===========
% Build loop
% ===========
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
  build_inputs.cccirc = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 ...
      11 12 13];
  
  circ = nstxu2016_circ(tok_data_struct);
  build_inputs.vvcirc = circ.vvcirc;
  
  % custom plasma resistance
  load('fit_Rp.mat');
  build_inputs.Rp = fit_Rp(time_ms/1000);  
  % build_inputs.Rp = Rp(itime);

  nstxu_sys = build_tokamak_system(build_inputs); 
  delete('NSTXU_netlist.dat')
  
  
  
  
  
%   [rzrig_data] = rzrig(eq,'NSTXU',tok_data_struct, 0, 0, [], 0, 0);
%   
%   % Now we do some modification from the GA codes to consider rcur and zcur 
%   % as part of the control inputs. 
%   
%   dfxdr = [rzrig_data.dfsdr; rzrig_data.dfpdrC];
%   dfxdz = [rzrig_data.dfsdz; rzrig_data.dfpdzC];
%   rxx = nstxu_sys.rxx;
%   Pxx = nstxu_sys.Pxx;
%   mxx = nstxu_sys.mxx;
%   mstar = Pxx'*mxx*Pxx;
%   mstari = inv(mstar);
%   
%   A = -mstari * Pxx'*rxx*Pxx;
%   B = [mstari(:,1:nstxu_sys.ncx) -mstari*Pxx'*dfxdr -mstari*Pxx'*dfxdz];
%       
% %   lstar = sys.lstar;
% %   rxx = sys.Pxx' * sys.rxx * sys.Pxx;
% %   drdix = sys.Pxx \ [sys.gspert_data.drdis sys.gspert_data.drdip]';
% %   dzdix = sys.Pxx \ [sys.gspert_data.dzdis sys.gspert_data.dzdip]';
% %   
% %   Ip0 = eq.cpasma;
% %   lstar(:,end+1) = 1./drdix;
% %   lstar(:,end+1) = 1./dzdix;
% %   lstar(end+1,:) = [drdix' -1 0];
% %   lstar(end+1,:) = [dzdix' 0 -1];
% %   rxx = blkdiag(rxx, zeros(2));
% %   
% %   A = -inv(lstar) * rxx;
% %   B = inv(lstar);
% %   B(:,sys.ncx+1:end) = [];
  

  A = nstxu_sys.amat;
  B = nstxu_sys.bmat;

  % Remove unconnected coils from the circuit equations
  ccnames = {'OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', ...
    'PF5', 'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL'};
  if exist('remove_coils', 'var') && ~isempty('remove_coils')
    iremove = zeros(size(ccnames));
    for i = 1:length(remove_coils)
      iremove = iremove | strcmp(ccnames, remove_coils{i});
    end
    A(iremove,:) = [];
    A(:,iremove) = [];
    B(:,iremove) = [];
    B(iremove,:) = [];
    ccnames = {ccnames{~iremove}};
  end        
  
  if saveit
    sys.A = A;
    sys.As = removeFirstEigval(A);
    sys.B = B;
    sys.inputs = ccnames;
    sys.states = [ccnames cellstr(circ.vvnames)' {'Ip', 'rcur', 'zcur'}];    
    fn = [savedir num2str(shot) '_' num2str(time_ms) '_sys.mat'];
    if ~exist(savedir, 'dir'), mkdir(savedir); end
    save(fn, 'sys')
  end
end















