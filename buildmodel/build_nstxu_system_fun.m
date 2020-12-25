% EXAMPLE:
% time_ms = 90;
% load('eq204660_90.mat')
% load('nstxu_obj_config2016_6565.mat')
% opts.islimited = time_ms < 220;
% opts.Te_res = min( time_ms/0.3, 4000);
% build_nstxu_system_fun(eq, tok_data_struct, opts)

function sys = build_nstxu_system_fun(eq, tok_data_struct, opts)

% some coils turned off for previous campaign, remove from circuit eqns
remove_coils = {'PF1BU', 'PF1BL', 'PF1CU', 'PF1CL', 'PF4'};  

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

% specify that plasma resistance is much larger when limited
if isfield(opts, 'islimited') && opts.islimited
  build_inputs.Rp = 2.44e-8 * 120;
end
if isfield(opts, 'islimited') && ~opts.islimited && isfield(opts, 'Te_res')
  build_inputs.Te_res = opts.Te_res;
end

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
  nstxu_sys.bmat(iremove,:) = [];
  nstxu_sys.amat(:,iremove) = [];
  nstxu_sys.bmat(:,iremove) = [];
  ccnames = ccnames(~iremove);
end
  
% Ai = c2d(nstxu_sys.amat, nstxu_sys.bmat, .01);
% Ai(end,end)

% write to outputs
sys.A = nstxu_sys.amat;
sys.As = removeFirstEigval(nstxu_sys.amat);
sys.B = nstxu_sys.bmat;
sys.inputs = ccnames;
sys.states = [ccnames cellstr(vvnames)' {'IP'}];
end















