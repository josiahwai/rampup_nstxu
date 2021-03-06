% Import an equilibrium for nstxu from the mdsplus tree
%
% This builds off the function read_eq.m included with Toksys, but modifies
% a few fields so that the coil current vectors are consistent with
% geometry. 
%
% Josiah Wai, 3/1/2021

% EXAMPLE:
% shot = 204660;
% time = 0.200;
% tree = 'EFIT01';
% tokamak = 'nstxu';
% server = 'skylark.pppl.gov:8501';
% opts.save_local_copy = 0;
% opts.save_dir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk_import2';
% opts.plotit = 1;

function eq = read_eq_nstxu(shot, time, tree, tokamak, server, opts)


eq = read_eq(shot,time,tree,tokamak,server);

% The coil current vector eq.cc is confusing and dependent on geometry in
% arcane ways. Instead of inverting eq.cc to eq.ic, we will load ic directly
% from mds, define our geometry here, and calculate what eq.cc is for
% our geometry.  

tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);

signal = mds_fetch_signal(shot, tree, time, '.RESULTS.AEQDSK:ECCURT', opts.plotit);  % OH
ecefit = signal.sigs;

signal = mds_fetch_signal(shot, tree, time, '.RESULTS.AEQDSK:CCBRSP', opts.plotit); % PF coils + vessel 
ccefit = signal.sigs;

icx = zeros(13,1);
ivx = zeros(40,1);

icx(1) = ecefit;    % OH
icx(2) = ccefit(1); % PF1AU
icx(3) = ccefit(2); % PF1BU
icx(4) = ccefit(3); % PF1CU
icx(5) = ccefit(4); % PF2U
icx(6) = ccefit(5); % PF3U
icx(7) = (ccefit(6)+ccefit(9))/2.0;  % PF4
icx(8) = (ccefit(7)+ccefit(8))/2.0;  % PF5
icx(9) = ccefit(10);  % PF3L
icx(10) = ccefit(11); % PF2L
icx(11) = ccefit(12); % PF1CL
icx(12) = ccefit(13); % PF1BL
icx(13) = ccefit(14); % PF1AL

ivx(:) = ccefit(15:end);

ic = circ.Pcc * icx;
iv = circ.Pvv * ivx;

% copy circuit information into eq
eq.ecturn = tok_data_struct.ecnturn;
eq.ecid   = ones(size(eq.ecturn));

eq.turnfc = tok_data_struct.fcnturn';
eq.fcturn = circ.fcfrac;
eq.fcid = circ.fccirc';

% Convert from ic to cc:
% The ic vector already represents coil currents in toksys format. Use a
% hack to convert from ic to efit format, so that gs codes can convert
% back from efit format to toksys.
idiot = eq;
ncc = length(ic);
for j = 1:ncc
  idiot.cc = zeros(ncc,1);
  idiot.cc(j) = 1;
  equil_I = cc_efit_to_tok(tok_data_struct,idiot);
  iccc(:,j) = equil_I.cc0t;
end
piccc = pinv(iccc);
cc = piccc * ic;

eq = append2struct(eq, cc, ic, iv, icx, ivx);

% error check - convert cc back to ic
equil_I = cc_efit_to_tok(tok_data_struct, eq);
cc = equil_I.cc0t;

if max(abs(cc - ic)) > sqrt(eps)
  error('Something is wrong with coil current vectors.')
end


if opts.save_local_copy
  fn = [opts.save_dir '/eq' num2str(shot) '_' num2str(floor(time*1000)) '.mat'];
  save(fn, 'eq')
end

if opts.plotit
  figure
  plot_eq(eq)
end
  


























