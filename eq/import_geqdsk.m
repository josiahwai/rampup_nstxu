% Import gfile and solve for the gsdesign equilibrium that best matches the 
% gfile equilibrium. 

% EXAMPLE
% shot = 204660;
% time_ms = 120;
% eq = import_gfile(shot, time_ms, 0, 1)

function eq = import_geqdsk(shot, time_ms, shotdir, savedir, saveit, plotit)

load('nstxu_obj_config2016_6565.mat')
config = tok_data_struct;

% load gfile
dum = read_eq(shot, time_ms/1000, shotdir, 'NSTX');
eq = dum.gdata;
    
% copy circuit information into eq
eq.ecturn = tok_data_struct.ecnturn;
eq.ecid   = ones(size(eq.ecturn));

circ = nstxu2016_circ(tok_data_struct);
eq.turnfc = tok_data_struct.fcnturn';
eq.fcturn = circ.fcfrac;
eq.fcid = circ.fccirc;

% import coil currents
load(['coils' num2str(shot) '.mat'])
i = find( floor(coils.t*1000) == time_ms);
if isempty(i)
  warning('Could not import coil currents')
end
ic = coils.ic(i,:)';
iv = coils.iv(i,:)';
ic = ic(circ.cccirc);
iv = iv(circ.vvcirc) .* circ.vvfrac(:); 

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
eq.cc = piccc * ic;
eq.ic = ic;
eq.iv = iv;

% Analyze and copy some equilibrium properties
init = eq;
gs_configure
gs_initialize
gs_eq_analysis
eq.rcur = rcur;
eq.zcur = zcur;
eq.betap = betap;
eq.li = li;
eq.rbdef = rbdef;
eq.zbdef = zbdef;


if plotit, figure; plot_eq(eq); end

if saveit
  fn = [savedir 'eq' num2str(shot) '_' num2str(time_ms) '.mat'];
  save(fn, 'eq');
end


end















