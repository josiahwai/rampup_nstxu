
function [spec, init, config] = make_gsdesign_inputs2(x, tok_data_struct, eq, circ, opts)
clear spec init config 

if ~exist('opts','var')
  opts.pres = eq.pres;
  opts.fpol = eq.fpol;
end

  
init = eq;

config = tok_data_struct;
config.max_iterations = 12;
config.constraints = 1;
config.nkn = config.nr-1;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

spec.targets.rsep = init.rbbbs;
spec.targets.zsep = init.zbbbs;
spec.weights.sep = ones(length(spec.targets.rsep),1) * 1e-5;

[rx1, zx1] = isoflux_xpFinder(eq.psizr, 0.6, -1.05, eq.rg, eq.zg);
% spec.targets.rbdef = rx1;
% spec.targets.zbdef = zx1;
% spec.weights.bdef = 1e-2;


spec.targets.psibry = init.psibry;
spec.weights.psibry = 1e-3;

config.mpp = load('nstxu_obj_config2016_6565.mat').tok_data_struct.mpp;

gs_configure
gs_initialize
gs_eq_analysis

spec.targets.li = li;
spec.weights.li = 0;

% spec.targets.betap = betap;
% spec.weights.betap = 0.01;

config.pres0 = init.pres;
config.fpol0 = init.fpol;


% spec.targets.pres = load('pres.mat').pres;
% spec.targets.fpol = load('fpol.mat').fpol;
% spec.weights.pres = ones(size(pres)) * 0.01;
% spec.weights.fpol = ones(size(fpol)) * 0.01;


spec.targets.pres = opts.pres;
spec.weights.pres = ones(size(opts.pres)) * 1e2;
spec.targets.fpol = opts.fpol;
spec.weights.fpol = ones(size(opts.fpol)) * 1e3;


is = circ.Pxx * x;
ic = is(circ.iic);
iv = is(circ.iiv);
ip = is(circ.iip);

spec.cccirc = circ.cccirc;
spec.locks.ic = ic;
spec.locks.iv = iv;
spec.locks.cpasma = ip;


















