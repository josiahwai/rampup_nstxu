function [spec,init,config] = make_gs_inputs(x, profiles, init, tok_data_struct)

clear spec config

circ = nstxu2016_circ(tok_data_struct);

ic = circ.Pcc * x(circ.iicx); % coil currents
iv = circ.Pvv * x(circ.iivx); % vessel currents
ip = x(circ.iipx);  % plasma current

config = tok_data_struct;
config.max_iterations = 10;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

% Constrain profiles to three degrees of freedom
config.constraints = 1;

% Take profile details from init
config.pres0 = profiles.pres;
config.fpol0 = profiles.fpol;

spec.targets.pres = profiles.pres;
spec.weights.pres(1:length(spec.targets.pres)) = 0.1;
spec.targets.fpol = profiles.fpol;
spec.weights.fpol(1:length(spec.targets.fpol)) = 100;

% A generous number of knots ensures the profile can be made in gsdesign
config.psikn = (0:config.nr-1)/(config.nr-1);

% load circuit connection info
circ = nstxu2016_circ(tok_data_struct);
spec.cccirc = circ.cccirc;

% apply limits
spec.limits.ic = circ.limits.ic;

% lock coil currents
spec.locks.ic = ic;
spec.locks.iv = iv;

% lock Ip
spec.locks.cpasma = ip;

spec.targets.rsep = init.rbbbs;
spec.targets.zsep = init.zbbbs;
spec.weights.sep = 1e-2;




















