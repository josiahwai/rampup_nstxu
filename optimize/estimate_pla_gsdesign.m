
function eq = estimate_pla_gsdesign(target, tok_data_struct, init, opts)    

clear spec config 

circ = nstxu2016_circ(tok_data_struct);

config = tok_data_struct;
config.max_iterations = opts.max_iterations;
config.constraints = 1;
config.nkn = config.nr-1;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;
config.mpp = load('nstxu_obj_config2016_6565.mat').tok_data_struct.mpp;

spec.targets.rsep = target.rcp;
spec.targets.zsep = target.zcp;
spec.weights.sep = ones(length(spec.targets.rsep),1) * 1e2;

spec.targets.rsep(end+1) = target.rbdef;
spec.targets.zsep(end+1) = target.zbdef;
spec.weights.sep(end+1) = 1e2;

if ~target.islimited
  spec.targets.rx = target.rbdef;
  spec.targets.zx = target.zbdef;
  spec.weights.x = 1e2;
  spec.targets.rbdef = target.rbdef;
  spec.targets.zbdef = target.zbdef;
end


spec.targets.li = target.li;
spec.weights.li = 1e2;

spec.targets.betap = target.betap;
spec.weights.betap = 1e2;

config.pres0 = init.pres;
% config.fpol0 = init.fpol;

spec.cccirc = circ.cccirc;

% spec.targets.ic = init.ic;    % nan(circ.nc,1); 
% spec.targets.iv = init.iv;
% spec.locks.iv = nan(circ.nv,1);
% spec.targets.iv = zeros(circ.nv,1);
% spec.weights.iv = ones(circ.nv,1)*1;

spec.locks.cpasma = target.ip;


eq = gsdesign(spec,init,config);

% figure
% plot(eq.ffprim)
% hold on
% plot(init.ffprim)

% figure
% hold on
% [~,cs] = contour(eq.psizr, 20, '--');
% contour(init.psizr, cs.LevelList, '-');






















