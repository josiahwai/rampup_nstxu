function eq = gsdesign_fit(ref, tok_data_struct)

clear spec init config gsdesign

circ = nstxu2016_circ(tok_data_struct);

ROOT = getenv('RAMPROOT');
opts.cache_dir = [ROOT '/fetch/cache/'];
if ref.islimited
  init = fetch_eqs_nstxu(204660, 0.085, 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts);
else
  init = fetch_eqs_nstxu(204660, 0.5, 'EFIT01', 'nstxu', 'skylark.pppl.gov:8501', opts);
end
init = struct_fields_to_double(init.gdata);

config = tok_data_struct;
config.max_iterations = 10;
config.constraints = 1;
config.nkn = config.nr-1;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

spec.targets.rsep = ref.rbbbs{:};
spec.targets.zsep = ref.zbbbs{:};
spec.weights.sep = ones(length(spec.targets.rsep),1) * 10;


spec.targets.li = ref.li;
spec.weights.li = 1;

spec.targets.Wth = ref.wmhd;
spec.weights.Wth = 1e-2;

config.pres0 = init.pres;
config.fpol0 = init.fpol;


spec.locks.ci = ref.icx;
% spec.weights.ci = ones(size(spec.targets.ci)) * 1e-4;

spec.locks.ci(circ.iremove) = 0;

spec.cccirc = circ.cccirc;

if ~any(isnan(ref.ivx))
  spec.locks.iv = circ.Pvv * ref.ivx(:);
end

spec.limits.ic = circ.limits.ic;

spec.targets.cpasma = ref.ip;
spec.weights.cpasma = 1;

eq = gsdesign(spec, init, config);

eq.icx = pinv(circ.Pcc) * eq.ic;
eq.ivx = pinv(circ.Pvv) * eq.iv;

















