
shot = 204660;
time = 0.4;

eq = fetch_eq_nstxu(shot, time);

vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;

tok_data_struct = vac_sys.build_inputs.tok_data_struct;
tok_data_struct.imks = 1;

circ = nstxu2016_circ(tok_data_struct);

ibad = (eq.rbbbs == 0 & eq.zbbbs == 0);
eq.rbbbs(ibad) = [];
eq.zbbbs(ibad) = [];

targ_geo.cp.n = 10;
[targ_geo.cp.r, targ_geo.cp.z, s] = interparc(eq.rbbbs, eq.zbbbs, targ_geo.cp.n, 1, 0);

clf
plot_eq(eq)
scatter(targ_geo.cp.r, targ_geo.cp.z, 100, 'filled')


cmat_data = build_cmat_data(eq, circ, tok_data_struct, targ_geo)


























