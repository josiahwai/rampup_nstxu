
function y = read_isoflux(eq, target, tok_data_struct)

rg = tok_data_struct.rg;
zg = tok_data_struct.zg;
circ = nstxu2016_circ(tok_data_struct);

y.psibry = eq.psibry;
y.psicp = bicubicHermite(rg, zg, eq.psizr, target.rcp(:), target.zcp(:));
[y.target_bdef_psi, y.target_bdef_psi_r, y.target_bdef_psi_z] = bicubicHermite(rg, zg, eq.psizr, target.rbdef, target.zbdef);
y.ip = eq.cpasma;
y.icx = pinv(circ.Pcc) * eq.ic;
% y.zcur = eq.zcur;

y.psi_in_gap = bicubicHermite(rg, zg, eq.psizr, target.r_ingap, target.z_ingap);
[y.psix, y.psix_r, y.psix_z] = bicubicHermite(rg, zg, eq.psizr, target.rx, target.zx);








































