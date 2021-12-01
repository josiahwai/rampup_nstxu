
function y = read_isoflux(eq, target, tok_data_struct)

rg = tok_data_struct.rg;
zg = tok_data_struct.zg;
circ = nstxu2016_circ(tok_data_struct);

y.psibry = eq.psibry;
y.psicp = bicubicHermite(rg, zg, eq.psizr, target.rcp(:), target.zcp(:));
[y.target_bdef_psi, y.target_bdef_psi_r, y.target_bdef_psi_z] = bicubicHermite(rg, zg, eq.psizr, target.rbdef, target.zbdef);
y.ip = eq.cpasma;
y.icx = eq.icx;
% y.zcur = eq.zcur;











































