function y = measure_y(psizr, currents, targets, tok_data_struct)

struct_to_ws(tok_data_struct);

psicp = bicubicHermite(rg, zg, psizr, targets.rcp, targets.zcp);  
[psibry, psibry_r, psibry_z] = bicubicHermite(rg, zg, psizr, targets.rbdef, targets.zbdef);  

y = [currents; psicp(:) - psibry; psibry_r; psibry_z];


















