function y = measure_y2(psizr, init, targets, tok_data_struct)

struct_to_ws(tok_data_struct);

psicp = bicubicHermite(rg, zg, psizr, targets.rcp(:), targets.zcp(:));  
[psibry, psibry_r, psibry_z] = bicubicHermite(rg, zg, psizr, targets.rbdef, targets.zbdef);  

psitouch = bicubicHermite(rg, zg, psizr, targets.rtouch, targets.ztouch);  
[psixlo, psixlo_r, psixlo_z] = bicubicHermite(rg, zg, psizr, targets.rx_lo, targets.zx_lo);  
[psixup, psixup_r, psixup_z] = bicubicHermite(rg, zg, psizr, targets.rx_up, targets.zx_up);  


icx = init.icx;
ivx = init.ivx;
ip = init.cpasma;
diff_psicp_psitouch = psicp - psitouch;
diff_psicp_psixlo = psicp - psixlo;
diff_psicp_psixup = psicp - psixup;


y = variables2struct(icx, ivx, ip, psicp, diff_psicp_psitouch, diff_psicp_psixlo, ...
  diff_psicp_psixup, psixlo_r, psixlo_z, psixup_r, psixup_z, psitouch, psixlo, psixup);

















