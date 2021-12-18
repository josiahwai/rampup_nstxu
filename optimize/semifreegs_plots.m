
function semifreegs_plots(eq, init, target, tok_data_struct)

struct_to_ws(tok_data_struct);
psin = linspace(0,1,nr);

% plot profiles
try 
  figure
  hold on
  plot(psin, eq.pprime)
  plot(psin, init.pprime)
  legend('Estimate', 'EFIT', 'fontsize', 18)
catch
end


try
  figure
  hold on
  plot(psin, eq.pres)
  plot(psin, init.pres)
  legend('Estimate', 'EFIT', 'fontsize', 18)
catch
end


try
  figure
  hold on
  plot(psin, eq.ffprim)
  ylabel("FF'", 'fontsize', 14, 'fontweight', 'bold')
  plot(psin, init.ffprim),
  legend('Estimate', 'EFIT', 'fontsize', 18)
catch
end


try
  figure
  subplot(121)
  [~,cs] = contourf(eq.jphi, 10);
  colorbar
  title('Estimate')
  cax = caxis;
  subplot(122)
  contourf(init.jphi, cs.LevelList)
  colorbar
  title('EFIT')
  caxis(cax);
catch
end


try    
  
  figure
  hold on

  % psibry = bicubicHermite(rg, zg, init.psizr, target.rbdef, target.zbdef);
  contour(rg, zg, init.psizr, [init.psibry init.psibry], 'b')

  %psibry = bicubicHermite(rg, zg, eq.psizr_pla + psi_app, target.rbdef, target.zbdef);
  
  psi_app = reshape(tok_data_struct.mpc*init.ic + tok_data_struct.mpv*init.iv, nr, nz);
  contour(rg, zg, eq.psizr_pla + psi_app, [init.psibry init.psibry], '--r')

  contour(rg, zg, eq.psizr, [eq.psibry eq.psibry], 'g')
  
  plot(tok_data_struct.limdata(2,:), tok_data_struct.limdata(1,:), 'k')
  
  
  axis equal
  set(gcf, 'Position', [992 183 431 622])
  scatter(target.rcp, target.zcp, 'k', 'filled')

catch
end












