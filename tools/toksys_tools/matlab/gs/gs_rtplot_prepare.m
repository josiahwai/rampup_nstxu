%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_rtplot_prepare
%
%  PURPOSE: Prepare data for gs_rtplot
%           
%  INPUTS:  
%	
%  OUTPUTS: Data that will be in the rtplot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander ON 2015-03-08
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_counter = plot_counter+1;
while plot_counter < length(plot_times) && ...
      time >= plot_times(plot_counter+1)
  plot_counter = plot_counter+1;
end
if isfield(plot_settings,'flux_contours') && ...
  plot_settings.flux_contours.LineWidth
  psizr(:) = lae.psizr(:) + dpsizrdx*dx;
end
if plasma % Update quantities that gs_rtplot may need
  if equilibrium_update_is_complete
    rbbbs = lae.rbbbs + drbdx*dx;
    zbbbs = lae.zbbbs + dzbdx*dx;
    rbdef = lae.rbdef + drbdefdx*dx;
    zbdef = lae.zbdef + dzbdefdx*dx;
    psimag = lae.psimag + dpsimagdx*dx;
    psibry = lae.psibry + dpsibrydx*dx;
    cpasma = lae.cpasma + dcpasmadx*dx;
    li = lae.li + dlidx*dx;
    betap = lae.betap + dbetapdx*dx;
  else % Contend with last analyzed equilibrium
    rbbbs = lae.rbbbs;
    zbbbs = lae.zbbbs;
    rbdef = lae.rbdef;
    zbdef = lae.zbdef;
    psimag = lae.psimag;
    psibry = lae.psibry;
    cpasma = lae.cpasma;
    li = lae.li;
    betap = lae.betap;
    if isfield(plot_settings,'flux_contours') && ...
      plot_settings.flux_contours.LineWidth
      psizr = lae.psizr;
    end
  end
  if plasma_is_tiny
    a0 = cpasma/lae.cpasma*lae.aminor;
    v = linspace(0,2*pi,nbbbs_max)';
    rbbbs = rbdef + a0*vlim(1) + a0*cos(v);
    zbbbs = zbdef + a0*vlim(2) + a0*sin(v);
    nbbbs = nbbbs_max;
  end
end
