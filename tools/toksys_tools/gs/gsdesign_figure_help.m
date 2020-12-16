function gsdesign_figure_help

%  USAGE:   gsdesign_figure_help
%
%  PURPOSE: Explain the gsdesign figure
%
%  INPUTS: ha, structure with information used by gsdesign_plot_progress
%
%  OUTPUTS: Explanations of each subplot
%	
	
%
%  WRITTEN BY:  Anders Welander ON 2015-10-14
%
%  MODIFICATION HISTORY:			
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ha.figure = figure;
fpos = get(ha.figure,'position');
fpos(4) = fpos(3)/2;
fpos(3:4) = [1000 500];
set(ha.figure,'position',fpos)
if isa(ha.figure,'matlab.ui.Figure')
  ha.figure = ha.figure.Number;
end
zoom on

ha.running = 1; % Flag that code is running, for stop button
ha.h1 = []; % handles to objects with font size fs1
ha.h2 = []; % handles to objects with font size fs2
ha.f1 = 0.8/60;  % fs1/window width in pixels
ha.f2 = 0.8/75; % fs2/window width in pixels
ha.fs = 0.8/60;  % fs1/window width in pixels

% CONDUCTORS
ha.hs(1) = subplot(2,4,1,'Position',[0.05 0.6 0.2 0.3],...
  'XTick',[],'YTick',[],'Box','on');
hold on
plot(1,1,'yo--','MarkerSize',18,'MarkerFaceColor','yellow','LineWidth',3)
plot(1,1,'bx','MarkerSize',12,'LineWidth',3)
plot(1,1,'go','MarkerSize',4,'MarkerFaceColor','green','LineWidth',3)
ha.h2(end+1) = text(1,1,'   Coil at lower limit');
plot(1,2,'bx','MarkerSize',12,'LineWidth',3)
ha.h2(end+1) = text(1,2,'   Lower limit');
plot(1,3,'go','MarkerSize',4,'MarkerFaceColor','green','LineWidth',3)
ha.h2(end+1) = text(1,3,'   Coils (vessel if ~=0)');
plot(1,4,'rx','MarkerSize',12,'LineWidth',3)
ha.h2(end+1) = text(1,4,'   Upper limit');
plot(1,5,'yo--','MarkerSize',18,'MarkerFaceColor','yellow','LineWidth',3)
plot(1,5,'rx','MarkerSize',12,'LineWidth',3)
plot(1,5,'go','MarkerSize',4,'MarkerFaceColor','green','LineWidth',3)
ha.h2(end+1) = text(1,5,'   Coil at upper limit');
ha.h1(end+1) = title('[coils; vessel]');
ha.h1(end+1) = ylabel('Current [A]');
ha.h1(end+1) = xlabel('Index in [ic;iv]',...
  'VerticalAlignment','top');
axis([.9 2 0 6])

% pprime, ffprim
ha.hs(2) = subplot(2,4,2,'Position',[0.3 0.6 0.2 0.3],...
  'XTick',[],'YTick',[],'Box','on');
try % needed for matlab2016b since it behaves badly but won't work for 2014
ha.h2(end+1) = text(ha.hs(2),1,6,'  \psi = poloidal flux');
ha.h2(end+1) = text(ha.hs(2),1,5,'  \psi_a = \psi on axis');
ha.h2(end+1) = text(ha.hs(2),1,4,'  \psi_b = \psi at boundary');
ha.h2(end+1) = text(ha.hs(2),1,3,'  p'' = 2\pi dp/d\psi');
ha.h2(end+1) = text(ha.hs(2),1,2,'  ff'' = 2\pi d(f^2/2)/d\psi');
ha.h2(end+1) = text(ha.hs(2),1,1,'  R = Axis radius');
catch
ha.h2(end+1) = text(1,6,'  \psi = poloidal flux');
ha.h2(end+1) = text(1,5,'  \psi_a = \psi on axis');
ha.h2(end+1) = text(1,4,'  \psi_b = \psi at boundary');
ha.h2(end+1) = text(1,3,'  p'' = 2\pi dp/d\psi');
ha.h2(end+1) = text(1,2,'  ff'' = 2\pi d(f^2/2)/d\psi');
ha.h2(end+1) = text(1,1,'  R = Axis radius');
end
ha.h1(end+1) = xlabel(ha.hs(2),'(\psi-\psi_a)/(\psi_b-\psi_a)',...
  'VerticalAlignment','top');
ha.h1(end+1) = ylabel(ha.hs(2),'j_\phi [A/m^2]');
ha.h1(end+1) = title(ha.hs(2),'Profiles: {\color{blue}Rp''}, {\color{red}ff''/\mu_0/R}');
axis(ha.hs(2),[.9 2 0 7])

% MACHINE GEOMETRY
ha.hs(3) = subplot(2,4,4,'Position',[0.55 0.1 0.4 0.8],...
  'XTick',[],'YTick',[],'Color',[.8 .8 .8],'Box','on');
hold on
goldenrod = [218 165 32]/255;
plot(0.95,11,'bo','MarkerSize',6,'LineWidth',2);
plot(1.05,11,'ro','MarkerSize',6,'LineWidth',2);
plot(1.00,11,'mo','MarkerSize',6,'LineWidth',2);
ha.h2(end+1) = text(1.05,11,[' {\color{blue}Limit 1},', ...
' {\color{magenta}target}, {\color{red}limit 2} for separatrix point']);
plot(1,10,'o','MarkerSize',11,'LineWidth',4, 'Color',goldenrod);
ha.h2(end+1) = text(1,10,'      Locked separatrix point');
plot(1,9,'yx','MarkerSize',18,'LineWidth',3);
plot(1+[-1 -1 1 1 -1]/20,9+[-1 1 1 -1 -1]/2,'--y','LineWidth',3);
ha.h2(end+1) = text(1,9,'      Target and polygon-limits for x-point');
plot(1,8,'x','MarkerSize',18,'LineWidth',3, 'Color',goldenrod);
ha.h2(end+1) = text(1,8,'      Locked x-point');
plot(1,7,'bx','MarkerSize',12,'LineWidth',2);
ha.h2(end+1) = text(1,7,'      x-point');
plot(1+[-1 -1 1 1 -1]/20,6+[-1 1 1 -1 -1]/2,'--r','LineWidth',3);
ha.h2(end+1) = text(1,6,'      polygon-limits for no-x region');
plot(1,5,'w*','MarkerSize',24,'LineWidth',4);
ha.h2(end+1) = text(1,5,'      Target snowflake');
plot(1,4,'*','MarkerSize',24,'LineWidth',4, 'Color',goldenrod);
ha.h2(end+1) = text(1,4,'      Locked snowflake');
plot(1,3,'y+','MarkerSize',12,'LineWidth',2);
plot(1,3,'yo','MarkerSize',18,'LineWidth',3);
ha.h2(end+1) = text(1,3,'      Target boundary-defining point');
plot(1,2,'+','MarkerSize',12,'LineWidth',2, 'Color',goldenrod);
plot(1,2,'o','MarkerSize',18,'LineWidth',3, 'Color',goldenrod);
ha.h2(end+1) = text(1,2,'      Locked boundary-defining point');
plot(1,1,'yx','MarkerSize',24,'LineWidth',5);
plot(1,1,'rx','MarkerSize',21,'LineWidth',3);
ha.h2(end+1) = text(1,1,'      Boundary-defining point');
plot(1,0,'ro','MarkerSize',14,'LineWidth',1,'MarkerFaceColor','red');
plot(1,0,'ws','MarkerSize',4,'MarkerFaceColor','white');
ha.h2(end+1) = text(1,0,'      Wrong way (if bdef target or lock)');
plot(1,-1,'kp','MarkerSize',18,'LineWidth',3,'MarkerFaceColor','y');
ha.h2(end+1) = text(1,-1,'      Target flux expansion point');
plot(1,-2,'kp','MarkerSize',18,'LineWidth',3,'MarkerFaceColor',goldenrod);
ha.h2(end+1) = text(1,-2,'      Locked flux expansion point');
ha.h1(end+1) = title('Symbols');
axis([.9 2 -3 12])

% ERROR VECTOR
ha.hs(5) = subplot(2,4,5,'Position',[0.05 0.1 0.2 0.3],...
  'XTick',[],'YTick',[],'Box','on');
hold on
plot(1,5,'go','LineWidth',3);
ha.h2(end+1) = text(1,5,'   weight\cdot(value-target)');
plot(1,4,'ro','LineWidth',3);
ha.h2(end+1) = text(1,4,'   biggest positive error');
plot(1,3,'bo','LineWidth',3);
ha.h2(end+1) = text(1,3,'   biggest negative error');
ha.h2(end+1) = text(1,1,'Zoom in close on green\newlinesymbols to see names');
ha.h1(end+1) = xlabel('Index in error vector',...
  'VerticalAlignment','top');
ha.h1(end+1) = title('Error vector');
axis([.9 2 0 6])

% FLUX ERROR
ha.hs(6) = subplot(2,4,6,'Position',[0.3 0.1 0.2 0.3],...
  'XTick',[],'YTick',[],'Box','on');
ha.h2(end+1) = text(1,4,'j_\phi(r,z) derived from \psi, p'', ff''');
ha.h2(end+1) = text(1,3,'\psi_p_l_a = flux from plasma, j_\phi');
ha.h2(end+1) = text(1,2,'\psi_a_p_p = flux from conductors');
ha.h1(end+1) = text(1,0.75,'\psi_e_r_r = \psi - \psi_p_l_a - \psi_a_p_p');
ha.h1(end+1) = ylabel('\psi_e_r_r/(\psi_b-\psi_a)');
ha.h1(end+1) = xlabel('Indices 1:nz',...
  'VerticalAlignment','top');
ha.h1(end+1) = title('Numerical error, \psi_e_r_r');
axis([.975 2 0 5])

% RESIZE FUNCTION
set(ha.figure,'UserData',ha)
set(ha.figure,'ResizeFcn','gs_resize_fonts(get(gcbo,''UserData''))')
gs_resize_fonts(ha)

set(ha.figure,'Name','Explanations')

% In matlab2016b the subplots need to be put back into their positions
set(ha.hs(1),'Position',[0.05 0.6 0.2 0.3])
set(ha.hs(2),'Position',[0.3 0.6 0.2 0.3])
set(ha.hs(3),'Position',[0.55 0.1 0.4 0.8])
set(ha.hs(5),'Position',[0.05 0.1 0.2 0.3])
set(ha.hs(6),'Position',[0.3 0.1 0.2 0.3])
