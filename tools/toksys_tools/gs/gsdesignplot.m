%  USAGE:   gsdesignplot
%
%  PURPOSE: Plot gsdesign iterations
%
%  INPUTS: pprime, ffprim
%          cpasma, total plasma current
%          betap, poloidal beta
%          li, normalized inductance
%          psimag, the flux at the magnetix axis
%          psibry, the flux at the boundary
%          nulls, a structure with nulls that may affect boundary
%          psizr, the flux at nz vertical positions * nr radial positions
%          rbdef, zbdef, point that defines the boundary (x or touch)
%          ic, iv, currents in coils and vessel
%
%  OUTPUTS: plots of conductor currents, pprime & ffprim, scalars,
%             contours, weights*errors, flux error
%	
%  FEATURES: A resize function (ResizeFcn) scales fonts with window size 
%            Zooming in on errors reveals their individual names 
	
%
%  WRITTEN BY:  Anders Welander  ON	4/12/13
%
%  MODIFICATION HISTORY:			
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Configure plot_settings
if ~exist('plot_settings','var')
  plot_settings = [];
end

% Maximum number of x-points to show and report positions for
if ~isfield(plot_settings,'nxptmax')
  plot_settings.nxptmax = 9;
end

% Number of flux surface contours from axis to boundary
if ~isfield(plot_settings,'nflux')
  plot_settings.nflux = 8;
end

if ~isfield(plot_settings,'SOL')
  plot_settings.SOL = [];
end

% Number of contours to show for scrape-off layer
if ~isfield(plot_settings.SOL,'n')
  plot_settings.SOL.n = 0;
end

% Distance between SOL contours in outboard midplane
if ~isfield(plot_settings.SOL,'d')
  plot_settings.SOL.d = 0.001;
end

if ~isfield(plot_settings,'SymbolSize')
  plot_settings.SymbolSize = [];
end

if isempty(plot_settings.SymbolSize)
  plot_settings.SymbolSize = 1;
end

if ~isstruct(plot_settings.SymbolSize)
  dum = plot_settings.SymbolSize;
  plot_settings = rmfield(plot_settings,'SymbolSize');
  plot_settings.SymbolSize.locks.x = dum;
  plot_settings.SymbolSize.locks.sep = dum;
  plot_settings.SymbolSize.locks.snf = dum;
  plot_settings.SymbolSize.locks.fluxexp = dum;
  plot_settings.SymbolSize.locks.bdef = dum;
  plot_settings.SymbolSize.limits.x = dum;
  plot_settings.SymbolSize.limits.nox = dum;
  plot_settings.SymbolSize.limits.sep = dum;
  plot_settings.SymbolSize.limits.snf = dum;
  plot_settings.SymbolSize.limits.bdef = dum;
  plot_settings.SymbolSize.targets.x = dum;
  plot_settings.SymbolSize.targets.sep = dum;
  plot_settings.SymbolSize.targets.snf = dum;
  plot_settings.SymbolSize.targets.fluxexp = dum;
  plot_settings.SymbolSize.targets.bdef = dum;
  plot_settings.SymbolSize.x = dum;
  plot_settings.SymbolSize.bdef = dum;
  plot_settings.SymbolSize.wrongway = dum;
end

if ~isfield(plot_settings.SymbolSize,'locks')
  plot_settings.SymbolSize.locks = [];
end
if ~isfield(plot_settings.SymbolSize.locks,'x')
  plot_settings.SymbolSize.locks.x = 1;
end
if ~isfield(plot_settings.SymbolSize.locks,'sep')
  plot_settings.SymbolSize.locks.sep = 1;
end
if ~isfield(plot_settings.SymbolSize.locks,'snf')
  plot_settings.SymbolSize.locks.snf = 1;
end
if ~isfield(plot_settings.SymbolSize.locks,'fluxexp')
  plot_settings.SymbolSize.locks.fluxexp = 1;
end
if ~isfield(plot_settings.SymbolSize.locks,'bdef')
  plot_settings.SymbolSize.locks.bdef = 1;
end

if ~isfield(plot_settings.SymbolSize,'limits')
  plot_settings.SymbolSize.limits = [];
end
if ~isfield(plot_settings.SymbolSize.limits,'x')
  plot_settings.SymbolSize.limits.x = 1;
end
if ~isfield(plot_settings.SymbolSize.limits,'nox')
  plot_settings.SymbolSize.limits.nox = 1;
end
if ~isfield(plot_settings.SymbolSize.limits,'sep')
  plot_settings.SymbolSize.limits.sep = 1;
end
if ~isfield(plot_settings.SymbolSize.limits,'snf')
  plot_settings.SymbolSize.limits.snf = 1;
end
if ~isfield(plot_settings.SymbolSize.limits,'bdef')
  plot_settings.SymbolSize.limits.bdef = 1;
end

if ~isfield(plot_settings.SymbolSize,'targets')
  plot_settings.SymbolSize.targets = [];
end
if ~isfield(plot_settings.SymbolSize.targets,'x')
  plot_settings.SymbolSize.targets.x = 1;
end
if ~isfield(plot_settings.SymbolSize.targets,'sep')
  plot_settings.SymbolSize.targets.sep = 1;
end
if ~isfield(plot_settings.SymbolSize.targets,'snf')
  plot_settings.SymbolSize.targets.snf = 1;
end
if ~isfield(plot_settings.SymbolSize.targets,'fluxexp')
  plot_settings.SymbolSize.targets.fluxexp = 1;
end
if ~isfield(plot_settings.SymbolSize.targets,'bdef')
  plot_settings.SymbolSize.targets.bdef = 1;
end

if ~isfield(plot_settings.SymbolSize,'wrongway')
  plot_settings.SymbolSize.wrongway = 1;
end
if ~isfield(plot_settings.SymbolSize,'bdef')
  plot_settings.SymbolSize.bdef = 1;
end
if ~isfield(plot_settings.SymbolSize,'x')
  plot_settings.SymbolSize.x = 1;
end

% Maximum number of x-points to show and report positions for
nxptmax = plot_settings.nxptmax;

% Number of flux surface contours from axis to boundary
nflux = plot_settings.nflux;

% Number of SOL contours
nSOL = plot_settings.SOL.n;

% Create new plot window, if needed or desired
if ~exist('ha','var')                       | ...
   isempty(ha)                              | ...
   ~ishandle(ha.figure)                     | ...
   ~all(ishandle(ha.hs))                    | ...
   iteration == 1 & fig < 0 & ~done | ...
   fig > 0 & ha.figure ~= fig

  newhandles = true;
  if fig > 0
    if ishandle(fig)
      close(fig)
    end
    ha.figure = figure(fig);
  else
    ha.figure = figure;
  end
  fpos = get(ha.figure,'position');
  fpos(4) = fpos(3)/2;
  fpos(3:4) = [1000 500];
  set(ha.figure,'position',fpos)
  if isa(ha.figure,'matlab.ui.Figure')
    ha.figure = ha.figure.Number;
  end
  zoom on
 
  ha.h1 = []; % handles to objects with font size fs1
  ha.h2 = []; % handles to objects with font size fs2
  ha.f1 = 0.8/60;  % fs1/window width in pixels
  ha.f2 = 0.8/100; % fs2/window width in pixels
  ha.fs = 0.8./[80 80 80 80 80 80];  % fs1/window width in pixels
  
  % CONDUCTORS
  ha.hs(1) = subplot(2,4,1,'Position',[0.08 0.6 0.2 0.3]);
  hold on
  ha.conductors.limits.ic = [];
  if exist('limits','var') & isfield(limits,'ic')
    for i = 1:nc
      ha.conductors.limits.ic(i) = plot(i,nan,'yo--', ...
        'MarkerSize',18,'MarkerFaceColor','yellow','LineWidth',3);
    end
    plot(limits.ic(:,1),'b--')
    plot(limits.ic(:,2),'r--')
    plot(limits.ic(:,1),'bx','MarkerSize',12,'LineWidth',3)
    plot(limits.ic(:,2),'rx','MarkerSize',12,'LineWidth',3)
  end
  ha.conductors.dashed = plot([ic; iv],'g--');
  ha.conductors.rings =  plot([ic; iv],'go', ...
    'MarkerSize',4,'MarkerFaceColor','green','LineWidth',3);
  ha.h1(end+1) = title('Conductors');
  grid on
 
  % pprime, ffprim
  ha.hs(2) = subplot(2,4,2,'Position',[0.33 0.6 0.2 0.3]);
  hold(ha.hs(2),'on')
  ha.profiles.pprime = plot(ha.hs(2),c.psibarp,c.psibarp,'blue','LineWidth',3);
  ha.profiles.ffprim = plot(ha.hs(2),c.psibarp,c.psibarp,'red', 'LineWidth',3);
  ha.h1(end+1) = title(ha.hs(2),'{\color{blue}Rp''}, {\color{red}ff''/\mu_0/R}');
  grid(ha.hs(2),'on')
  
  % LIST OF SCALAR VALUES
  ha.hs(3) = subplot(2,4,3,'Position',[0.55 0.1 0.175 0.8]);
  plot([0 1],[0 1],'w')
  axis(ha.hs(3),[0 1 0 1])
  set(ha.hs(3),'Xtick',[])
  set(ha.hs(3),'Ytick',[])
  ha.h1(end+1) = text(0,0,'I_p');
  ha.scalars.cpasma.name = ha.h1(end);
  ha.h1(end+1) = text(0,0,'0');
  ha.scalars.cpasma.value = ha.h1(end);
  ha.h1(end+1) = text(0,0,'l_i');
  ha.scalars.li.name = ha.h1(end);
  ha.h1(end+1) = text(0,0,'0');
  ha.scalars.li.value = ha.h1(end);
  ha.h1(end+1) = text(0,0,'\beta_p');
  ha.scalars.betap.name = ha.h1(end);
  ha.h1(end+1) = text(0,0,'0');
  ha.scalars.betap.value = ha.h1(end);
  ha.h1(end+1) = text(0,0,'\beta_n');
  ha.scalars.betan.name = ha.h1(end);
  ha.h1(end+1) = text(0,0,'0');
  ha.scalars.betan.value = ha.h1(end);
  ha.h1(end+1) = text(0,0,'\psi_a');
  ha.scalars.psimag.name = ha.h1(end);
  ha.h1(end+1) = text(0,0,'0');
  ha.scalars.psimag.value = ha.h1(end);
  ha.h1(end+1) = text(0,0,'\psi_b');
  ha.scalars.psibry.name = ha.h1(end);
  ha.h1(end+1) = text(0,0,'0');
  ha.scalars.psibry.value = ha.h1(end);
  ha.h1(end+1) = text(0,0,'\psi_p');
  ha.scalars.psipla.name = ha.h1(end);
  ha.h1(end+1) = text(0,0,'0');
  ha.scalars.psipla.value = ha.h1(end);
  ha.scalars.rx.name = [];
  ha.scalars.zx.name = [];
  ha.scalars.rx.value = [];
  ha.scalars.zx.value = [];
  for i = 1:nxptmax
    ha.h1(end+1) = text(0,0,' ');
    ha.scalars.rx.name(i) = ha.h1(end);
    ha.h1(end+1) = text(0,0,'0');
    ha.scalars.rx.value(i) = ha.h1(end);
    ha.h1(end+1) = text(0,0,' ');
    ha.scalars.zx.name(i) = ha.h1(end);
    ha.h1(end+1) = text(0,0,'0');
    ha.scalars.zx.value(i) = ha.h1(end);
  end
  ha.h1(end+1) = title(' ');
  ha.scalars.title = ha.h1(end);
  
  % MACHINE GEOMETRY
  ha.hs(4) = subplot(2,4,4,'Position',[0.765 0.1 0.225 0.8]);
  hold on
  patch(c.Rlim,c.Zlim,[0.8 0.8 0.8])
  if exist('showgrid','var') & showgrid
    plot(rg(1)+ones(2,1)*linspace(-0.5,nr-0.5,nr+1)*dr,...
       [zg(1)-dz/2; zg(nz)+dz/2]*ones(1,nr+1),'--','color',[.9 .9 .9]);
    plot([rg(1)-dr/2; rg(nr)+dr/2]*ones(1,nz+1),...
       zg(1)+ones(2,1)*linspace(-0.5,nz-0.5,nz+1)*dz,'--','color',[.9 .9 .9]);
    plot(rgg,zgg,'+','color',[.85 .85 .85]);
  end
  plot(c.Rlim,c.Zlim,'k','LineWidth',4)
  ha.tokamak.patch = patch(nan,nan,[1 1 0.9]);
  goldenrod = [218 165 32]/255;
  ha.tokamak.wrongway(1) = plot(nan,nan,'ro','MarkerFaceColor','red');
  ha.tokamak.wrongway(2) = plot(nan,nan,'ws','MarkerFaceColor','white');
  ha.tokamak.targets.snf = plot(nan,nan,'w*');
  ha.tokamak.targets.fluxexp = plot(nan,nan,'kp','MarkerFaceColor','yellow');
  ha.tokamak.locks.snf = plot(nan,nan,'*','Color',goldenrod);
  ha.tokamak.locks.fluxexp = plot(nan,nan,'kp','MarkerFaceColor',goldenrod);
  ha.tokamak.targets.x = plot(nan,nan,'yx');
  ha.tokamak.limits.x = plot(nan,nan,'--y');
  ha.tokamak.limits.nox = plot(nan,nan,'--r');
  ha.tokamak.locks.x = plot(nan,nan,'x','Color',goldenrod);
  ha.tokamak.targets.sep = plot(nan,nan,'mo');
  ha.tokamak.limits.sep(1) = plot(nan,nan,'bo');
  ha.tokamak.limits.sep(2) = plot(nan,nan,'ro');
  ha.tokamak.locks.sep = plot(nan,nan,'o','Color',goldenrod);
  ha.tokamak.bbbs = plot(nan,nan,'LineWidth',3,'color',[0,0, 0.5625]);
  ha.tokamak.targets.bdef(1) = plot(nan,nan,'y+');
  ha.tokamak.targets.bdef(2) = plot(nan,nan,'yo');
  ha.tokamak.locks.bdef(1) = plot(nan,nan,'+','Color',goldenrod);
  ha.tokamak.locks.bdef(2) = plot(nan,nan,'o','Color',goldenrod);
  ha.tokamak.bdef(1) = plot(nan,nan,'yx');
  ha.tokamak.bdef(2) = plot(nan,nan,'rx');
  map = colormap;
  n = size(map,1);
  ha.tokamak.contours = [];
  for i = 1:nflux-1
    ha.tokamak.contours(i) = plot(nan,nan,'Color',map(round(i*n/nflux),:));
  end
  ha.tokamak.sol = [];
  for i = 1:nSOL
    ha.tokamak.sol(i) = plot(nan,nan,'Color',map(round(i*n/nSOL),:));
  end
  ha.tokamak.x = [];
  for i = 1:nxptmax
    ha.tokamak.x(i) = plot(nan,nan,'bx');
  end
  axis('image')
  axis([c.rg(1) c.rg(c.nr) c.zg(1) c.zg(c.nz)])
  ha.h1(end+1) = title(upper(c.tokamak));
  
  % ERROR VECTOR
  ha.hs(5) = subplot(2,4,5,'Position',[0.08 0.1 0.2 0.3]);
  hold on
  ha.ev.value = plot(ev,'go','LineWidth',3);
  ha.h1(end+1) = title('Error vector');
  ha.h1(end+1) = text(1,0,' ','Color','red');
  ha.ev.max.name = ha.h1(end);
  ha.ev.max.value = plot(nan,nan,'ro','LineWidth',3);
  ha.h1(end+1) = text(1,0,' ','Color','blue');
  ha.ev.min.name = ha.h1(end);
  ha.ev.min.value = plot(nan,nan,'bo','LineWidth',3);
  for i = 1:28
    text(0,0,'')
  end
  
  % FLUX ERROR
  ha.hs(6) = subplot(2,4,6,'Position',[0.33 0.1 0.2 0.3]);
  ha.fluxerror = plot(nan(c.nz,c.nr));
  ha.h1(end+1) = title('\psi error');
  grid on

  % RESIZE FUNCTION
  set(ha.figure,'UserData',ha)
  set(ha.figure,'ResizeFcn','gs_resize_fonts(get(gcbo,''UserData''))')
  gs_resize_fonts(ha)
  
  % In matlab2016b the subplots need to be put back into their positions
  set(ha.hs(1),'Position',[0.08 0.6 0.2 0.3])
  set(ha.hs(2),'Position',[0.33 0.6 0.2 0.3])
  set(ha.hs(3),'Position',[0.55 0.1 0.175 0.8])
  set(ha.hs(4),'Position',[0.765 0.1 0.225 0.8])
  set(ha.hs(5),'Position',[0.08 0.1 0.2 0.3])
  set(ha.hs(6),'Position',[0.33 0.1 0.2 0.3])
  
  % BUTTONS
  ha.toolbar = findall(ha.figure,'Type','uitoolbar');
  % STOP
  cdata = zeros(20,20,3);
  for j = 1:20
    for k = 1:20
      cdata(j,k,1:3) = 0.9;
      if (j-10.5)^2+(k-10.5)^2<100
	cdata(j,k,1) = 1;
	cdata(j,k,2:3) = 0;
      end
      if j>5 & j<16 & k>5 & k<16
	cdata(j,k,2:3) = 1;
      end
    end
  end
  ha.button.stop = uitoggletool(ha.toolbar, ...
    'CData',cdata, ...
    'Userdata',ha, ...
    'Enable','on', ...
    'TooltipString','Stop iterating', ...
    'Separator','on', ...
    'HandleVisibility','off');
    % If window recreated between iterations, onCallback must be reactived
    set(ha.button.stop,'OnCallback', [...
	'set(getfield(get(gcbo,''UserData''),''figure''),''Name'',', ...
	'[get(getfield(get(gcbo,''UserData''),''figure''),''Name'') ', ...
	''' - stopping'']),drawnow'])
  % HELP
  cdata = zeros(20,20,3);
  for j = 1:20
    for k = 1:20
      cdata(j,k,1:3) = 0.9;
      if (j-10.5)^2+(k-10.5)^2<100
	cdata(j,k,1) = 0;
	cdata(j,k,2) = 0.9;
	cdata(j,k,3) = 0.9;
      end
    end
  end
  cdata(3,9:12,1:3) = 1;
  cdata(4,8:13,1:3) = 1;
  cdata(5,[7:9 12:14],1:3) = 1;
  cdata(6,[6:8 13:15],1:3) = 1;
  cdata(7,[6:7 14:15],1:3) = 1;
  cdata(8,13:15,1:3) = 1;
  cdata(9,12:14,1:3) = 1;
  cdata(10,11:13,1:3) = 1;
  cdata(11,10:12,1:3) = 1;
  cdata(12:13,10:11,1:3) = 1;
  cdata(16:18,10:11,1:3) = 1;
  ha.button.help = uitoggletool(ha.toolbar, ...
    'CData',cdata, ...
    'Userdata',ha, ...
    'Enable','on', ...
    'TooltipString','Explain figure', ...
    'Separator','on', ...
    'HandleVisibility','off', ...
    'OnCallback','set(gcbo,''State'',''off'');gsdesign_figure_help');
else
  newhandles = false;
end

if newhandles | iteration == 1
  if ishandle(ha.tokamak.wrongway(1))
    ms = max(0.001,abs(plot_settings.SymbolSize.wrongway*14));
    lw = max(0.001,abs(plot_settings.SymbolSize.wrongway*1));
    set(ha.tokamak.wrongway(1),'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.wrongway(2))
    ms = max(0.001,abs(plot_settings.SymbolSize.wrongway*4));
    set(ha.tokamak.wrongway(2),'MarkerSize',ms)
  end
  if ishandle(ha.tokamak.targets.snf)
    ms = max(0.001,abs(plot_settings.SymbolSize.targets.snf*24));
    lw = max(0.001,abs(plot_settings.SymbolSize.targets.snf*4));
    set(ha.tokamak.targets.snf,'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.targets.fluxexp)
    ms = max(0.001,abs(plot_settings.SymbolSize.targets.fluxexp*18));
    lw = max(0.001,abs(plot_settings.SymbolSize.targets.fluxexp*3));
    set(ha.tokamak.targets.fluxexp,'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.locks.snf)
    ms = max(0.001,abs(plot_settings.SymbolSize.locks.snf*24));
    lw = max(0.001,abs(plot_settings.SymbolSize.locks.snf*4));
    set(ha.tokamak.locks.snf,'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.locks.fluxexp)
    ms = max(0.001,abs(plot_settings.SymbolSize.locks.fluxexp*18));
    lw = max(0.001,abs(plot_settings.SymbolSize.locks.fluxexp*3));
    set(ha.tokamak.locks.fluxexp,'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.targets.x)
    ms = max(0.001,abs(plot_settings.SymbolSize.targets.x*18));
    lw = max(0.001,abs(plot_settings.SymbolSize.targets.x*3));
    set(ha.tokamak.targets.x,'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.limits.x)
    lw = max(0.001,abs(plot_settings.SymbolSize.limits.x*3));
    set(ha.tokamak.limits.x,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.limits.nox)
    lw = max(0.001,abs(plot_settings.SymbolSize.limits.nox*3));
    set(ha.tokamak.limits.nox,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.locks.x)
    ms = max(0.001,abs(plot_settings.SymbolSize.locks.x*18));
    lw = max(0.001,abs(plot_settings.SymbolSize.locks.x*3));
    set(ha.tokamak.locks.x,'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.targets.sep)
    ms = max(0.001,abs(plot_settings.SymbolSize.targets.sep*6));
    lw = max(0.001,abs(plot_settings.SymbolSize.targets.sep*2));
    set(ha.tokamak.targets.sep,'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.limits.sep(1))
    ms = max(0.001,abs(plot_settings.SymbolSize.limits.sep*6));
    lw = max(0.001,abs(plot_settings.SymbolSize.limits.sep*2));
    set(ha.tokamak.limits.sep(1),'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.limits.sep(2))
    ms = max(0.001,abs(plot_settings.SymbolSize.limits.sep*6));
    lw = max(0.001,abs(plot_settings.SymbolSize.limits.sep*2));
    set(ha.tokamak.limits.sep(2),'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.locks.sep)
    ms = max(0.001,abs(plot_settings.SymbolSize.locks.sep*11));
    lw = max(0.001,abs(plot_settings.SymbolSize.locks.sep*4));
    set(ha.tokamak.locks.sep,'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.targets.bdef(1))
    ms = max(0.001,abs(plot_settings.SymbolSize.targets.bdef*12));
    lw = max(0.001,abs(plot_settings.SymbolSize.targets.bdef*2));
    set(ha.tokamak.targets.bdef(1),'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.targets.bdef(2))
    ms = max(0.001,abs(plot_settings.SymbolSize.targets.bdef*18));
    lw = max(0.001,abs(plot_settings.SymbolSize.targets.bdef*3));
    set(ha.tokamak.targets.bdef(2),'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.locks.bdef(1))
    ms = max(0.001,abs(plot_settings.SymbolSize.locks.bdef*12));
    lw = max(0.001,abs(plot_settings.SymbolSize.locks.bdef*2));
    set(ha.tokamak.locks.bdef(1),'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.locks.bdef(2))
    ms = max(0.001,abs(plot_settings.SymbolSize.locks.bdef*18));
    lw = max(0.001,abs(plot_settings.SymbolSize.locks.bdef*3));
    set(ha.tokamak.locks.bdef(2),'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.bdef(1))
    ms = max(0.001,abs(plot_settings.SymbolSize.bdef*24));
    lw = max(0.001,abs(plot_settings.SymbolSize.bdef*5));
    set(ha.tokamak.bdef(1),'MarkerSize',ms,'LineWidth',lw)
  end
  if ishandle(ha.tokamak.bdef(2))
    ms = max(0.001,abs(plot_settings.SymbolSize.bdef*21));
    lw = max(0.001,abs(plot_settings.SymbolSize.bdef*3));
    set(ha.tokamak.bdef(2),'MarkerSize',ms,'LineWidth',lw)
  end
  ms = max(0.001,abs(plot_settings.SymbolSize.x*12));
  lw = max(0.001,abs(plot_settings.SymbolSize.x*2));
  for i = 1:nxptmax
    if ishandle(ha.tokamak.x(i))
      set(ha.tokamak.x(i),'MarkerSize',ms,'LineWidth',lw)
    end
  end
end

if iteration == 1
  set(ha.button.stop,'State','off')
  set(ha.button.stop,'OnCallback', [...
      'set(getfield(get(gcbo,''UserData''),''figure''),''Name'',', ...
      '[get(getfield(get(gcbo,''UserData''),''figure''),''Name'') ', ...
      ''' - stopping'']),drawnow'])
end

if exist('figure_message','var')
  set(ha.figure,'Name',figure_message)
elseif exist('iteration','var')
  str = ['Iteration ' num2str(iteration)];
  set(ha.figure,'Name',str)
end

if done
  set(ha.button.stop,'OnCallback','')
else
  if strcmp(get(ha.button.stop,'State'),'on')
    str = get(ha.figure,'Name');
    str = [str ' - stopping'];
    set(ha.figure,'Name',str)
  end
end

% UPDATE CONDUCTOR CURRENTS
if exist('limits','var') & isfield(limits,'ic')
  for i = 1:nc
    if ishandle(ha.conductors.limits.ic(i))
      if ic(i) < limits.ic(i,1)+1e-3 | ic(i) > limits.ic(i,2)-1e-3
        set(ha.conductors.limits.ic(i),'YData',ic(i))
      else
        set(ha.conductors.limits.ic(i),'YData',nan)
      end
    end
  end
end
if ishandle(ha.conductors.dashed)
  set(ha.conductors.dashed,'XData', 1:nc+nv);
  set(ha.conductors.dashed,'YData', [ic; iv]);
end
if ishandle(ha.conductors.rings)
  set(ha.conductors.rings,'XData', 1:nc+nv);
  set(ha.conductors.rings,'YData', [ic; iv]);
end
if all(iv == 0)
  xmax = nc+1/2;
  ymin = 1.2*min(ic)-0.2*max(ic);
  ymax = 1.2*max(ic)-0.2*min(ic);
else
  xmax = nc+nv+1/2;
  ymin = 1.2*min([ic; iv])-0.2*max([ic; iv]);
  ymax = 1.2*max([ic; iv])-0.2*min([ic; iv]);
end
if ymin == ymax
  ymin = ymin-0.5;
  ymax = ymax+0.5;
end
axis(ha.hs(1),[0.5 xmax ymin ymax])

% UPDATE pprime, ffprim
dumnr = e.rmaxis*e.pprime;
if min(dumnr)*max(dumnr) < 0 % min must be negative and max positive
  if min(dumnr)/max(dumnr) > -1e-9
    dumnr(dumnr < 0) = 0;
  elseif max(dumnr)/min(dumnr) > -1e-9
    dumnr(dumnr > 0) = 0;
  end
end
if ishandle(ha.profiles.pprime)
  set(ha.profiles.pprime,'XData',c.psibarp)
  set(ha.profiles.pprime,'YData',dumnr)
end
dumnr = e.ffprim/mu0/e.rmaxis;
if min(dumnr)*max(dumnr) < 0 % min must be negative and max positive
  if min(dumnr)/max(dumnr) > -1e-6
    dumnr(dumnr < 0) = 0;
  elseif max(dumnr)/min(dumnr) > -1e-6
    dumnr(dumnr > 0) = 0;
  end
end
if ishandle(ha.profiles.ffprim)
  set(ha.profiles.ffprim,'XData',c.psibarp)
  set(ha.profiles.ffprim,'YData',dumnr)
end

% UPDATE LIST OF SCALAR VALUES
ix = 0; % index to x-point that defines boundary
rx = nan(1,nxptmax);
zx = nan(1,nxptmax);
nxpt = 0; % number of x-points
for i = 1:b.nulls.count
  if b.nulls.type(i) == 'X'
    nxpt = nxpt+1;
    if b.nulls.r(i) == b.bd.r(1) & b.nulls.z(i) == b.bd.z(1)
      ix = nxpt;
    end
    rx(nxpt) = c.rg(1) + (b.nulls.r(i)-1)*c.dr;
    zx(nxpt) = c.zg(1) + (b.nulls.z(i)-1)*c.dz;    
  end
  if nxpt == nxptmax
    break
  end
end
nrow = 7+2*nxpt; % Total number of rows in the list
set(ha.scalars.cpasma.name,'Position', [0.1,0.985-1/nrow*0.9])
set(ha.scalars.cpasma.value,'Position', [0.4,0.990-1/nrow*0.9])
set(ha.scalars.cpasma.value,'String', sprintf('%8d',round(e.cpasma)))
set(ha.scalars.li.name,'Position', [0.1,0.985-2/nrow*0.9])
set(ha.scalars.li.value,'Position', [0.4,0.990-2/nrow*0.9])
set(ha.scalars.li.value,'String', sprintf('%7.5f',e.li))
set(ha.scalars.betap.name,'Position', [0.1,0.985-3/nrow*0.9])
set(ha.scalars.betap.value,'Position', [0.4,0.990-3/nrow*0.9])
set(ha.scalars.betap.value,'String', sprintf('%7.5f',e.betap))
set(ha.scalars.betan.name,'Position', [0.1,0.985-4/nrow*0.9])
set(ha.scalars.betan.value,'Position', [0.4,0.990-4/nrow*0.9])
set(ha.scalars.betan.value,'String', sprintf('%7.5f',e.betan))
set(ha.scalars.psimag.name,'Position', [0.1,0.985-5/nrow*0.9])
set(ha.scalars.psimag.value,'Position', [0.4,0.990-5/nrow*0.9])
set(ha.scalars.psimag.value,'String', sprintf('%7.5f',e.psimag))
set(ha.scalars.psibry.name,'Position', [0.1,0.985-6/nrow*0.9])
set(ha.scalars.psibry.value,'Position', [0.4,0.990-6/nrow*0.9])
set(ha.scalars.psibry.value,'String', sprintf('%7.5f',e.psibry))
set(ha.scalars.psipla.name,'Position', [0.1,0.985-7/nrow*0.9])
set(ha.scalars.psipla.value,'Position', [0.4,0.990-7/nrow*0.9])
set(ha.scalars.psipla.value,'String', sprintf('%7.5f',e.psipla))
for i = 1:nxptmax
  if i <= nxpt
    set(ha.scalars.rx.name(i),'String', ['R_x_' num2str(i)])
    set(ha.scalars.zx.name(i),'String', ['Z_x_' num2str(i)])
    set(ha.scalars.rx.name(i),'Position', [0.1,0.985-(6+2*i)/nrow*0.9])
    set(ha.scalars.zx.name(i),'Position', [0.1,0.985-(7+2*i)/nrow*0.9])
    set(ha.scalars.rx.value(i),'String', sprintf('%7.5f',rx(i)))
    set(ha.scalars.zx.value(i),'String', sprintf('%7.5f',zx(i)))
    set(ha.scalars.rx.value(i),'Position', [0.4,0.990-(6+2*i)/nrow*0.9])
    set(ha.scalars.zx.value(i),'Position', [0.4,0.990-(7+2*i)/nrow*0.9])
    if i == ix
      set(ha.scalars.rx.name(i),'Color', 'red')
      set(ha.scalars.zx.name(i),'Color', 'red')
      set(ha.scalars.rx.value(i),'Color', 'red')
      set(ha.scalars.zx.value(i),'Color', 'red')
    else
      set(ha.scalars.rx.name(i),'Color', 'blue')
      set(ha.scalars.zx.name(i),'Color', 'blue')
      set(ha.scalars.rx.value(i),'Color', 'blue')
      set(ha.scalars.zx.value(i),'Color', 'blue')
    end
  else
    set(ha.scalars.rx.name(i),'String', '     ')
    set(ha.scalars.zx.name(i),'String', '     ')
    set(ha.scalars.rx.value(i),'String', '       ')
    set(ha.scalars.zx.value(i),'String', '       ')
  end
end
if ishandle(ha.scalars.title)
  if exist('time','var') & ~isempty(time) & ~isnan(time)
    set(ha.scalars.title,'String',['time = ' num2str(time)]);
  else
    set(ha.scalars.title,'String',' ');
  end
end

% UPDATE PLASMA GEOMETRY
if all(ishandle(ha.tokamak.wrongway))
  set(ha.tokamak.wrongway,'XData',wrongway.R)
  set(ha.tokamak.wrongway,'YData',wrongway.Z)
end
if ishandle(ha.tokamak.patch)
  set(ha.tokamak.patch,'XData',b.R(~isnan(b.R)))
  set(ha.tokamak.patch,'YData',b.Z(~isnan(b.R)))
end
plasma = b.plasma;
if plasma
  psia = e.psimag;
  psib = e.psibry;
else
  igg = inpolygon(c.rgg,c.zgg,c.Rlim,c.Zlim);
  psia = max(e.psizr(igg(:)));
  psib = min(e.psizr(igg(:)));
end
if exist('nhires','var') & nhires > 1
  nzhr = 1+nhires*(c.nz-1);
  nrhr = 1+nhires*(c.nr-1);
  rghr = linspace(c.rg(1),c.rg(c.nr),nrhr);
  zghr = linspace(c.zg(1),c.zg(c.nz),nzhr)';
  px = gs_interp2(c.rg,c.zg,e.psizr,ones(nzhr,1)*rghr, zghr*ones(1,nrhr));
end
for i = 1:nflux-1
  if exist('nhires','var') & nhires > 1
    ccc = contourc(rghr,zghr,1-(px-psia)/(psib-psia),i/nflux+[0 0]);
  else
    ccc = contourc(c.rg,c.zg,1-(e.psizr-psia)/(psib-psia),i/nflux+[0 0]);
  end
  j = 1;
  n = size(ccc,2);
  while j < n
    k = ccc(2,j);
    ccc(:,j) = nan;
    j = j+k+1;
  end
  if ishandle(ha.tokamak.contours(i))
    set(ha.tokamak.contours(i),'XData',ccc(1,:))
    set(ha.tokamak.contours(i),'YData',ccc(2,:))
  end
end
if length(ha.tokamak.sol) ~= nSOL | any(~ishandle(ha.tokamak.sol))
  delete(ha.tokamak.sol(ishandle(ha.tokamak.sol)))
  n = size(map,1);
  ha.tokamak.sol = [];
  for i = 1:nSOL
    ha.tokamak.sol(i) = plot(ha.hs(4),0,0,'Color',map(round(i*n/nSOL),:));
  end
end
set(ha.tokamak.sol,'XData',nan)
set(ha.tokamak.sol,'YData',nan)
if plasma & nSOL > 0
  rsol = b.R8(1) + plot_settings.SOL.d*(1:plot_settings.SOL.n);
  psisol = gs_interp2(rg,zg,e.psizr,rsol,b.Z8(1));
  for i = 1:nSOL
    if exist('nhires','var') & nhires > 1
      ccc = contourc(rghr,zghr,1-(px-psia)/(psib-psia), ...
        1-(psisol(i)-psia)/(psib-psia)+[0 0]);
    else
      ccc = contourc(rg,zg,1-(e.psizr-psia)/(psib-psia), ...
        1-(psisol(i)-psia)/(psib-psia)+[0 0]);
    end
    j = 1;
    n = size(ccc,2);
    while j < n
      k = ccc(2,j);
      ccc(:,j) = nan;
      j = j+k+1;
    end
    if ishandle(ha.tokamak.sol(i))
      set(ha.tokamak.sol(i),'XData',ccc(1,:))
      set(ha.tokamak.sol(i),'YData',ccc(2,:))
    end
  end
end
for i = nxptmax:-1:1
  if ishandle(ha.tokamak.x(i))
    if i ~= ix
      set(ha.tokamak.x(i),'XData',rx(i))
      set(ha.tokamak.x(i),'YData',zx(i))
    else
      set(ha.tokamak.x(i),'XData',nan)
      set(ha.tokamak.x(i),'YData',nan)
    end
  end
end
if exist('locks','var')
  if isfield(locks,'rsnf') & isfield(locks,'zsnf')
    if ishandle(ha.tokamak.locks.snf)
      set(ha.tokamak.locks.snf,'XData',locks.rsnf)
      set(ha.tokamak.locks.snf,'YData',locks.zsnf)
    end
  end
  if isfield(locks,'rfluxexp') & isfield(locks,'zfluxexp')
    if ishandle(ha.tokamak.locks.fluxexp)
      set(ha.tokamak.locks.fluxexp,'XData',locks.rfluxexp)
      set(ha.tokamak.locks.fluxexp,'YData',locks.zfluxexp)
    end
  end
  if isfield(locks,'rx') & isfield(locks,'zx')
    if ishandle(ha.tokamak.locks.x)
      set(ha.tokamak.locks.x,'XData',locks.rx)
      set(ha.tokamak.locks.x,'YData',locks.zx)
    end
  end
  if isfield(locks,'rsep') & isfield(locks,'zsep')
    if ishandle(ha.tokamak.locks.sep)
      set(ha.tokamak.locks.sep,'XData',locks.rsep)
      set(ha.tokamak.locks.sep,'YData',locks.zsep)
    end
  end
  if isfield(locks,'rbdef') & isfield(locks,'zbdef')
    if all(ishandle(ha.tokamak.locks.bdef))
      set(ha.tokamak.locks.bdef,'XData',locks.rbdef)
      set(ha.tokamak.locks.bdef,'YData',locks.zbdef)
    end
  end
end
if ishandle(ha.tokamak.targets.snf)
  set(ha.tokamak.targets.snf,'XData',nan)
  set(ha.tokamak.targets.snf,'YData',nan)
end
if ishandle(ha.tokamak.targets.fluxexp)
  set(ha.tokamak.targets.fluxexp,'XData',nan)
  set(ha.tokamak.targets.fluxexp,'YData',nan)
end
if ishandle(ha.tokamak.targets.x)
  set(ha.tokamak.targets.x,'XData',nan)
  set(ha.tokamak.targets.x,'YData',nan)
end
if ishandle(ha.tokamak.targets.sep)
  set(ha.tokamak.targets.sep,'XData',nan)
  set(ha.tokamak.targets.sep,'YData',nan)
end
if all(ishandle(ha.tokamak.targets.bdef))
  set(ha.tokamak.targets.bdef,'XData',nan)
  set(ha.tokamak.targets.bdef,'YData',nan)
end
if exist('targets','var')
  if isfield(targets,'rsnf') & isfield(targets,'zsnf')
    if ishandle(ha.tokamak.targets.snf)
      set(ha.tokamak.targets.snf,'XData',targets.rsnf)
      set(ha.tokamak.targets.snf,'YData',targets.zsnf)
    end
  end
  if isfield(targets,'rfluxexp') & isfield(targets,'zfluxexp')
    if ishandle(ha.tokamak.targets.fluxexp)
      set(ha.tokamak.targets.fluxexp,'XData',targets.rfluxexp)
      set(ha.tokamak.targets.fluxexp,'YData',targets.zfluxexp)
    end
  end
  if isfield(targets,'rx') & isfield(targets,'zx')
    if ishandle(ha.tokamak.targets.x)
      set(ha.tokamak.targets.x,'XData',targets.rx)
      set(ha.tokamak.targets.x,'YData',targets.zx)
    end
  end
  if isfield(targets,'rsep') & isfield(targets,'zsep')
    if ishandle(ha.tokamak.targets.sep)
      set(ha.tokamak.targets.sep,'XData',targets.rsep)
      set(ha.tokamak.targets.sep,'YData',targets.zsep)
    end
  end
  if isfield(targets,'rbdef') & isfield(targets,'zbdef')
    if all(ishandle(ha.tokamak.targets.bdef))
      set(ha.tokamak.targets.bdef,'XData',targets.rbdef)
      set(ha.tokamak.targets.bdef,'YData',targets.zbdef)
    end
  end
end
if ishandle(ha.tokamak.limits.sep(1))
  set(ha.tokamak.limits.sep(1),'XData',nan)
  set(ha.tokamak.limits.sep(1),'YData',nan)
end
if ishandle(ha.tokamak.limits.sep(2))
  set(ha.tokamak.limits.sep(2),'XData',nan)
  set(ha.tokamak.limits.sep(2),'YData',nan)
end
if ishandle(ha.tokamak.limits.x)
  set(ha.tokamak.limits.x,'XData',nan)
  set(ha.tokamak.limits.x,'YData',nan)
end
if ishandle(ha.tokamak.limits.nox)
  set(ha.tokamak.limits.nox,'XData',nan)
  set(ha.tokamak.limits.nox,'YData',nan)
end
if exist('limits','var')
  if isfield(limits,'rsep') & isfield(limits,'zsep')
    if ishandle(ha.tokamak.limits.sep(1))
      set(ha.tokamak.limits.sep(1),'XData',limits.rsep(:,1))
      set(ha.tokamak.limits.sep(1),'YData',limits.zsep(:,1))
    end
    if ishandle(ha.tokamak.limits.sep(2))
      set(ha.tokamak.limits.sep(2),'XData',limits.rsep(:,2))
      set(ha.tokamak.limits.sep(2),'YData',limits.zsep(:,2))
    end
  end
  if isfield(limits,'rx') & isfield(limits,'zx')
    if ishandle(ha.tokamak.limits.x)
      rs = limits.rx';
      zs = limits.zx';
      set(ha.tokamak.limits.x,'XData',rs(:))
      set(ha.tokamak.limits.x,'YData',zs(:))
    end
  end
  if isfield(limits,'rnox') & isfield(limits,'znox')
    if ishandle(ha.tokamak.limits.nox)
      rs = limits.rnox';
      zs = limits.znox';
      set(ha.tokamak.limits.nox,'XData',rs(:))
      set(ha.tokamak.limits.nox,'YData',zs(:))
    end
  end
end
if ishandle(ha.tokamak.bbbs)
  set(ha.tokamak.bbbs,'XData',b.R)
  set(ha.tokamak.bbbs,'YData',b.Z)
end
if all(ishandle(ha.tokamak.bdef))
  set(ha.tokamak.bdef,'XData',e.rbdef)
  set(ha.tokamak.bdef,'YData',e.zbdef)
end

% UPDATE ERROR VECTOR
set(ha.ev.value,'YData',ev)
if ~isempty(ev)
  [zmin,kmin] = min(ev);
  [zmax,kmax] = max(ev);
  if zmax > zmin
    label_min_and_max_error = true;
  else
    label_min_and_max_error = false;
    if zmin == 0
      zmin = zmin-0.5;
      zmax = zmax+0.5;
    else
      zmin = -2*abs(zmin);
      zmax = +2*abs(zmax);
    end
  end
  set(ha.ev.max.value,'XData',kmax)
  set(ha.ev.max.value,'YData',zmax)
  set(ha.ev.max.name,'String',['  ' deblank(strev(kmax,:))])
  set(ha.ev.max.name,'Position',[kmax zmax*1.1-zmin*0.1])
  set(ha.ev.min.value,'XData',kmin)
  set(ha.ev.min.value,'YData',zmin)
  set(ha.ev.min.name,'String',['  ' deblank(strev(kmin,:))])
  set(ha.ev.min.name,'Position',[kmin zmin*1.1-zmax*0.1])
  if ~label_min_and_max_error
    set(ha.ev.max.value,'YData',nan)    
    set(ha.ev.min.value,'YData',nan)    
    set(ha.ev.max.name,'String','')
    set(ha.ev.min.name,'String','')
  end
  axis(ha.hs(5),[0.5 length(ev)+0.5 zmin*1.3-zmax*0.3-1e-16 zmax*1.3-zmin*0.3])
end
% To avoid clobbering, zoom will be a dummy variable
stringevnames = ['zoom.evnames = ''', ...
  reshape(strev',1,prod(size(strev))) ''';' 10];
stringev = ['zoom.ev = [' num2str(ev(:)') '];' 10];
stringcallback = [...
  'try % The try avoids errors when figure used for something else' 10 ...
  'zoom.n = ' num2str(size(strev,1)) ';', 10 ...
  'zoom.m = ' num2str(size(strev,2)) ';', 10 ...
  'zoom.h = get(' num2str(ha.figure) ',''Children'');', 10 ...
  'zoom.xx = get(zoom.h(2),''XLim'');', 10 ...
  'zoom.yy = get(zoom.h(2),''YLim'');' 10, ...
  'zoom.nah = 0; % number of available handles' 10, ...
  'zoom.h2ind = ones(1,28);' 10, ...
  'zoom.h2 = get(zoom.h(2),''Children'');', 10 ...
  'zoom.i = 0;', 10 ...
  'while zoom.i < length(zoom.h2)', 10 ...
  '  zoom.i = zoom.i + 1;', 10 ...
  '  if strcmp(get(zoom.h2(zoom.i),''Type''),''text'')', 10 ...
  '    zoom.col = get(zoom.h2(zoom.i),''Color'');', 10 ...
  '    if all(zoom.col == [0 0 0])', 10 ...
  '       set(zoom.h2(zoom.i),''String'','''')', 10 ...
  '       zoom.nah = zoom.nah + 1;' 10, ...
  '       zoom.h2ind(zoom.nah) = zoom.i;' 10, ...
  '    end', 10 ...
  '  end', 10 ...
  'end', 10 ...
  'if zoom.xx(2)-zoom.xx(1) < 25', 10 ...
  '  zoom.fs = min(18,28-(zoom.xx(2)-zoom.xx(1)));', 10 ...
  '  zoom.i1 = max(1,ceil(zoom.xx(1)));', 10 ...
  '  zoom.i2 = min(length(zoom.ev),floor(zoom.xx(2)));', 10 ...
  '  zoom.i = zoom.i1-1;', 10 ...
  '  while zoom.i < zoom.i2', 10 ...
  '    zoom.i = zoom.i + 1;', 10 ...
  '    if zoom.yy(1) < zoom.ev(zoom.i) & zoom.ev(zoom.i) < zoom.yy(2) & zoom.i+1-zoom.i1 <= zoom.nah', 10 ...
  '      zoom.k = zoom.h2ind(zoom.i+1-zoom.i1);', 10 ...
  '      set(zoom.h2(zoom.k), ''Position'',[zoom.i,zoom.ev(zoom.i)]);', 10 ...
  '      set(zoom.h2(zoom.k), ''String'',zoom.evnames(zoom.i*zoom.m+[1-zoom.m:0]));', 10 ...
  '      set(zoom.h2(zoom.k), ''FontSize'',zoom.fs);', 10 ...
  '    end', 10 ...
  '  end', 10 ...
  'end', 10 ...
  'end', 10 ...
  'clear zoom'];
set(zoom(ha.hs(5)),'ActionPostCallback',...
  [stringevnames stringev stringcallback]);

% UPDATE FLUX ERROR
dumzr = e.psizr_err/(e.psibry-e.psimag);
if ~all(ishandle(ha.fluxerror)) | numel(ha.fluxerror) ~= c.nr
  ha.fluxerror = plot(ha.hs(6),dumzr);
else
  for j = 1:c.nr
    set(ha.fluxerror(j),'YData',dumzr(:,j))
  end
end
axis(ha.hs(6),[0.5 c.nz+0.5 min(dumzr(:))-1e-16 max(dumzr(:))+1e-16])

drawnow

if exist('png','var')
  if png
    i = 1;
    plotfig = ['fig', num2str(i), '.png'];
    while exist(plotfig,'file')
      i = i+1;
      plotfig = ['fig', num2str(i), '.png'];
    end
    ffff = getframe(gcf);
    gggg = frame2im(ffff);
    imwrite(gggg,plotfig,'png');
  end
end


