%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_rtplot
%
%  PURPOSE: Real-time plot of plasma during simulation
%
%  INPUTS:  plot_settings (optional)
%           Rlim, Zlim, vvdata, fcdata, ecdata (optional)
%           time-dependent data produced by gs_rtplot_prepare:
%             rbbbs, zbbbs, psizr, cpasma, li, betap
%
%  OUTPUTS: Fast plots in figure window which is current on first call
%           Other figure windows can be used while gs_rtplot is running
%
%  METHOD:  Handles to all objects remembered to allow speedy updates 
	
%
%  WRITTEN BY: Anders Welander ON 6/28/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Configuration of gs_rtplot

if ~exist('plot_settings','var') | ...
   isempty(plot_settings) | ...
   plot_settings.refresh
  
  plot_settings.refresh = false;

  % Boundary
  if ~isfield(plot_settings,'boundary')
    plot_settings.boundary = [];
  end
  if ~isfield(plot_settings.boundary,'LineWidth')
    plot_settings.boundary.LineWidth = 3;
  end
  if ~isfield(plot_settings.boundary,'Color')
    plot_settings.boundary.Color = [0.56 0 1];
  end
  if ~isfield(plot_settings.boundary,'patch')
    plot_settings.boundary.patch = false;
  end
  if ~isfield(plot_settings.boundary,'EdgeColor')
    plot_settings.boundary.EdgeColor = [.8 .6 .1];
  end
  if ~isfield(plot_settings.boundary,'FaceColor')
    plot_settings.boundary.FaceColor = [1 1 .8];
  end

  % Boundary-defining point
  if ~isfield(plot_settings,'bdef')
    plot_settings.bdef = [];
  end
  if ~isfield(plot_settings.bdef,'LineWidth')
    plot_settings.bdef.LineWidth = 3;
  end
  if ~isfield(plot_settings.bdef,'Color')
    plot_settings.bdef.Color = [1 .4 0];
  end
  if ~isfield(plot_settings.bdef,'Marker')
    plot_settings.bdef.Marker = 'x';
  end
  if ~isfield(plot_settings.bdef,'MarkerSize')
    plot_settings.bdef.MarkerSize = 18;
  end

  % Touch point
  if ~isfield(plot_settings,'touch')
    plot_settings.touch = [];
  end
  if ~isfield(plot_settings.touch,'LineWidth')
    plot_settings.touch.LineWidth = 3;
  end
  if ~isfield(plot_settings.touch,'Color')
    plot_settings.touch.Color = [1 .8 0];
  end
  if ~isfield(plot_settings.touch,'Marker')
    plot_settings.touch.Marker = '*';
  end
  if ~isfield(plot_settings.touch,'MarkerSize')
    plot_settings.touch.MarkerSize = 24;
  end

  % Diagnostic points (dp)
  if ~isfield(plot_settings,'dp')
    plot_settings.dp = [];
  end
  if ~isfield(plot_settings.dp,'LineWidth')
    plot_settings.dp.LineWidth = 0;
  end
  if ~isfield(plot_settings.dp,'Color')
    plot_settings.dp.Color = [0 1 0];
  end
  if ~isfield(plot_settings.dp,'Marker')
    plot_settings.dp.Marker = 'x';
  end
  if ~isfield(plot_settings.dp,'MarkerSize')
    plot_settings.dp.MarkerSize = 10;
  end

  % Flux contours
  if ~isfield(plot_settings,'flux_contours')
    plot_settings.flux_contours = [];
  end
  if ~isfield(plot_settings.flux_contours,'psibar')
    plot_settings.flux_contours.psibar = [];
  end
  if ~isfield(plot_settings.flux_contours,'psivalues')
    plot_settings.flux_contours.psivalues = [];
  end
  if ~isfield(plot_settings.flux_contours,'qvalues')
    plot_settings.flux_contours.qvalues = [];
  end
  if ~isfield(plot_settings.flux_contours,'nvac')
    plot_settings.flux_contours.nvac = 33;
  end
  if ~isfield(plot_settings.flux_contours,'npla')
    plot_settings.flux_contours.npla = 33;
  end
  if ~isfield(plot_settings.flux_contours,'LineWidth')
    plot_settings.flux_contours.LineWidth = 0;
  end
  if ~isfield(plot_settings.flux_contours,'fill')
    plot_settings.flux_contours.fill = false;
  end
  if plot_settings.flux_contours.LineWidth == 0
    plot_settings.flux_contours.linemax = 0;
  else
    plot_settings.flux_contours.linemax = ...
      length(plot_settings.flux_contours.psibar) + ...
      length(plot_settings.flux_contours.psivalues) + ...
      length(plot_settings.flux_contours.qvalues) + ...
      max(plot_settings.flux_contours.nvac, ...
          plot_settings.flux_contours.npla);
  end
  
  % Camera
  if ~isfield(plot_settings,'camera')
    plot_settings.camera = [];
  end
  if ~isfield(plot_settings.camera,'showit')
    plot_settings.camera.showit = 0;
  end
  if ~isfield(plot_settings.camera,'animate')
    plot_settings.camera.animate = 1;
  end
  if ~isfield(plot_settings.camera,'nxpix')
    plot_settings.camera.nxpix = 64;
  end
  if ~isfield(plot_settings.camera,'nzpix')
    plot_settings.camera.nzpix = 128;
  end

  % Limiter
  if ~isfield(plot_settings,'limiter')
    plot_settings.limiter = [];
  end
  if ~isfield(plot_settings.limiter,'LineWidth')
    if plot_settings.camera.showit
      plot_settings.limiter.LineWidth = 6;
    else
      plot_settings.limiter.LineWidth = 3;
    end
  end
  if ~isfield(plot_settings.limiter,'Color')
    if plot_settings.camera.showit
      plot_settings.limiter.Color = [0 1 0];
    else
      plot_settings.limiter.Color = [0.5, 0.5, 0.5];
    end
  end
  if ~isfield(plot_settings.limiter,'patch')
    plot_settings.limiter.patch = false;
  end
  if ~isfield(plot_settings.limiter,'EdgeColor')
    plot_settings.limiter.EdgeColor = [0.5 0.5 0.5];
  end
  if ~isfield(plot_settings.limiter,'FaceColor')
    plot_settings.limiter.FaceColor = [0.7 0.7 0.7];
  end
  if ~isfield(plot_settings.limiter,'GlowColor')
    plot_settings.limiter.GlowColor = [0.9 0.8 1.0];
  end
  if length(plot_settings.limiter.GlowColor) ~= 3 || ...
    any(plot_settings.limiter.GlowColor < 0 | ...
    plot_settings.limiter.GlowColor > 1)
    disp(['Ouch in gs_rtplot: ', ...
      'plot_settings.limiter.GlowColor must be 3 values between 0 and 1'])
    disp('Using default values...')
    plot_settings.limiter.GlowColor = [0.9 0.8 1.0];
  end

  % Vessel
  if ~isfield(plot_settings,'vessel')
    plot_settings.vessel = [];
  end
  if ~isfield(plot_settings.vessel,'LineWidth')
    plot_settings.vessel.LineWidth = 0;
  end
  if ~isfield(plot_settings.vessel,'Color')
    plot_settings.vessel.Color = [.1 .1 .8];
  end
  if ~isfield(plot_settings.vessel,'patch')
    plot_settings.vessel.patch = false;
  end
  if ~isfield(plot_settings.vessel,'EdgeColor')
    plot_settings.vessel.EdgeColor = [0 1 1];
  end
  if ~isfield(plot_settings.vessel,'FaceColor')
    plot_settings.vessel.FaceColor = [0 0 1];
  end

  % Coils
  if ~isfield(plot_settings,'coils')
    plot_settings.coils = [];
  end
  if ~isfield(plot_settings.coils,'LineWidth')
    plot_settings.coils.LineWidth = 0;
  end
  if ~isfield(plot_settings.coils,'Color')
    plot_settings.coils.Color = [1 .1 0];
  end
  if ~isfield(plot_settings.coils,'patch')
    plot_settings.coils.patch = false;
  end
  if ~isfield(plot_settings.coils,'EdgeColor')
    plot_settings.coils.EdgeColor = [0 1 1];
  end
  if ~isfield(plot_settings.coils,'FaceColor')
    plot_settings.coils.FaceColor = [0 0 1];
  end

  % Text
  if ~isfield(plot_settings,'text')
    plot_settings.text = [];
  end
  if ~isfield(plot_settings.text,'FontSize')
    plot_settings.text.FontSize = 15;
  end

  % Profiles
  if ~isfield(plot_settings,'profiles')
    plot_settings.profiles = 0;
  end

  % Movie making
  if ~isfield(plot_settings,'save_frames')
    plot_settings.save_frames = 0;
  end
  if isfield(plot_settings,'axis')
    if length(plot_settings.axis) ~= 4 || ...
      plot_settings.axis(1) >= plot_settings.axis(2) || ...
      plot_settings.axis(3) >= plot_settings.axis(4)
      plot_settings = rmfield(plot_settings,'axis');
      disp('Warning gs_rtplot: Invalid plot_settings.axis removed')
      disp('Format should be [xmin xmax ymin ymax] with xmin<xmax & ymin<ymax')
    end
  end
  if ~isfield(plot_settings,'frame_name_prefix')
    plot_settings.frame_name_prefix = 'gsev';
  end
end


% Handles to all displayed objects

if ~exist('rtplot_handles','var') | ...
   isempty(rtplot_handles) | ...
   ~isfield(rtplot_handles,'figure') | ...
   ~ishandle(rtplot_handles.figure) | ...
   ~isfield(rtplot_handles,'refresh') | ...
   rtplot_handles.refresh
  rtplot_handles.figure = gcf;
  if plot_settings.profiles & exist('index_in_y','var') & ...
    isfield(index_in_y,'T') & isfield(index_in_y,'jpar') & ...
    isfield(index_in_y,'vind') & isfield(index_in_y,'vres')
    rtplot_handles.profiles = true;
  else
    rtplot_handles.profiles = false;
  end
  clf
  if rtplot_handles.profiles
    subplot(1,3,1)
    hold on
    rtplot_handles.vind = plot(nan(1,nr),'b');
    rtplot_handles.vres = plot(nan(1,nr),'r');
    title('{\color{blue}-V_i_n_d}, {\color{red}V_r_e_s}');
    xlabel('\rho')
    subplot(1,3,2)
    hold on
    rtplot_handles.jpar = plot(nan(1,nr),'color',[0,0.6,0]);
    title('{\color[rgb]{0,0.6,0}j_p_a_r}');
    xlabel('\rho')
    subplot(1,3,3)
  end
  hold on
  rtplot_handles.axes = get(rtplot_handles.figure,'Children');
  if ishandle(rtplot_handles.axes(1))
    if rtplot_handles.profiles
      set(rtplot_handles.axes(1),'Position',[0.65 0.11 0.25 0.815])
    else
      set(rtplot_handles.axes(1),'Position',[0.05 0.11 0.75 0.815])
    end
  end
  rtplot_handles.boundary = nan;
  rtplot_handles.touch = nan;
  rtplot_handles.dp = nan;
  rtplot_handles.bdef = nan;
  rtplot_handles.flux_contours = nan(1,plot_settings.flux_contours.linemax);
  rtplot_handles.camera = nan;
  rtplot_handles.limiter = nan;
  rtplot_handles.vessel = nan;
  rtplot_handles.coils = nan;
  rtplot_handles.text = nan(1,2);
  rtplot_handles.title = nan;
  rtplot_handles.refresh = false;
  set(rtplot_handles.figure,'Name','gsevolve simulation')
end

% PROFILES
if rtplot_handles.profiles
  T = y(index_in_y.T);
  rhot = sqrt(T/T(end));
  if ishandle(rtplot_handles.vind)
    set(rtplot_handles.vind,'XData',rhot)
    set(rtplot_handles.vind,'YData',-y(index_in_y.vind))
  end
  if ishandle(rtplot_handles.vres)
    set(rtplot_handles.vres,'XData',rhot)
    set(rtplot_handles.vres,'YData',y(index_in_y.vres))
  end
  if ishandle(rtplot_handles.jpar)
    set(rtplot_handles.jpar,'XData',rhot)
    set(rtplot_handles.jpar,'YData',y(index_in_y.jpar))
  end
end

% CAMERA
if plot_settings.camera.showit
  nxpix = plot_settings.camera.nxpix;
  nzpix = plot_settings.camera.nzpix;
  imrgb = zeros(nzpix,nxpix,3);
  if ~ishandle(rtplot_handles.camera)
    figure(rtplot_handles.figure)
    hold on
    gs_limiter_radiation % Calculate Plg
    gs_render_limiter
    rtplot_handles.camera = image(xpix,zpix,imrgb);
  end
  if plasma
    im = Pcg(:,iplasma)*(1-psibarzr(iplasma));
    dum = max(im);
    if dum > 0 && min(im)>= 0
      im = im/dum;
    else % psibarzr(iplasma) messed up
      im = zeros(nzpix,nxpix);
    end
    imr = reshape(im,nzpix,nxpix)*plot_settings.limiter.GlowColor(1);
    img = reshape(im,nzpix,nxpix)*plot_settings.limiter.GlowColor(2);
    imb = reshape(im,nzpix,nxpix)*plot_settings.limiter.GlowColor(3);
  else
    im = ones(nzpix,nxpix)*0.2;
    imr = im;
    img = im;
    imb = im;
  end
  imrgb(:,:,1) = imr.*fpix + 1-fpix;
  imrgb(:,:,2) = img.*fpix + 1-fpix;
  imrgb(:,:,3) = imb.*fpix + 1-fpix;
  set(rtplot_handles.camera,'CData',imrgb)
end

% FLUX CONTOURS
if plot_settings.flux_contours.LineWidth
  if length(plot_settings.flux_contours.psivalues) > 1
    psia = max(plot_settings.flux_contours.psivalues);
    psib = min(plot_settings.flux_contours.psivalues);
  else
    psia = max(psizr(:));
    psib = min(psizr(:));
  end
  map = colormap;
  maplevels = linspace(psib,psia,size(map,1));
  rcont2plot =   nan(plot_settings.flux_contours.linemax, nbbbs_max);
  zcont2plot =   nan(plot_settings.flux_contours.linemax, nbbbs_max);
  color2plot = zeros(plot_settings.flux_contours.linemax, 3);
  if length(plot_settings.flux_contours.psivalues) > 1
    dumcontour = contourc(rg, zg, psizr, ...
      plot_settings.flux_contours.psivalues);
  else
    if plasma
      n = plot_settings.flux_contours.npla;
    else
      n = plot_settings.flux_contours.nvac;
    end
    if n > 0 && psia > psib
      dumcontour = contourc(rg, zg, psizr, n);
    else
      dumcontour = zeros(2,0);
    end
  end
  i = 0; % index to last filled row of rcont2plot, zcont2plot, color2plot
  j = 1;
  while i < plot_settings.flux_contours.linemax & j < size(dumcontour,2)
    i = i+1;
    dum = dumcontour(1,j);
    color2plot(i,:) = interp1(maplevels,map,dum);
    n = 0;
    while j < size(dumcontour,2) && dum == dumcontour(1,j) && n < nbbbs_max
      k = min(n+dumcontour(2,j),nbbbs_max)-n;
      rcont2plot(i,n+(1:k)) = dumcontour(1,j+(1:k));
      zcont2plot(i,n+(1:k)) = dumcontour(2,j+(1:k));
      j = j + 1 + dumcontour(2,j);
      n = n + 1 + k;
    end
  end
  if any(~ishandle(rtplot_handles.flux_contours))
    figure(rtplot_handles.figure)
    hold on
    if ishandle(rtplot_handles.boundary) % Boundary underneath contours
      delete(rtplot_handles.boundary) % Will be regenerated on top of contours
    end
  end
  for i = 1:length(rtplot_handles.flux_contours)
    if ishandle(rtplot_handles.flux_contours(i))
      set(rtplot_handles.flux_contours(i), ...
	'XData', rcont2plot(i,:), ...
	'YData', zcont2plot(i,:), ...
	'Color', color2plot(i,:) ...
	)
    else
      rtplot_handles.flux_contours(i) = plot(...
	rcont2plot(i,:), zcont2plot(i,:), ...
	'LineWidth', plot_settings.flux_contours.LineWidth, ...
	'Color', color2plot(i,:));      
    end
  end
end

% LIMITER
if ~ishandle(rtplot_handles.limiter) & plot_settings.limiter.LineWidth
  figure(rtplot_handles.figure) 
  hold on 
  if plot_settings.limiter.patch
    rtplot_handles.limiter = patch(rl, zl, ...
      plot_settings.limiter.FaceColor, ...
      'LineWidth', plot_settings.limiter.LineWidth, ...
      'EdgeColor', plot_settings.limiter.EdgeColor);  
  else
    rtplot_handles.limiter = plot(rl, zl, ...
      'LineWidth', plot_settings.limiter.LineWidth, ...
      'Color', plot_settings.limiter.Color);
  end
end

% DIAGNOSTIC POINTS
if plot_settings.dp.LineWidth
  if ~ishandle(rtplot_handles.dp)
    figure(rtplot_handles.figure)
    hold on
    rtplot_handles.dp = plot(rdp, zdp, plot_settings.dp.Marker, ...
      'MarkerSize', plot_settings.dp.MarkerSize, ...
      'LineWidth', plot_settings.dp.LineWidth, ...
      'Color', plot_settings.dp.Color);
    axis image
    zoom on
  end
end

% BOUNDARY
if plot_settings.boundary.LineWidth
  if ~ishandle(rtplot_handles.boundary)
    figure(rtplot_handles.figure)
    hold on
    if plot_settings.boundary.patch
      rtplot_handles.boundary = patch(nan, nan, ...
      	plot_settings.boundary.FaceColor, ...
	'LineWidth', plot_settings.boundary.LineWidth, ...
	'EdgeColor', plot_settings.boundary.EdgeColor);
    else
      rtplot_handles.boundary = plot(nan, nan, ...
	'LineWidth', plot_settings.boundary.LineWidth, ...
	'Color', plot_settings.boundary.Color);
    end
    axis image
    zoom on
  end
  set(rtplot_handles.boundary, ...
    'XData', rbbbs(1:nbbbs), ...
    'YData', zbbbs(1:nbbbs));
end

% TOUCH POINT
if plot_settings.touch.LineWidth
  if ~ishandle(rtplot_handles.touch)
    figure(rtplot_handles.figure)
    hold on
    rtplot_handles.touch = plot(rbdef, zbdef);
    axis image
    zoom on
  end
  if exist('shape','var') & strcmp(shape,'LIM')
    set(rtplot_handles.touch, ...
      'XData', rbdef, ...
      'YData', zbdef, ...
      'Marker', plot_settings.touch.Marker, ...
      'MarkerSize', plot_settings.touch.MarkerSize, ...
      'LineWidth', plot_settings.touch.LineWidth, ...
      'Color', plot_settings.touch.Color);
  else % Here bdef is not a touch point, so plot it the same as bdef
    set(rtplot_handles.touch, ...
      'XData', rbdef, ...
      'YData', zbdef, ...
      'Marker', plot_settings.bdef.Marker, ...
      'MarkerSize', plot_settings.bdef.MarkerSize, ...
      'LineWidth', plot_settings.bdef.LineWidth, ...
      'Color', plot_settings.bdef.Color);
  end
end

% BOUNDARY DEFINING POINT
if plot_settings.bdef.LineWidth
  if ~ishandle(rtplot_handles.bdef)
    figure(rtplot_handles.figure)
    hold on
    rtplot_handles.bdef = plot(rbdef, zbdef, ...
      'Marker', plot_settings.bdef.Marker, ...
      'MarkerSize', plot_settings.bdef.MarkerSize, ...
      'LineWidth', plot_settings.bdef.LineWidth, ...
      'Color', plot_settings.bdef.Color);
    axis image
    zoom on
  else
    set(rtplot_handles.bdef, ...
      'XData', rbdef, ...
      'YData', zbdef)
  end
end

% VESSEL
if plot_settings.vessel.LineWidth
  if ~ishandle(rtplot_handles.vessel)
    if exist('vvdata','var') & ~isempty(vvdata)
      figure(rtplot_handles.figure)
      hold on
      if plot_settings.vessel.patch
      else
	rtplot_handles.vessel = ...
          plot_boxx(vvdata(2,:),vvdata(1,:),vvdata(4,:),vvdata(3,:),...
	  'g',vvdata(5,:),vvdata(6,:));

	set(rtplot_handles.vessel, ...
	  'LineWidth', plot_settings.vessel.LineWidth, ...
	  'Color', plot_settings.vessel.Color);
      end
      axis image
      zoom on
    end
  end
end

% COILS
if plot_settings.coils.LineWidth
  if ~ishandle(rtplot_handles.coils)
    if exist('ecdata','var') & ~isempty(ecdata)
      figure(rtplot_handles.figure)
      hold on
      if plot_settings.coils.patch
      else
	rtplot_handles.coils = ...
          plot_boxx(ecdata(2,:),ecdata(1,:),ecdata(4,:),ecdata(3,:));

	set(rtplot_handles.coils, ...
	  'LineWidth', plot_settings.coils.LineWidth, ...
	  'Color', plot_settings.coils.Color);
      end
    end
    if exist('fcdata','var') & ~isempty(fcdata)
      figure(rtplot_handles.figure)
      hold on
      if plot_settings.coils.patch
      else
	rtplot_handles.coils = ...
          plot_boxx(fcdata(2,:),fcdata(1,:),fcdata(4,:),fcdata(3,:),...
	  'g',fcdata(5,:),fcdata(6,:));

	set(rtplot_handles.coils, ...
	  'LineWidth', plot_settings.coils.LineWidth, ...
	  'Color', plot_settings.coils.Color);
      end
    end
    axis image
    zoom on
  end
end

if ~ishandle(rtplot_handles.title)
  rtplot_handles.title = title('','FontSize',18);
end
if exist('time','var') & time > -inf
  if ~exist('titstrlen','var') || isempty(titstrlen)
    titstrlen = 0;
  end
  titstr = ['time = ' num2str(time)];
  if length(titstr) >= titstrlen
    titstrlen = length(titstr);
  else
    titstr(end+1:titstrlen) = '0';
  end
  % Present the time counter in a nice way
  if dtplot > 0
    titstr = sprintf(title_time_format,time);
  end
else
  titstr = '';
end
set(rtplot_handles.title,'String',titstr)

% First text handle gives values of Ip, li, betap or the text vacuum
if ~ishandle(rtplot_handles.text(1))
  rtplot_handles.text(1) = text(0,0,'VACUUM',...
    'HorizontalAlignment','Left',...
    'VerticalAlignment','Top',...
    'BackgroundColor',[1 1 .9],...
    'EdgeColor',[0 0.8 0],...
    'LineWidth',2,...
    'FontWeight','Bold',...
    'FontSize',plot_settings.text.FontSize);
end
if plasma
  strIp = sprintf(' %2.0f A ',cpasma);
  strli = sprintf(' l_i = %2.2f ',li);
  strbp = sprintf(' \\beta_p = %2.2f ',betap);
  strIlbv = [strIp '\newline' strli '\newline' strbp];
else
  strIlbv = 'VACUUM';
end
set(rtplot_handles.text(1), 'String',strIlbv)

% Second text shows gamma and if plasma is artificially stabilized
if ~ishandle(rtplot_handles.text(2))
  rtplot_handles.text(2) = text(0,0,'Announcement',...
    'HorizontalAlignment','Left',...
    'LineWidth',2,...
    'FontWeight','Bold',...
    'FontSize',plot_settings.text.FontSize);
  % Store handles so that the zoom feature below always works
  set(rtplot_handles.figure,'UserData',rtplot_handles.text)
  % A zoom callback function
  if isa(rtplot_handles.figure,'matlab.ui.Figure')
    figno = num2str(rtplot_handles.figure.Number);
  else
    figno = num2str(rtplot_handles.figure);
  end
  set(zoom(gca),'ActionPostCallback',...
    [...
    'if gcf == ', figno, 10, ...
    '  zoom.text = get(' figno ',''UserData'');', 10, ...
    '  if length(zoom.text) == 2 & all(ishandle(zoom.text))', 10, ...
    '    zoom.axis = axis;', 10, ...
    '    zoom.x = 1.02*zoom.axis(2)-0.02*zoom.axis(1);', 10, ...
    '    zoom.y1 = zoom.axis(4);', 10, ...
    '    set(zoom.text(1), ''Position'',[zoom.x zoom.y1])', 10, ...
    '    zoom.Extent = get(zoom.text(1),''Extent'');', 10, ...
    '    zoom.y2 = zoom.y1-zoom.Extent(4)*1.1;', 10, ...
    '    set(zoom.text(2), ''Position'',[zoom.x zoom.y2])', 10, ...
    '  end', 10, ...  
    '  clear zoom', 10, ...  
    'end', 10, ...
    ]);
end
if stabilized
  set(rtplot_handles.text(2),...
    'String','Artificial\newlinestability',...
    'HorizontalAlignment','Left',...
    'VerticalAlignment','Top',...
    'BackgroundColor',[1 1 0],...
    'EdgeColor',[1 0 0]);
else
  set(rtplot_handles.text(2),...
    'String',['\gamma_z = ' num2str(real(gamma))],...
    'HorizontalAlignment','Left',...
    'VerticalAlignment','Top',...
    'BackgroundColor',[1 1 .9],...
    'Color',[0 0 0],...
    'EdgeColor',[0 0.8 0.8]);
end
% Position text
if gcf == rtplot_handles.figure
  dum4 = axis;
  xtext = 1.02*dum4(2)-0.02*dum4(1);
  ytext1 = dum4(4);
  set(rtplot_handles.text(1), 'Position',[xtext ytext1])
  dum4b = get(rtplot_handles.text(1),'Extent');
  ytext2 = ytext1-dum4b(4)*1.1;
  set(rtplot_handles.text(2), 'Position',[xtext ytext2])
end

if plot_settings.save_frames
  set(rtplot_handles.figure,'PaperPositionMode','auto')
  if isfield(plot_settings,'axis')
    axis(plot_settings.axis)
  end
  drawnow
  print(rtplot_handles.figure,'-dpng',...
    [plot_settings.frame_name_prefix num2str(plot_counter-1) '.png'])
else
  drawnow
end
