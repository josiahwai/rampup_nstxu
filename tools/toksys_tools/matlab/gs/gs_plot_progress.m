%  USAGE:   gs_plot_progress
%
%  PURPOSE: Plot progress of Grad-Shafranov (gs) calculations
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
%  OUTPUTS: plots of pprime, ffprim, scalar values,
%             contours, conductor currents, flux error
%	
%  METHOD:  A resize function (ResizeFcn) scales fonts with window size 
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	4/12/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('ResizeFcn') ~= 1
  ResizeFcn = [...
    'data4ResizeFcn.fpos = get(gcf,''position'');', 10, ...
    'data4ResizeFcn.fs1 = round(data4ResizeFcn.fpos(3)/60);', 10, ...
    'data4ResizeFcn.fs2 = round(data4ResizeFcn.fpos(3)/100);', 10, ...
    'data4ResizeFcn.h = get(gcf,''UserData'');', 10, ...
    'if ~isempty(data4ResizeFcn.h)', 10, ...
    '  data4ResizeFcn.j = 1;', 10, ...
    '  while data4ResizeFcn.j <= data4ResizeFcn.h(1)', 10, ...
    '    data4ResizeFcn.j = data4ResizeFcn.j+1;', 10, ...
    '    if ishandle(data4ResizeFcn.h(data4ResizeFcn.j))', 10, ...
    '      set(data4ResizeFcn.h(data4ResizeFcn.j),''fontsize'',data4ResizeFcn.fs1);', 10, ...
    '    end', 10, ...
    '  end', 10, ...
    '  while data4ResizeFcn.j < length(data4ResizeFcn.h)', 10, ...
    '    data4ResizeFcn.j = data4ResizeFcn.j+1;', 10, ...
    '    if ishandle(data4ResizeFcn.h(data4ResizeFcn.j))', 10, ...
    '      set(data4ResizeFcn.h(data4ResizeFcn.j),''fontsize'',data4ResizeFcn.fs2);', 10, ...
    '    end', 10, ...
    '  end', 10, ...
    'end', 10, ...
    'clear data4ResizeFcn'];
end
hs = NaN*ones(1,10); % handles to subplots

h1 = []; % handles to objects with font size fs1
h2 = []; % handles to objects with font size fs2

if exist('htt','var')
  if ishandle(htt)
    if ~isempty(get(htt,'UserData'))
      not_done = 0;
    end
  end
end

if isempty(get(0,'CurrentFigure'))
  figure_handle = figure;
  fpos = get(figure_handle,'position');
  fpos(4) = fpos(3)/2;
  fpos(3:4) = [1000 500];
  set(figure_handle,'position',fpos)
end
if ~exist('figure_handle','var')
  figure_handle = gcf;
end
if ~ishandle(figure_handle)
  figure(figure_handle)
end

fpos = get(figure_handle,'position');
if 0
if ~exist('ht','var') | ~ishandle(ht)
  ht = uitoolbar(figure_handle);
  htt = uitoggletool(ht);
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
  set(htt,'Cdata',cdata);
  set(htt,'ClickedCallback','set(gcbo,''Userdata'',1);')
  set(htt,'Enable','on');
  set(htt,'TooltipString','Stop iterating');
  set(figure_handle,'position',fpos);
end
end

if exist('figure_message','var')
  set(figure_handle,'Name',figure_message)
elseif exist('iteration_counter','var')
  set(figure_handle,'Name',['Iteration ' num2str(iteration_counter)])
end

ix = 0; % index of x-point that defines boundary
for i = 1:nulls.count
  if nulls.r(i) == rbdef & nulls.z(i) == zbdef
    ix = i;
  end
end
strix = ['_x_' num2str(ix) ' '];

if ishandle(hs(1))
  delete(get(hs(1),'Children'))
end
hs(1) = subplot(2,4,1);
set(hs(1),'position',[0.05 0.6 0.2 0.3])
plot(sqrt(psibar),pprime,'linew',3)
if max(pprime) > min(pprime)
  axis([0 1 min(pprime) max(pprime)])
end
grid(hs(1),'on')
h1(end+1) = title('p''(\rho)');

if ishandle(hs(2))
  delete(get(hs(2),'Children'))
end
hs(2) = subplot(2,4,2);
set(hs(2),'position',[0.3 0.6 0.2 0.3])
plot(sqrt(psibar),ffprim,'linew',3)
if max(ffprim) > min(ffprim)
  axis([0 1 min(ffprim) max(ffprim)])
end
h1(end+1) = title('ff''(\rho)');
grid(hs(2),'on')

if ishandle(hs(3))
  delete(hs(3))
end
hs(3) = subplot(2,4,3);
set(hs(3),'position',[0.55 0.1 0.175 0.8])
plot([0 1],[0 1],'w')
axis(hs(3),[0 1 0 1])
set(hs(3),'Xtick',[])
set(hs(3),'Ytick',[])
report.name = strvcat('I_p',          '\beta_p','l_i','\psi_a','\psi_b');
report.expr = strvcat('round(cpasma)','betap',  'li', 'psimag', 'psibry');
report.form = strvcat('%8d','%7.5f',  '%7.5f','%7.5f','%7.5f','%7.5f');
for i = 1:nulls.count
  report.name = strvcat(report.name,['R_x_' num2str(i)],['Z_x_' num2str(i)]);
  report.expr = strvcat(report.expr,['nulls.r(' num2str(i) ')'],['nulls.z(' num2str(i) ')']);
  report.form = strvcat(report.form,'%7.5f','%7.5f');
end
nreport = size(report.name,1);
for j = 1:nreport
  h1(end+1) = text(0.1,0.985-j/nreport*0.9,report.name(j,:));
  if strcmp(report.name(j,2:4),'_x_')
    set(h1(end),'color','blue')
  end
  if strcmp(report.name(j,2:1+length(strix)),strix)
    set(h1(end),'color','red')
  end
  val = eval(report.expr(j,:));
  h1(end+1) = text(0.4,0.990-j/nreport*0.9,sprintf(report.form(j,:),val));
  if strcmp(report.name(j,2:4),'_x_')
    set(h1(end),'color','blue')
  end
  if strcmp(report.name(j,2:1+length(strix)),strix)
    set(h1(end),'color','red')
  end
end
if exist('time','var')
  h1(end+1) = title(['time = ' num2str(time)]);
end

if ishandle(hs(4))
  delete(hs(4))
end
hs(4) = subplot(2,4,4);
set(hs(4),'position',[0.765 0.1 0.225 0.8])
hold on
%plot(rg(1)+ones(2,1)*linspace(-0.5,nr-0.5,nr+1)*dr,[zg(1)-dz/2; zg(nz)+dz/2]*ones(1,nr+1),'--','color',[.7 .7 .7])
%plot([rg(1)-dr/2; rg(nr)+dr/2]*ones(1,nz+1),zg(1)+ones(2,1)*linspace(-0.5,nz-0.5,nz+1)*dz,'--','color',[.7 .7 .7])
%plot(rgg,zgg,'rx')
plot(rl,zl,'k','linew',6)
patch(Rlim,Zlim,[0.8 0.8 0.8])
patch(rbbbs(1:nbbbs),zbbbs(1:nbbbs),[1 1 0.9])
if ~plasma
  psia = max(psizr(ilimgg>-1));
  psib = min(psizr(ilimgg>-1));
else
  psia = psimag;
  psib = psibry;
end
if exist('nhires','var') & nhires > 1
  nzhr = 1+nhires*(nz-1);
  nrhr = 1+nhires*(nr-1);
  rghr = linspace(rg(1),rg(nr),nrhr);
  zghr = linspace(zg(1),zg(nz),nzhr)';
  px = gs_interp2(rg,zg,psizr,ones(nzhr,1)*rghr, zghr*ones(1,nrhr));
  contour(rghr,zghr,1-(px-psia)/(psib-psia),[0:7]/8,'linew',3)
else
  contour(rg,zg,1-(psizr-psia)/(psib-psia),[0:7]/8,'linew',3)
end
for i = nulls.count:-1:1
  if i ~= ix
    plot(nulls.r(i),nulls.z(i),'bx','markersize',12,'linew',2)
  end
end
if exist('targets','var')
  if isfield(targets,'rsnf') & isfield(targets,'zsnf')
    plot(targets.rsnf,targets.zsnf,'w*','markersize',24,'linew',4)
  end
  if isfield(targets,'rx') & isfield(targets,'zx')
    plot(targets.rx,targets.zx,'yx','markersize',18,'linew',3)
  end
  if isfield(targets,'rsep') & isfield(targets,'zsep')
    plot(targets.rsep,targets.zsep,'mo','markersize',9,'linew',2)
  end
end
plot(rbbbs(1:nbbbs),zbbbs(1:nbbbs),'linew',3,'color',[0,0, 0.5625])
plot(rbdef,zbdef,'yx','markers',24,'linew',5)
plot(rbdef,zbdef,'rx','markers',21,'linew',3)
axis('image')
axis([rg(1) rg(nr) zg(1) zg(nz)])
h1(end+1) = title('Flux contours');


if ishandle(hs(5))
  delete(hs(5))
end
hs(5) = subplot(2,4,5);
set(hs(5),'position',[0.05 0.1 0.2 0.3])
hold on
if exist('limits','var') & isfield(limits,'ic')
  for i = 1:nc
    if ic(i) < limits.ic(i,1)+1e-9 | ic(i) > limits.ic(i,2)-1e-9
      plot(i,ic(i),'yo--','MarkerSize',18,'MarkerFaceColor','yellow','linew',3)
    end
  end
  plot(limits.ic(:,1),'b--')
  plot(limits.ic(:,2),'r--')
  plot(limits.ic(:,1),'bx','MarkerSize',12,'linew',3)
  plot(limits.ic(:,2),'rx','MarkerSize',12,'linew',3)
end
plot([ic; iv],'g--')
plot([ic; iv],'go','MarkerSize',4,'MarkerFaceColor','green','linew',3)
h1(end+1) = title('Conductors');
if all(iv == 0)
  xmax = nc+1/2;
  ymin = 1.2*min(ic)-0.2*max(ic);
  ymax = 1.2*max(ic)-0.2*min(ic);
else
  xmax = nc+nv+1/2;
  ymin = 1.2*min([ic; iv])-0.2*max([ic; iv]);
  ymax = 1.2*max([ic; iv])-0.2*min([ic; iv]);
end
axis(hs(5),[0.5 xmax ymin ymax])
grid(hs(5),'on')


if ishandle(hs(6))
  delete(hs(6))
end
hs(6) = subplot(2,4,6);
set(hs(6),'position',[0.3 0.1 0.2 0.3])
if circle_model
  plot(psiherr/(psibry-psimag))
else
  plot(psizr_err/(psibry-psimag))
end
h1(end+1) = title('\psi error');
ax = axis;
axis([0.5 nz+0.5 ax(3:4)])
grid(hs(6),'on')

h1 = [h1 hs([1 2 3 5 6])];
h2 = [h2 hs(4)];

set(figure_handle,'UserData',[length(h1) h1 h2])
eval(ResizeFcn)
set(figure_handle,'ResizeFcn',ResizeFcn)
zoom on

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
