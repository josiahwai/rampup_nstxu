%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_configure
%
%  PURPOSE: Configure parameters that are used by the gs_*.m codes
%           Reconfiguration possible by invoking again with different config
%
%  INPUTS: config, structure with TokSys tokamak info (tok_data_struct)
%
%  OUTPUTS: Parameters for analysis of equilibrium and its response
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%  WRITTEN BY:  Anders Welander  ON	10/19/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If gs_configure is invoked for the first time by a function, all variables
% may exist but be empty, if they were declared persistent by the function.
% To ensure that the code runs the same from such a function or otherwise,
% all actions are independent of whether a variable is empty or non-existent.

% Variables that are specified in config will be updated.
% Variables that aren't specified in config and don't exist (or are empty) 
% will be assigned defaults.
% Variables that aren't specified in config and do exist will be checked for
% consistency with the configuration and assigned new defaults if needed.
% All persistent variables have their final fixed sizes after configuration

% Verbose if set will report about every recognized field in config
% and warn about unrecognized fields
if isfield(config,'verbose')
  verbose = config.verbose;
else
  verbose = 0;
end

mu0 = 0.4e-6*pi;
twopi = 2*pi;

% Update psizr if grid is being reconfigured
if isfield(config,'rg') & isfield(config,'zg')
  if exist('rg','var') & exist('zg','var') & exist('psizr','var') & ...
    length(rg) == size(psizr,2) & length(zg) == size(psizr,1) & ...
    ~isempty(psizr)
    if verbose > 0
      disp('psizr for old grid being interpolated for new grid')
    end
    psizr = gs_interp2(rg,zg,psizr,...
      ones(length(config.zg),1)*config.rg(:)', ...
      config.zg(:)*ones(1,length(config.rg)));
  end
end

if isfield(config,'tokamak')
  tokamak = config.tokamak;
end

if ~exist('tokamak','var') | isempty(tokamak)
  tokamak = 'Unknown';
  if verbose > 0
    disp('The name of this tokamak is not known, config.tokamak is missing')
  end
end

if isfield(config,'min_iterations')
  min_iterations = config.min_iterations;
end
if ~exist('min_iterations','var') | isempty(min_iterations)
  min_iterations = 1;
end

if isfield(config,'max_iterations')
  max_iterations = config.max_iterations;
end
if ~exist('max_iterations','var') | isempty(max_iterations)
  max_iterations = 9;
end

if isfield(config,'constraints')
  constraints = config.constraints;
end
if ~exist('constraints','var') | isempty(constraints)
  if exist('xc','var')
    if any([isfield(xc,'sp') isfield(xc,'sf')])
      constraints = 0;
      disp('constraints = 0, xc.sp, xc.sf specify pressure and current profiles')
    elseif any([isfield(xc,'ip') isfield(xc,'li') isfield(xc,'betap')])
      constraints = 1;
      disp('constraints = 1, xc.ip, xc.li, xc.betap specify current and pressure profiles')
    elseif any([isfield(xc,'ip') isfield(xc,'li') isfield(xc,'betap')])
      constraints = 2;
      disp('constraints = 2, xc.ip, xc.li, xc.Wth specify pressure and current profiles')
    else
      constraints = 1;
      disp('constraints = 1, xc.ip, xc.li, xc.betap specify current and pressure profiles')
    end
  end
end
if ~exist('constraints','var') | isempty(constraints)
  constraints = 1;
end
if verbose > 0
  if constraints == 1
    disp('constraints = 1, only 3 degrees of freedom for internal profiles')
    disp('described by cpasma, li, betap')
  elseif constraints == 2
    disp('constraints = 2, only 3 degrees of freedom for internal profiles')
    disp('described by cpasma, li, Wth')
  end
end

% Read the grid from config
if isfield(config,'rg')
  rg = config.rg(:)';
end
if ~exist('rg','var') | isempty(rg)
  error('Grid info required in first call. Field rg (and maybe others) are missing from config.')
end
nr = length(rg);
dr = (rg(nr)-rg(1))/(nr-1);
rge = [rg(1)-dr/2 rg+dr/2]; % positions of nr+1 edges of grid rectangles
if isfield(config,'zg')
  zg = config.zg(:);
end
if ~exist('rg','var') | isempty(rg)
  error('Grid info required in first call. Field zg (and maybe others) are missing from config.')
end
nz = length(zg);
dz = (zg(nz)-zg(1))/(nz-1);
zge = [zg(1)-dz/2; zg+dz/2]; % positions of nz+1 edges of grid rectangles
ngg = nr*nz;
rgg = ones(nz,1)*rg;
zgg = zg*ones(1,nr);
Ag = dr*dz;
RA = rgg*Ag; % Surface integral of R over cell
AR = (log(rgg+dr/2)-log(rgg-dr/2))*dz; % Surface integral of 1/R over cell
if ~exist('pprime','var') || isempty(pprime)
  pprime = zeros(nr,1);
end
if length(pprime) ~= nr % Allow reconfiguration that preserve the equilibrium
  pprime = zeros(nr,1);
end
if ~exist('rmaxis','var') | isempty(rmaxis)
  rmaxis = mean(rg);
end
if ~exist('zmaxis','var') | isempty(zmaxis)
  zmaxis = mean(zg);
end
if ~exist('cpasma','var') | isempty(cpasma)
  cpasma = 0;
end
if ~exist('plasma','var') | isempty(plasma)
  plasma = 0;
end
if ~exist('lae','var') | isempty(lae)
  lae.cpasma = 0;
end
if ~exist('gamma','var') | isempty(gamma)
  gamma = 0;
end
    
% It's important that the plasma response is a *continuous* function of the independent
% variables. In order to ensure that the plasma moves smoothly across the grid, cubic
% interpolation is used. It is also important that all derivatives are completely accurate.

% Variables for cubic interpolation on the grid
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % Used for interpolation
ir4 = (-1:2)'*nz; % For finding flux on horizontal grid lines
iz4 = (-1:2)';    % For finding flux on vertical grid lines
ir8 = [iz4-nz;iz4+nz]; % For finding derivatives w.r.t. r on vertical grid lines
iz8 = [ir4- 1;ir4+ 1]; % For finding derivatives w.r.t. z on horizontal grid lines
neighbors = reshape(iz4*ones(1,4)+ones(4,1)*ir4',1,16);

% pres and fpol^2/2-fvac^2/2 described by splined third degree polynomials
% where fvac = rzero*bzero
if isfield(config,'psikn')
  psikn = config.psikn; % psibar for knot junctions
elseif isfield(config,'nkn')
  psikn = linspace(0,1,config.nkn+1); % default is equally spaced values
end
if ~exist('psikn','var') | isempty(psikn)
  if exist('nkn','var') & ~isempty(nkn) & nkn > 0
    psikn = linspace(0,1,nkn+1);
  else
    psikn = [0 1];
  end
end
nkn = length(psikn) - 1; % Number of knots, i.e. number of third degree polynomials
if verbose > 0
  disp(['Spline knots are at psikn = ' num2str(psikn)])
end

if isfield(config,'no_edge_current')
  no_edge_current = config.no_edge_current;
end
if ~exist('no_edge_current','var') | isempty(no_edge_current)
  no_edge_current = 0;
end
if isfield(config,'no_edge_gradient')
  no_edge_gradient = config.no_edge_gradient;
end
if ~exist('no_edge_gradient','var') | isempty(no_edge_gradient)
  no_edge_gradient = 0;
end

if constraints & no_edge_current & no_edge_gradient & nkn == 1
  disp(['Warning gs_configure: ', ...
  'only 1 degree of freedom to vary pres and 1 to vary fpol', 10 ... 
  'since no_edge_current & no_edge_gradient & psikn = [0 1]'])
  if constraints
    disp(['Setting new value: psikn = [0 0.5 1];', 10 ...
    'since choice of constraints requires at least 3 degrees of freedom'])
    psikn = [0 0.5 1];
    nkn = 2;
  end
end

if ~exist('sp','var') | length(sp) ~= nkn+2
  sp = zeros(nkn+2,1);
end

if ~exist('sf','var') | length(sf) ~= nkn+2
  sf = zeros(nkn+2,1);
end

% The number of coefficients that can be chosen for a spline is nkn+2.
% The rest are determined by requiring continuous values and derivatives at knots
% Here spline vectors are
% sp(1:nkn) are values of c3 for each interval
% sp(nkn+1) is value of c2 for last interval
% sp(nkn+2) is value of c1 for last interval
c0 = zeros(nkn+1,nkn+2);
c1 = zeros(nkn+1,nkn+2);
c2 = zeros(nkn+1,nkn+2);
c3 = zeros(nkn+1,nkn+2);
c3(nkn,nkn) = 1;
c2(nkn,nkn+1) = 1;
c1(nkn,nkn+2) = 1;
c0(nkn,nkn+[0:2]) = -1;
for i = nkn-1:-1:1
  y = psikn(i+1);
  c3(i,i) = 1;
  % Match second derivative at knot i, i.e. 2*c2+6*c3*y
  dc3 = c3(i+1,:)-c3(i,:);
  c2(i,:) = c2(i+1,:) + 3*y*dc3;
  % Match first derivative at knot i, i.e. c1+2*c2*y+3*c3*y^2
  dc2 = c2(i+1,:)-c2(i,:);
  c1(i,:) = c1(i+1,:) + 2*y*dc2 + 3*y^2*dc3;
  % Match value at knot i, i.e. c0+c1*y+c2*y^2+c3*y^3
  dc1 = c1(i+1,:)-c1(i,:);
  c0(i,:) = c0(i+1,:) + y*dc1 + y^2*dc2 + y^3*dc3;
end

% The d-matrices are used for splines that are 0 at the axis rather than boundary
d0 = zeros(nkn+1,nkn+2); d0(nkn+1,nkn+2) = 1;
d1 = zeros(nkn+1,nkn+2);
d2 = zeros(nkn+1,nkn+2);
d3 = zeros(nkn+1,nkn+2);
c33 = inv([1 0 0; 0 2 0; psikn(2) psikn(2)^2 psikn(2)^3]);
d1(1,1:3) = c33(1,:);
d2(1,1:3) = c33(2,:);
d3(1,1:3) = c33(3,:);
for j = 2:nkn
  cs = inv([1 psikn(j  )   psikn(j  )^2   psikn(j  )^3; ...
            0 1          2*psikn(j  )   3*psikn(j  )^2; ...
	    0 0          2              6*psikn(j  )  ; ...
            1 psikn(j+1)   psikn(j+1)^2   psikn(j+1)^3]) * ...
         [0 1 psikn(j  )   psikn(j  )^2   psikn(j  )^3; ...
          0 0 1          2*psikn(j  )   3*psikn(j  )^2; ...
	  0 0 0          2              6*psikn(j  )  ; ...
	  1 0 0          0              0             ];
              
  d0(j,1:j+2) = [cs(1,2)*d0(j-1,1:j+1)+cs(1,3)*d1(j-1,1:j+1)+...
                 cs(1,4)*d2(j-1,1:j+1)+cs(1,5)*d3(j-1,1:j+1) cs(1,1)];
       
  d1(j,1:j+2) = [cs(2,2)*d0(j-1,1:j+1)+cs(2,3)*d1(j-1,1:j+1)+...
                 cs(2,4)*d2(j-1,1:j+1)+cs(2,5)*d3(j-1,1:j+1) cs(2,1)];
       
  d2(j,1:j+2) = [cs(3,2)*d0(j-1,1:j+1)+cs(3,3)*d1(j-1,1:j+1)+...
                 cs(3,4)*d2(j-1,1:j+1)+cs(3,5)*d3(j-1,1:j+1) cs(3,1)];
       
  d3(j,1:j+2) = [cs(4,2)*d0(j-1,1:j+1)+cs(4,3)*d1(j-1,1:j+1)+...
                 cs(4,4)*d2(j-1,1:j+1)+cs(4,5)*d3(j-1,1:j+1) cs(4,1)];
       
end

% Coefficients for the 3:rd degree polynomials are found by multiplying c0(j,:), etc by
% nkn+2 long spline vector, where j is index to a region, starting at axis
ns = 2*nkn+4; % total number of spline parameters

% Normalized flux for nr equally spaced values from axis to boundary
psibar = linspace(0,1,nr)';

% Index of spline region at each of the psibar values
iknotg = ones(nr,1);
for j = 2:nkn
  iknotg(psibar > psikn(j)) = j;
end

% Construct sp0, sf0, sg0 from template or default or keep existing
ll = logical(ones(nkn+2,1));
e1 = c1;
e2 = c2;
e3 = c3;
pp1 = eye(nkn+2);
if no_edge_gradient % 2*sp(nkn+1)+6*sp(nkn) = 0
  ll(nkn+1) = false;
  pp1(nkn+1,nkn) = -3;
  pp1(nkn+1,nkn+1) = 0;
end
pp2 = eye(nkn+2);
if no_edge_current % sp(nkn+2)+2*sp(nkn+1)+3*sp(nkn) = 0
  ll(nkn+2) = false;
  pp2(nkn+2,nkn) = -3;
  pp2(nkn+2,nkn+1) = -2;
  pp2(nkn+2,nkn+2) = 0;
end
pp = pp2*pp1(:,ll);
if isfield(config,'pres0')
  pres0 = config.pres0;
  psibar0 = linspace(0,1,length(pres0(:)))';
  pprime0 = (spline(psibar0,pres0,psibar0+1e-9) - ...
             spline(psibar0,pres0,psibar0-1e-9))/2e-9;
elseif isfield(config,'pprime0')
  pprime0 = config.pprime0(:);
end
if ~exist('pprime0','var') | isempty(pprime0)
  pprime0 = ones(nr,1);
end
psibar0 = linspace(0,1,length(pprime0))';
ikn = ones(length(pprime0),1);
for j = 2:nkn
  ikn(psibar0 > psikn(j)) = j;
end
pbnkn = psibar0*ones(1,nkn+2);
sp0 = pp*pinv([...
  c1(ikn,:)*pp + ...
  c2(ikn,:)*2.*pbnkn*pp + ...
  c3(ikn,:)*3.*pbnkn.^2*pp])*pprime0;
if isfield(config,'fpol0')
  fpol0 = config.fpol0;
  psibar0 = linspace(0,1,length(fpol0))';
  fprim0 = (spline(psibar0,fpol0,psibar0+1e-9) - ...
            spline(psibar0,fpol0,psibar0-1e-9))/2e-9;
  ffprim0 = fpol0.*fprim0;
elseif isfield(config,'ffprim0')
  ffprim0 = config.ffprim0(:);
end
if ~exist('ffprim0','var') | isempty(ffprim0)
  ffprim0 = ones(nr,1);
end
psibar0 = linspace(0,1,length(ffprim0))';
ikn = ones(length(ffprim0),1);
for j = 2:nkn
  ikn(psibar0 > psikn(j)) = j;
end
pbnkn = psibar0*ones(1,nkn+2);
sf0 = pp*pinv([...
  c1(ikn,:)*pp + ...
  c2(ikn,:)*2.*pbnkn*pp + ...
  c3(ikn,:)*3.*pbnkn.^2*pp])*ffprim0;
sg0 = pp*pinv([...
  c1(ikn,:)*pp + ...
  c2(ikn,:)*2.*pbnkn*pp + ...
  c3(ikn,:)*3.*pbnkn.^2*pp])*(1-psibar0);

if ~exist('er','var') | isempty(er)
  er = 1;
end

% For contouring
if isfield(config,'npola')
  npola = config.npola;
end
if ~exist('npola','var') | isempty(npola)
   npola = 2*nr-1;
end
if isfield(config,'psibarc')
  psibarc = config.psibarc;
end
if ~exist('psibarc','var') | isempty(psibarc)
   psibarc = psibar;
end
ncont = length(psibarc);
% Number of heating/current-drive actuators
if constraints
  ncd = 3;
else
  ncd = 2*nkn+4;
end
if isfield(config,'ncd')
  ncd = config.ncd;
end

% Read in conductor system
% Automatic calculation of mutuals may get implemented, if geometry is supplied
if isfield(config,'mcc')
  mcc = config.mcc;
end
if ~exist('mcc','var') | isempty(mcc)
  error(['Field mcc (and maybe others) are missing from config ', ...
         'but are required in the first call.'])
end
if isfield(config,'mcv')
  mcv = config.mcv;
end
if ~exist('mcv','var') | isempty(mcv)
  error(['Field mcv (and maybe others) are missing from config ', ...
         'but are required in the first call.'])
end
if isfield(config,'mvv')
  mvv = config.mvv;
end
if ~exist('mvv','var') | isempty(mvv)
  error(['Field mvv (and maybe others) are missing from config ', ...
         'but are required in the first call.'])
end
if isfield(config,'mpc')
  mpc = config.mpc;
end
if ~exist('mpc','var') | isempty(mpc)
  'room for code to calculate mpc because we do know rg and zg here'
end
if isfield(config,'mpv')
  mpv = config.mpv;
end
if ~exist('mpv','var') | isempty(mpv)
  'room for code to calculate mpv because we do know rg and zg here'
end
if isfield(config,'mpp')
  mpp = config.mpp;
end
if ~exist('mpp','var') | isempty(mpp)
  'room for code to calculate mpp because we do know rg and zg here'
end
nc = size(mcc,1);
nv = size(mvv,1);
if ~exist('need_mrmz')
  need_mrmz = false;
end
if isfield(config,'ress') % Resistances for all conductors
  ress = config.ress;
elseif isfield(config,'resc') & isfield(config,'resv')
  ress = [config.resc(:); config.resv(:)];
end
if ~exist('ress','var') | isempty(ress)
  error(['Conductor resistances resc and resv are missing from config ', ...
         'but are required in the first call.'])
end
rss = diag(ress);
mss = [mcc mcv; mcv' mvv];

% Number of inputs
nx = nc+nv+ns+1;
% Indices into x vector
indic = 1:nc;
indiv = indic(end)+(1:nv);
indis = [indic indiv];
indsp = indiv(end)+(1:nkn+2);
indsf = indsp(end)+(1:nkn+2);
inder = nx;

% if evolve_option is 0 (meaning don't evolve) then xc is in use, 
% its format decided by the flag constraints
if constraints == 0
  nxc = nx;
elseif constraints == 1
  nxc = nc+nv+4;
elseif constraints == 2
  nxc = nc+nv+4;
else
  nxc = nc+nv+2*nkn+5;
end
dxdxc = eye(nx,nxc); % Allocate memory


% Read the limiter from config. Tolerate many variations in the format
if size(config.limdata,2) < size(config.limdata,1)
  limdata = config.limdata;
  s1r = 'config.limdata(';
  s2r = ',1) = ';
  s1z = 'config.limdata(';
  s2z = ',2) = ';
else
  limdata = config.limdata';
  s1r = 'config.limdata(1,';
  s2r = ') = ';
  s1z = 'config.limdata(2,';
  s2z = ') = ';
end
if min(limdata(:,1)) < min(limdata(:,2))
  limdata = limdata(:,[2 1]);
  s1r = strrep(s1r,'1','2');
  s2r = strrep(s2r,'1','2');
  s1z = strrep(s1z,'2','1');
  s2z = strrep(s2z,'2','1');
end
if limdata(1,1) == limdata(end,1) & limdata(1,2) == limdata(end,2)
  Rlim = limdata(:,1);
  Zlim = limdata(:,2);
else
  Rlim = [limdata(:,1); limdata(1,1)];
  Zlim = [limdata(:,2); limdata(1,2)];
end
nlim = length(Rlim);
% The EAST limiter extends outside the grid as of March 2014, correct this kind of problem
ilim = Rlim < rg(2) | Rlim >= rg(nr-1) | Zlim < zg(2) | Zlim >= zg(nz-1);
if any(ilim) & 0
  disp('Warning gs_configure: Grid must cover limiter with 1 grid-step to spare!')
  disp('The following points from config.limdata will be automatically adjusted:')
  disp('index   old R, Z       after adjustment')
  h = figure;
  set(h,'NumberTitle','Off');
  set(h,'Name','Warning gs_configure: Grid must cover limiter with 1 grid-step to spare!');
  hold on
  for i = 1:nz
    plot([rg(1) rg(nr)],[zg(i) zg(i)],'color',[0.8 0.8 0.8])
  end
  for j = 1:nr
    plot([rg(j) rg(j)],[zg(1) zg(nz)],'color',[0.8 0.8 0.8])
  end
  h1(1) = plot(Rlim,Zlim,'b','linew',4);
  h1(2) = plot(Rlim(ilim),Zlim(ilim),'rx','linew',8);
  howtofix = 'Adjust config as follows to avoid this message in the future:';
  for i = 1:nlim
    r = Rlim(i);
    z = Zlim(i);
    if r < rg(2)
      Rlim(i) = rg(2);
      howtofix = char(howtofix,[s1r num2str(i) s2r num2str(ceil(rg(2)*1000)/1000) ';']);
    end
    if z < zg(2)
      Zlim(i) = zg(2);
      howtofix = char(howtofix,[s1z num2str(i) s2z num2str(ceil(zg(2)*1000)/1000) ';']);
    end
    if r >= rg(nr-1)
      Rlim(i) = rg(nr-1)-dr/1e9;
      howtofix = char(howtofix,[s1r num2str(i) s2r num2str(floor(rg(nr-1)*1000)/1000) ';']);
    end
    if z >= zg(nz-1)
      Zlim(i) = zg(nz-1)-dz/1e9;
      howtofix = char(howtofix,[s1z num2str(i) s2z num2str(floor(zg(nz-1)*1000)/1000) ';']);
    end
    if r ~= Rlim(i) | z ~= Zlim(i)
      fprintf('%3d   (%5.2f,%5.2f)     (%5.2f,%5.2f)\n',i,r,z,Rlim(i),Zlim(i))
    end
  end
  disp(howtofix)
  h1(3) = plot(Rlim,Zlim,'k--','linew',4);
  title('Warning gs\_configure: Grid must cover limiter!','color','red','fontw','bold')
  legend(h1(1:3),'old limiter','adjusted points','new limiter')
  axis([rg(1) 2*rg(nr)-rg(1) zg(1) zg(nz)])
  g = text(2*rg(nr)-rg(1),zg(ceil(nz/2)),' ', ...
    'FontSize',21,'HorizontalAlignment','Right','FontWeight','Bold');
  for i = 9:-1:2
    if ishandle(g)
      set(g,'String',['Figure closes   \newlinein ' num2str(i) ' seconds...  ']);
      pause(1)
    end
  end
  if ishandle(g)
    set(g,'String',['Figure closes   \newlinein 1 second...  ']);
    pause(1)
  end
  if ishandle(g)
    set(g,'String',['Figure about   \newlineto close...  ']);
  end
  drawnow
  if ishandle(h)
    close(h)
  end
end
ilim = logical(ones(1,nlim)); % Flag for whether to include point
% Remove points that lie on a straight line
for i = 1:nlim-2
  v1 = [Rlim(i+1)-Rlim(i  ), Zlim(i+1)-Zlim(i  )];
  v2 = [Rlim(i+2)-Rlim(i+1), Zlim(i+2)-Zlim(i+1)];
  dum = v1*v2'/norm(v1)/norm(v2);
  if abs(dum-1) < 1e-9
    ilim(i+1) = false;
  end
  if abs(dum+1) < 1e-9 % The limiter data for DIII-D has this error (as of 2014)
    ilim(i+2) = false;
  end
end
nlim2 = sum(ilim);

% Prepare isinpoly
isinpoly([],[],Rlim(ilim),Zlim(ilim))

% Create an additional limiter description (rl, zl) for finding if the plasma touches it
rvmin = zeros(nz,nr)+inf; % Lowest r that can be reached horizontally without hitting vessel
rvmax = zeros(nz,nr)-inf; % Highest r that can be reached horizontally without hitting vessel
zvmin = zeros(nz,nr)+inf; % Lowest z that can be reached vertically without hitting vessel
zvmax = zeros(nz,nr)-inf; % Highest z that can be reached vertically without hitting vessel
nl = nlim2; % nl will become the lengths of rl & zl
iRlim2 = (Rlim(ilim)-rg(1))/dr*2;
iZlim2 = (Zlim(ilim)-zg(1))/dz*2;
for j = 1:nlim2-1
  if iRlim2(j+1) > iRlim2(j)
    nl = nl + ceil(iRlim2(j+1))-floor(iRlim2(j))-1;
  elseif iRlim2(j+1) < iRlim2(j)
    nl = nl + ceil(iRlim2(j))-floor(iRlim2(j+1))-1;
  end
  if iZlim2(j+1) > iZlim2(j)
    nl = nl + ceil(iZlim2(j+1))-floor(iZlim2(j))-1;
  elseif iZlim2(j+1) < iZlim2(j)
    nl = nl + ceil(iZlim2(j))-floor(iZlim2(j+1))-1;
  end
end
rl = zeros(1,nl);
zl = zeros(1,nl);
% Now that memory has been allocated for rl, zl, fill in the values
rl(1) = Rlim(1);
zl(1) = Zlim(1);
nl = 1;
ir = iRlim2(1); % ir is in units of grid column number * 2
iz = iZlim2(1); % iz is in units of grid row number * 2
for j = 1:nlim2-1
  sz = sign(iZlim2(j+1)-iz); % At this point iz = iZlim(j)
  sr = sign(iRlim2(j+1)-ir); % At this point ir = iRlim(j)
  fz = 0; % To get while loop started
  fr = 0;
  while fr < 1 | fz < 1
    % Look for desired r values between Rlim(j) and Rlim(j+1)
    if sr > 0
      kr = floor(ir)+sr;
    elseif sr < 0
      kr = ceil(ir)+sr;
    else
      kr = 0;
    end
    if (kr - ir) * (kr - iRlim2(j+1)) < 0 & sr ~= 0
      fr = (kr-ir)/(iRlim2(j+1)-ir); % Grid line or cell edge located at this fraction of remaining distance to iRlim2(j+1)
    else
      fr = 2; % Any number > 1 will do
    end
    % Look for desired z values between Zlim(j) and Zlim(j+1)
    if sz > 0
      kz = floor(iz)+sz;
    elseif sz < 0
      kz = ceil(iz)+sz;
    else
      kz = 0;
    end
    if (kz - iz) * (kz - iZlim2(j+1)) < 0 & sz ~= 0
      fz = (kz-iz)/(iZlim2(j+1)-iz); % Grid line or cell edge located at this fraction of remaining distance to iZlim2(j+1)
    else
      fz = 2; % Any number > 1 will do
    end
    if fr < 1 & fr < fz % Moving to new grid line or cell edge at lower r if sr < 0 or higher r if sr > 0
      ir = kr;
      iz = iz + fr*(iZlim2(j+1)-iz);
      if mod(kr,2) == 0
	i = 1+round(ir/2);
	for k = 1+ceil(iz/2):nz
	  if isinf(zvmin(k,i))
	    zvmin(k,i) = iz*dz/2+zg(1);
	  else
	    zvmin(k,i) = max(zvmin(k,i),iz*dz/2+zg(1));
	  end
	end
	for k = 1:1+floor(iz/2)
	  if isinf(zvmax(k,i))
	    zvmax(k,i) = iz*dz/2+zg(1);
	  else
	    zvmax(k,i) = min(zvmax(k,i),iz*dz/2+zg(1));
	  end
	end
      end
    elseif fz < 1 & fz < fr
      iz = kz;
      ir = ir + fz*(iRlim2(j+1)-ir);
      if mod(kz,2) == 0
	i = 1+round(iz/2);
	for k = 1+ceil(ir/2):nr
	  if isinf(rvmin(i,k))
	    rvmin(i,k) = ir*dr/2+rg(1);
	  else
	    rvmin(i,k) = max(rvmin(i,k),ir*dr/2+rg(1));
	  end
	end
	for k = 1:1+floor(ir/2)
	  if isinf(rvmax(i,k))
	    rvmax(i,k) = ir*dr/2+rg(1);
	  else
	    rvmax(i,k) = min(rvmax(i,k),ir*dr/2+rg(1));
	  end
	end
      end
    else
      ir = iRlim2(j+1);
      iz = iZlim2(j+1);
    end
    nl = nl+1;
    rl(nl) = ir*dr/2+rg(1);
    zl(nl) = iz*dz/2+zg(1);
  end
end
% Remove any duplicates
il = [diff(rl).^2+diff(zl).^2 > (dr^2+dz^2)/1e6, 1==1];
rl = rl(il);
zl = zl(il);
nl = length(rl);

if 0 % Confirmation plots, useful if more debugging is needed above
  pb
  for j = 1:ngg
    h1 = plot(rgg(j),zgg(j),'go',rgg(j)+[0 0],[zvmin(j) zvmax(j)],'r-','linew',3);
    h2 = plot(rgg(j),zgg(j),'go',[rvmin(j) rvmax(j)],zgg(j)+[0 0],'r-','linew',3);
    drawnow
    delete(h1)
    delete(h2)
  end
end
% Prepare for calculating fluxes at rl, zl using the cubic interpolation method
iil = zeros(nl,16); % Will hold indices into the grid
wl  = zeros(nl,16); % Calculates flux at rl, zl
wlr = zeros(nl,16); % Calculates derivative in flux w.r.t. R
wlz = zeros(nl,16); % Calculates derivative in flux w.r.t. Z
wld = zeros(nl,16); % Calculates derivative in flux toward higher index
wlb = zeros(nl,16); % Calculates 2:nd derivative toward higher index
wlt = zeros(nl,16); % Calculates 3:rd derivative toward higher index
drl = diff([rl rl(2)]);
dzl = diff([zl zl(2)]);
dl = zeros(nl,1);
for i = 1:nl
  % Indices valid from i to i+1 (rl(i)+drl(i)/2,zl(i)+dzl(i)/2 is in this range)
  kr0 = min(nr-3,max(1,floor( (rl(i)+drl(i)/2 -rg(1))/dr)));
  kz1 = min(nz-2,max(2,floor( (zl(i)+dzl(i)/2 -zg(1))/dz)))+1;
  k = kr0*nz+kz1;
  iil(i,:) = k+neighbors;
  tr = (rl(i)-rgg(k))/dr;
  tz = (zl(i)-zgg(k))/dz;
  dl(i) = norm([drl(i) dzl(i)]);
  ur = drl(i)/dl(i);
  uz = dzl(i)/dl(i);
  wr0 = [1 tr tr^2 tr^3]*mx;
  wz0 = mx'*[1 tz tz^2 tz^3]';
  wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
  wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
  wr2 = [0 0 2 6*tr]*mx/dr^2;
  wz2 = mx'*[0 0 2 6*tz]'/dz^2;
  wr3 = [0 0 0 6]*mx/dr^3;
  wz3 = mx'*[0 0 0 6]'/dz^3;
  wl(i,:)  = reshape(wz0*wr0, 1, 16);
  wlr(i,:) = reshape(wz0*wr1, 1, 16);
  wlz(i,:) = reshape(wz1*wr0, 1, 16);
  wld(i,:) = reshape(wz0*wr1*ur + wz1*wr0*uz, 1, 16);
  wlb(i,:) = reshape(wz0*wr2*ur^2 + 2*wz1*wr1*ur*uz + wz2*wr0*uz^2, 1, 16);
  wlt(i,:) = reshape(wz0*wr3*ur^3 + 3*wz1*wr2*ur^2*uz + ...
                                    3*wz2*wr1*ur*uz^2 + wz3*wr0*uz^3, 1, 16);  
end

% Check if the left or right side is inside when walking along limiter
inside = inpolygon((rl(1:nl-1)+rl(2:nl))/2+diff(zl)/1e9, ...
                  (zl(1:nl-1)+zl(2:nl))/2-diff(rl)/1e9,rl,zl);
if sum(inside) == 0
  turnin = [0 -1;+1 0];
elseif sum(inside) == nl-1
  turnin = [0 +1;-1 0];
else
  if round(sum(inside)/(nl-1)) == 0 % Hopefully almost everyone is 0
    turnin = [0 -1;+1 0];
    disp(['Warning gs_configure: Limiter strange in region ', ...
      num2str(min(rl(inside))),' < R < ' num2str(max(zl(inside))) ', ', ...
      num2str(min(zl(inside))) ' < Z < ' num2str(max(zl(inside)))])
  else % Hopefully almost everyone is 1
    turnin = [0 +1;-1 0];
    disp(['Warning gs_configure: Limiter strange in region ', ...
      num2str(min(rl(~inside))),' < R < ' num2str(max(rl(~inside))) ', ', ...
      num2str(min(zl(~inside))) ' < Z < ' num2str(max(zl(~inside)))])
  end
end
anglel = angle(drl.*diff(rl([nl-1 1:nl]))+dzl.*diff(zl([nl-1 1:nl])) + ...
  1i*(turnin(1,2)*dzl.*diff(rl([nl-1 1:nl]))+...
      turnin(2,1)*drl.*diff(zl([nl-1 1:nl]))));
concavel = anglel > 0.001;
convexl = anglel < -0.001;

% ilimgg is flags for the grid: -1 = outside limiter, 0 = inside limiter, otherwise at the limiter
% Limiter considered "at grid element j" if within rgg(j)+[0 dr], zgg(j)+[0 dz].
ilimgg = zeros(nz,nr);
ilimgg( 1,:) = -1; % The edge of the grid assumed to be outside the limiter
ilimgg(nz,:) = -1;
ilimgg(:, 1) = -1;
ilimgg(:,nr) = -1;
% For grid elements at the limiter, use index of nearest Rlim, Zlim as flag
for j = 1:nlim-1
  kra = (Rlim(j+0)-rg(1))/dr;
  kza = (Zlim(j+0)-zg(1))/dz;
  krb = (Rlim(j+1)-rg(1))/dr;
  kzb = (Zlim(j+1)-zg(1))/dz;
  if krb > kra
    sr = 1;
  else
    sr = 0;
  end
  if kzb > kza
    sz = 1;
  else
    sz = 0;
  end

  kr = kra;
  kz = kza;
  while abs(kr - krb) > 1e-9 | abs(kz - kzb) > 1e-9
    if floor(kr) ~= floor(krb)
      fr = (floor(kr)+sr-kr)/(krb-kr); % Fraction of remaining distance to get to new column
    else
      fr = 1;
    end
    if floor(kz) ~= floor(kzb)
      fz = (floor(kz)+sz-kz)/(kzb-kz); % Fraction of remaining distance to get to new row
    else
      fz = 1;
    end
    if fr < fz & fr < 1
      kr = kr+(fr+1e-9*(fr==0))*(krb-kr);
      kz = kz+(fr+1e-9*(fr==0))*(kzb-kz);
    elseif fz < 1
      kr = kr+(fz+1e-9*(fz==0))*(krb-kr);
      kz = kz+(fz+1e-9*(fz==0))*(kzb-kz);
    else
      kr = krb;
      kz = kzb;
    end
    kr0 = min(nr-1,max(0,floor(kr)));
    kz1 = min(nz,max(1,ceil(kz)));
    ilimgg(kz1+kr0*nz) = j;
  end
end
not_found_all_vac_gg = true;
while not_found_all_vac_gg
  not_found_all_vac_gg = false;
  for j = 2:nz-1
    for k = 2:nr-1
      if ilimgg(j,k) == 0 & (ilimgg(j-1,k)==-1 | ilimgg(j+1,k)==-1 | ilimgg(j,k-1)==-1 | ilimgg(j,k+1)==-1)
	ilimgg(j,k) = -1;
	not_found_all_vac_gg = true;
      end
    end
  end      
end
for i = 1:nl
  ilimgg(iil(i,6)) = i;
end

dpsizrpladpcurrt = zeros(ngg);
for j = 1:nz
  izshift = (1+abs(j-(1:nz))-(1:nz))'*ones(1,nr);
  dpsizrpladpcurrt(:,j+(0:nr-1)*nz) = mpp((1:ngg)'+izshift(:),:);
end
M = zeros(ngg);

% Allocate memory for boundary points
nbbbs_max = 2*nr+2*nz;
rbbbs = zeros(nbbbs_max,1);
zbbbs = zeros(nbbbs_max,1);
nbbbs = 0;
r1 = zeros(nbbbs_max,1); % Coefficients for boundary interpolation
r2 = zeros(nbbbs_max,1); % Coefficients for boundary interpolation
r3 = zeros(nbbbs_max,1); % Coefficients for boundary interpolation
z1 = zeros(nbbbs_max,1); % Coefficients for boundary interpolation
z2 = zeros(nbbbs_max,1); % Coefficients for boundary interpolation
z3 = zeros(nbbbs_max,1); % Coefficients for boundary interpolation
thbbbs = zeros(nbbbs_max,1); % Angles
rhobbbs = zeros(nbbbs_max,1); % distances to the axis

% circle is a flag to use circular plasma model
if ~exist('circle','var') | length(circle) ~= 1
  circle = logical(0);
end
tiny_plasma = 0;

% Variables used to trace the boundary
a = zeros(4*nbbbs_max,1); % Coefficients for cubic term
b = zeros(4*nbbbs_max,1); % Coefficients for quadratic term
c = zeros(4*nbbbs_max,1); % Coefficients for linear term
d = zeros(4*nbbbs_max,1); % The rest
aa = zeros(4*nbbbs_max,1); % Calculation step toward solving cubic equation
bb = zeros(4*nbbbs_max,1); % Calculation step toward solving cubic equation
q1 = zeros(4*nbbbs_max,1); % Calculation step toward solving cubic equation
q2 = zeros(4*nbbbs_max,1); % Calculation step toward solving cubic equation
x1 = zeros(4*nbbbs_max,1); % First of 3 solutions
x2 = zeros(4*nbbbs_max,1); % Second of 3 solutions
x3 = zeros(4*nbbbs_max,1); % Third of 3 solutions
xx = zeros(4*nbbbs_max,1); % Chosen solution (value between 0 and 1)
re = zeros(4*nbbbs_max,1);
ze = zeros(4*nbbbs_max,1);
fe = zeros(4*nbbbs_max,1);
gbbbs = zeros(nbbbs_max,1); % boundary point flags
ibbbs = zeros(nbbbs_max,1); % Indices to grid point below and inside each boundary point
kk = zeros(nbbbs_max,1); % Indices of sorted thbbbs
dpsibbbsdr = zeros(nbbbs_max,1); % d(psi)/dR for bbbs points
dpsibbbsdz = zeros(nbbbs_max,1); % d(psi)/dZ for bbbs points
d2psibbbsdrdr = zeros(nbbbs_max,1); % d2(psi)/dR2 for bbbs points
d2psibbbsdrdz = zeros(nbbbs_max,1); % d2(psi)/dRdZ for bbbs points
d2psibbbsdzdz = zeros(nbbbs_max,1); % d2(psi)/dZ2 for bbbs points
drbbbs = zeros(nbbbs_max,1);
dzbbbs = zeros(nbbbs_max,1);
dzbbbs_max = 2*dz; % Used by gs_find_bdef, biggest possible dzbbbs since last analysis
redge = zeros(nbbbs_max,1);
zedge = zeros(nbbbs_max,1);
fedge = zeros(nbbbs_max,1);
gedge = zeros(nbbbs_max,1);
iedge = zeros(nbbbs_max,1);
xedge = zeros(nbbbs_max,1);
thedge = zeros(nbbbs_max,1);
rhoedge = zeros(nbbbs_max,1);
id4 = ones(nbbbs_max,4);
id8 = ones(nbbbs_max,8);
id16 = ones(nbbbs_max,16);
wbbbs = zeros(nbbbs_max,4); % Weights to calculate flux, used with ir4 and iz4
wds = zeros(nbbbs_max,4);
wbs = zeros(nbbbs_max,4);
wrs = zeros(nbbbs_max,4);
wzs = zeros(nbbbs_max,4);
dpsid4 = zeros(nbbbs_max,1);
dpsid8 = zeros(nbbbs_max,1);
dpsib4 = zeros(nbbbs_max,1);
dpsib8 = zeros(nbbbs_max,1);
dpsib16 = zeros(nbbbs_max,1);
dCldrbbbs = zeros(1,nbbbs_max-1);
dCldzbbbs = zeros(1,nbbbs_max-1);

drb3dpsia = zeros(nbbbs_max,16);
drb3dpsib = zeros(nbbbs_max,16);
drb3dpsip = zeros(nbbbs_max,16);
dzb3dpsia = zeros(nbbbs_max,16);
dzb3dpsib = zeros(nbbbs_max,16);
dzb3dpsip = zeros(nbbbs_max,16);

drbbbsdx   = zeros(nbbbs_max,nx);
dzbbbsdx   = zeros(nbbbs_max,nx);
drbdx   = zeros(nbbbs_max,nx);
dzbdx   = zeros(nbbbs_max,nx);
drsurfdx = zeros(1,nx);
daminordx = zeros(1,nx);

if ~exist('rbdef','var') | isempty(rbdef)
  rbdef = nan;
end
if ~exist('zbdef','var') | isempty(zbdef)
  zbdef = nan;
end
if ~exist('psizr','var') | isempty(psizr)
  psizr = zeros(nz,nr);
end
[n1,n2] = size(psizr);
if n1 ~= nz | n2 ~= nr
  psizr = zeros(nz,nr);  
end
if ~exist('maxis','var') | isempty(maxis)
  maxis = false;
end

% Circle - model for small plasmas
if ~exist('circle_model','var') | isempty(circle_model)
  circle_model = 0; % index to a circle model, 0 = no circle
end
% Flag for whether to use circle_model at all
if isfield(config,'use_circle_model')
  use_circle_model = config.use_circle_model;
end
if ~exist('use_circle_model','var') | isempty(use_circle_model)
  use_circle_model = 0;
end

% Simulation switches from circle to grid model when Atot/Ag > ncp2gr
if isfield(config,'ncp2gr')
  ncp2gr = config.ncp2gr;
end
if ~exist('ncp2gr','var') | isempty(ncp2gr)
  ncp2gr = 100;
end

% Simulation switches from grid to circle model when Atot/Ag < ngr2cp
if isfield(config,'ngr2cp')
  ngr2cp = config.ngr2cp;
end
if ~exist('ngr2cp','var') | isempty(ngr2cp)
  ngr2cp = 75;
end

% Halo model - model for current on open field lines
if ~exist('halo','var') | isempty(halo)
  halo = 0; % index to a halo model, 0 = no halo
end
geohzr = zeros(nz,nr);
for i = 1:nz
  for j = 1:nr
    if ilimgg(i,j) > -1
      geohzr = min((rl-rg(j)).^2+(zl-zg(i)).^2)/rg(j);
    end
  end
end
ih = zeros(nz,nr); % Halo current or current on open field lines [A]
ihalo = 0; % sum(ihzr(:))

dpsizr = zeros(nz,nr); % Total flux [Wb]
dpsizrx = zeros(nz,nr); % Flux change dummy variable [Wb]
dpsizr_err = zeros(nz,nr); % Total flux change for error correction [Wb]
iused = logical(zeros(nz,nr)); % True for points that are used in equations
iplasma = logical(zeros(nz,nr)); % True for points covered by GS plasma
dAcelldx = zeros(ngg,nx); % How covered area of cells respond to x
drmaxisdx = zeros(1,nx); % How radius of axis responds to x
dzmaxisdx = zeros(1,nx); % How height of axis responds to x

% Load alternative default for plotting if supplied
if isfield(config,'plotit')
  plotit = config.plotit;
  default_plotit = config.plotit;
end
if ~exist('default_plotit','var') | isempty(default_plotit)
  default_plotit = false;
end
if ~exist('plotit','var') | isempty(plotit)
  plotit = default_plotit;
end
if isfield(config,'plot_settings')
  plot_settings = config.plot_settings;
end
rtplot_handles.refresh = true;
plot_settings.refresh = true;
if isfield(config,'ecdata')
  ecdata = config.ecdata;
end
if isfield(config,'fcdata')
  fcdata = config.fcdata;
end
if isfield(config,'vvdata')
  vvdata = config.vvdata;
end

% Load alternative default for converging if supplied
if isfield(config,'converge')
  default_converge = config.converge;
  converge = config.converge;
end
if ~exist('default_converge','var') | isempty(default_converge)
  default_converge = false;
end
if ~exist('converge','var') | isempty(converge)
  converge = default_converge;
end

% Flag for artificial stabilization
if isfield(config,'stabilize')
  default_stabilize = config.stabilize;
  stabilize = default_stabilize;
end
if ~exist('default_stabilize','var') | isempty(default_stabilize)
  default_stabilize = 0;
end
if ~exist('stabilize','var') | isempty(stabilize)
  stabilize = default_stabilize;
end
% Flag for whether stabilization was actually done
if ~exist('stabilized','var') | isempty(stabilized)
  stabilized = false;
end

% Load alternative default for time if supplied
if isfield(config,'time')
  default_time = config.time;
end
if ~exist('default_time','var') | isempty(default_time)
  default_time = 0;
end
if ~exist('time','var') | isempty(time)
  time = default_time;
end

% For time-dependent & conditional plotting
if isfield(config,'plot_times')
  plot_times = config.plot_times;
end
if ~exist('plot_times','var') | isempty(plot_times)
  plot_times = inf; 
end
if ~exist('plot_counter','var') | isempty(plot_counter)
  plot_counter = 1; 
end
if isfield(config,'dtplot')
  dtplot = config.dtplot;
end
if ~exist('dtplot','var') | isempty(dtplot)
  dtplot = 0; 
end
if dtplot > 0
  dum = max(0,-floor(log10(dtplot)));
  title_time_format = ['time = %0.' num2str(dum) 'f'];
end
if ~exist('time_for_next_plot','var') | isempty(time_for_next_plot)
 time_for_next_plot  = -inf; 
end

if isfield(config,'plot_if_new_response')
  plot_if_new_response = config.plot_if_new_response;
end
if ~exist('plot_if_new_response','var') | isempty(plot_if_new_response)
  plot_if_new_response = false;
end

% Configure gs_dynamics
if isfield(config,'evolve_option')
  evolve_option = config.evolve_option;
end
if ~exist('evolve_option','var') | isempty(evolve_option)
  evolve_option = 0;
end
if isfield(config,'plares')
  if isstruct(config.plares)
    plares = config.plares;
  else
    plares.values = config.plares;
    plares.times = linspace(0,10,length(plares.values));
  end
end
if ~exist('plares','var') | isempty(plares)
  plares.values = [1 1]*1e-12;
  plares.times = [0 10];
end

% Configure gs_update_eq
if isfield(config,'dpsibar_limit_for_linear')
  dpsibar_limit_for_linear = config.dpsibar_limit_for_linear;
end
if ~exist('dpsibar_limit_for_linear','var') | isempty(dpsibar_limit_for_linear)
  dpsibar_limit_for_linear = 2e-3;
end

calculate_helical_voltage = evolve_option == 2; % May change down below
calculate_profile_responses = 0;
calculate_profiles = 0; % May change down below



if isfield(config,'limits')
  limits = config.limits;
end
if ~exist('limits','var')
  limits = [];
end
if ~isfield(limits,'ic')
  if strcmp(upper(tokamak),'D3D')
     % Ampere limits for E&F coils
    limits.ic = [2e3 2e3 ...
      59 59 59 59 54 124 124 59 89 ...
      59 59 59 59 54 124 124 59 89]'*[-100 100];
  end
  if strcmp(upper(tokamak),'EAST')
    % Limits for PF1 to PF14 and IC1, IC2, taken from PCS settings:
    limits.ic = [12.8*ones(6,1);11.6*ones(4,1);10.2*ones(4,1);5;5]*...
      [-1000 1000];
  end
  if strcmp(upper(tokamak),'KSTAR')
    % Limits taken from PCS settings for shot 5516 for coils in the order:
    % PF1-7U, PF1-7L, IVCU, IRCU, IVCL, IRCL
    limits.ic = [7.5 8.5 9 9 7 3.5 3.5 7.5 8.5 9 9 7 3.5 3.5 3 8 3 8]'*...
      [-1000 1000];
  end
  if strcmp(upper(tokamak),'ITER')
    % Limits taken from ITER_data_2010-v3.3.xls
    % PF1 PF2 PF3 PF4 PF5 PF6 CS3L CS2L CS1L CS1U CS2U CS3U VS3U VS3L
    % The high B limits, 0.4 K sub cooling PF6
    limits.bc = [65 50 50 50 60 70 130 130 130 130 130 130 nan nan]'*[-0.10 0.10];
    limits.ic = [41 50 50 50 33 41  40  40  40  40  40  40  10  10]'*[-1000 1000];
    % The low B limits, 0.4 K sub cooling PF6
    limits.bc = [64 48 48 48 57 68 126 126 126 126 126 126 nan nan]'*[-0.10 0.10];
    limits.ic = [48 55 55 55 52 52  45  45  45  45  45  45 10  10]'*[-1000 1000];
    % The high B limits
    limits.bc = [65 50 50 50 60 65 130 130 130 130 130 130 nan nan]'*[-0.10 0.10];
    limits.ic = [41 50 50 50 33 41  40  40  40  40  40  40 10  10]'*[-1000 1000];
    % The low B limits
    limits.bc = [64 48 48 48 57 64 126 126 126 126 126 126 nan nan]'*[-0.10 0.10];
    limits.ic = [48 55 55 55 52 48  45  45  45  45  45  45 10  10]'*[-1000 1000];
  end
end


% icci makes TokSys terminal currents ic from circuit currents ci
if isfield(config,'icci')
  if size(config.icci,1) > nc
    disp(['Warning gs_configure: Projection matrix icci has more coils ' ...
      'than available and will be cut after coil ' num2str(nc)])
  end
  if isfield(config,'cccirc')
    disp(['Warning gs_configure: Circuits were specified with icci ' ...
      'and cccirc, icci will be used, cccirc will be ignored'])
  end
  if isfield(config,'buscode')
    disp(['Warning gs_configure: Circuits were specified with icci ' ...
      'and buscode, icci will be used, buscode will be ignored'])
  end
  nci = size(config.icci,2);
  n = min(nc,size(config.icci,1));
  icci = [config.icci(1:n,:); zeros(nc-n,nci)];
elseif isfield(config,'cccirc')
  nci = max(abs(config.cccirc));
  icci = zeros(nc,nci);
  if length(config.cccirc(:)) ~= nc
    error(['The length of cccirc differs from nc. ' ...
      'Check the connections in cccirc.' 10 ...
      'The tokamak is ' upper(tokamak) ', and nc = ' num2str(nc) ...
      ' but length(config.cccirc) = ' num2str(length(config.cccirc))])
  end
  for j = 1:nc
    k = abs(config.cccirc(j));
    s = sign(config.cccirc(j));
    icci(j,k) = s;
  end
  if isfield(config,'buscode')
    disp(['Warning gs_configure: Circuits were specified with both ' ...
      'cccirc and buscode, cccirc will be used, buscode will be ignored'])
  end
elseif isfield(config,'buscode')
  icci = eye(nc);
  ind = find(config.buscode);
  icci(min(ind),ind) = -1;
  icci = icci(:,[1:ind(1)-1 ind(1)+1:nc]);
end
if ~exist('icci','var') | isempty(icci)
  icci = eye(nc);
end
if size(icci,1) ~= nc
  disp(['Warning gs_configure: The coil connection matrix, icci is for nc = ', ...
    num2str(size(icci,1)), ', BUT nc (derived from mcc) = ', num2str(nc)])
  disp('The coil connection matrix will be set to eye(nc)')
  icci = eye(nc);
end
% netlist can be used to describe the circuits
if isfield(config,'netlist')
  netlist = config.netlist;
  in = logical(ones(size(netlist.names,1),1));
  for i = 1:length(in)
    if strcmp(deblank(netlist.names(i,:)),'MIp')
      in(i) = false;
    end
  end
  netlist.names  = netlist.names(in,:);
  netlist.N1     = netlist.N1(in,:);
  netlist.N2     = netlist.N2(in,:);
  netlist.values = netlist.values(in,:);
  if isfield(config,'netlist_currents')
    netlist_currents = config.netlist_currents;
  else
    netlist_currents = '';
  end
  if isfield(config,'netlist_voltages')
    netlist_voltages = config.netlist_voltages;
  else
    netlist_voltages = [];
  end
  n = sum(netlist.names(:,1) == 'M');
  netlist_rxx = zeros(n);
  for i = 1:min(n,size(rss,1))
    for j = 1:min(n,size(rss,2))
      netlist_rxx(i,j) = rss(i,j);
    end
  end
  if exist('Rextra','var')
    for i = 1:min(n,size(Rextra,1))
      for j = 1:min(n,size(Rextra,2))
        netlist_rxx(i,j) = netlist_rxx(i,j)+Rextra(i,j);
      end
    end
  end
  netlist_model = model_from_netlist(...
    netlist,eye(n),netlist_rxx,netlist_currents,netlist_voltages);
  netlist_model.Lckt_extra = ...
    netlist_model.Mhat - netlist_model.Pxx'*netlist_model.Pxx;
else
  netlist_model = [];
end
nci = size(icci,2);
picci = pinv(icci);

if evolve_option == 99
  disp('evolve_option 99 has been changed to 9')
  evolve_option = 9
end
if evolve_option == 1 | evolve_option == 9
  % conductors & scalar plasma
  if isempty(netlist_model)
    Pxx = [icci zeros(nc,nv+3); zeros(nv+3,nci) eye(nv+3)];
    Lextra = zeros(nc+nv+3); % Extra inductances in coils and vessel
    if isfield(config,'Lextra')
      [n1,n2] = size(config.Lextra);
      Lextra(1:n1,1:n2) = config.Lextra;
    end
    rxx = [rss zeros(nc+nv,3); zeros(3,nc+nv+3)];
    Rextra = zeros(nc+nv+3); % Extra resistances in coils and vessel
    if isfield(config,'Rextra')
      [n1,n2] = size(config.Rextra);
      Rextra(1:n1,1:n2) = config.Rextra;
    end
    Rhat = (Pxx'*(rxx+Rextra)*Pxx);
    Vhat = eye(size(Pxx,2));
  else % It is assumed here that netlist_model.Pxx projects onto [ic;iv]
    [n1,n2] = size(netlist_model.Pxx);
    Pxx = [netlist_model.Pxx zeros(n1,3); zeros(3,n2) eye(3)];
    Lextra = zeros(nc+nv+3); % Extra inductances in coils and vessel
    [n1,n2] = size(netlist_model.Vhat);
    Vhat = [netlist_model.Vhat zeros(n1,3); zeros(3,n2) eye(3)];
    [n1,n2] = size(netlist_model.Rhat);
    Rhat = [netlist_model.Rhat zeros(n1,3); zeros(3,n2) eye(3)];
  end
  Pxxi = pinv(Pxx);
elseif evolve_option == 2
  %calculate_profiles = 1;
  % conductors & 1-D resistivity (nx-1 states, omitting er)
  if isfield(config,'eta_vs_rhot')
    eta_vs_rhot = config.eta_vs_rhot(:);
  end
  if ~exist('eta_vs_rhot','var') | isempty(eta_vs_rhot)
%    disp(['Warning gs_configure: evolve_option = 2 but ', ...
%      'no resistivity profile given, setting eta_vs_rhot to a flat 1e-8...'])
%    eta_vs_rhot = 1e-8*ones(nr,1);
    eta_vs_rhot = zeros(nr,1);
  end
  Pxx = [icci zeros(nc,nx-1-nc); zeros(nx-1-nc,nci) eye(nx-1-nc)];
  Lextra = zeros(nx-1); % Extra inductances in coils and vessel
  rxx = [rss zeros(nc+nv,nx-1-nc-nv); zeros(nx-1-nc-nv,nx-1)];
  Rextra = zeros(nx-1); % Extra resistances in coils and vessel
else
  Pxx = [icci zeros(nc,nv); zeros(nv,nci) eye(nv)];
  Lextra = zeros(nc+nv,nx); % Extra inductances in coils and vessel
  Rextra = zeros(nc+nv); % Extra resistances in coils and vessel
end
nxs = size(Pxx,2);
if isfield(config,'Lextra')
  [n1,n2] = size(config.Lextra);
  Lextra(1:n1,1:n2) = config.Lextra;
end
if isfield(config,'Rextra')
  [n1,n2] = size(config.Rextra);
  Rextra(1:n1,1:n2) = config.Rextra;
end
if isfield(config,'nps')
  nps = config.nps;
end
if ~exist('nps','var') | isempty(nps)
  nps = sum(any(Pxx(1:nc,:)));
end


% Correct size for all variables
if ~exist('psibarzr','var') || size(psibarzr,1) ~= nz || size(psibarzr,2) ~= nr
  psibarzr = zeros(nz,nr);
end
if ~exist('dpsimagdx','var') || size(dpsimagdx,1) ~= 1 || size(dpsimagdx,2) ~= nx
  dpsimagdx = zeros(1,nx);
end
if ~exist('dpsibrydx','var') || size(dpsibrydx,1) ~= 1 || size(dpsibrydx,2) ~= nx
  dpsibrydx = zeros(1,nx);
end
if ~exist('dpsizrdx','var') || size(dpsizrdx,1) ~= ngg || size(dpsizrdx,2) ~= nx
  dpsizrdx = zeros(ngg,nx);
end
if ~exist('wa','var') || size(wa,1) ~= 1 || size(wa,2) ~= 16
  wa = zeros(1,16);
end
if ~exist('wb','var') || size(wb,1) ~= 1 || size(wb,2) ~= 16
  wb = zeros(1,16);
end
if ~exist('iia','var') || size(iia,1) ~= 16 || size(iia,2) ~= 1
  iia = ones(16,1);
end
if ~exist('iib','var') || size(iib,1) ~= 16 || size(iib,2) ~= 1
  iib = ones(16,1);
end
nys = nc+nv;
if ~exist('ys','var') || size(ys,1) ~= nys || size(ys,2) ~= 1
  ys = zeros(nys,1);
end
if ~exist('dysdx','var') || size(dysdx,1) ~= nys || size(dysdx,2) ~= nx
  dysdx = zeros(nys,nx);
end


% Configure diagnostic points rdp, zdp
if isfield(config,'rdp') & isfield(config,'zdp')
  lrdp = length(config.rdp(:));
  lzdp = length(config.zdp(:));
  ndp = max(lrdp, lzdp)*(min(lrdp, lzdp) > 0);
  if lrdp == ndp
    rdp = config.rdp(:);
  else
    rdp = config.rdp(1)+zeros(ndp,1);
  end
  if lzdp == ndp
    zdp = config.zdp(:);
  else
    zdp = config.zdp(1)+zeros(ndp,1);
  end
  if lrdp > 1 & lzdp > 1 & lrdp ~= lzdp
    disp('Warning gs_configure: The sizes of diagnostic points rdp and zdp differ!')
    disp('They need to be two equally sized arrays or one scalar and one array.')
  end
else
  rdp = [];
  zdp = [];
end
ndp = length(rdp);
mdc = zeros(ndp,nc);
mdv = zeros(ndp,nv);
mpd = zeros(ngg,ndp);
grdc = zeros(ndp,nc);
grdv = zeros(ndp,nv);
grpd = zeros(ngg,ndp);
gzdc = zeros(ndp,nc);
gzdv = zeros(ndp,nv);
gzpd = zeros(ngg,ndp);
idp =  rdp > rg(1) & rdp < rg(nr) & zdp > zg(1) & zdp < zg(nz);
for j = 1:nc
  [mdc(idp,j) grdc(idp,j) gzdc(idp,j)] = gs_interp2(rg,zg,mpc(:,j),rdp(idp),zdp(idp));
end
for j = 1:nv
  [mdv(idp,j) grdv(idp,j) gzdv(idp,j)] = gs_interp2(rg,zg,mpv(:,j),rdp(idp),zdp(idp));
end
for j = 1:ngg
  [mpd(j,idp) grpd(j,idp) gzpd(j,idp)] = gs_interp2(rg,zg,dpsizrpladpcurrt(:,j),rdp(idp),zdp(idp));
end

% Number of times that a plasma response has been calculated
response_count = 0;

% Configure diagnostic vector of x-points
if isfield(config,'nxpoints')
  nxpoints = config.nxpoints;
end
if ~exist('nxpoints','var') | isempty(nxpoints)
  nxpoints = 2;
end
if nxpoints < 1
  nxpoints = 1;
end
if isfield(config,'xpointorder')
  xpointorder = config.xpointorder;
end
if ~exist('xpointorder','var') | isempty(xpointorder)
  xpointorder = 1;
end
if xpointorder < 1 || xpointorder > 4
  disp('Warning gs_configure: xpointorder must be between 1 and 4, setting it to 1.')
  xpointorder = 1;
end

% Configure gaps, defined by distance from r,z along gr,gz to separatrix
if isfield(config,'gapspec')
  gapspec = config.gapspec; % gapspec = [r z gr gz]
end
if ~exist('gapspec','var') | isempty(gapspec)
  gapspec = zeros(0,4);
end
ngaps = size(gapspec,1);

% Configure outputs
if isfield(config,'outputs')
  outputs = config.outputs;
end
if ~exist('outputs','var') | isempty(outputs)
  outputs = '';
end
outputs(outputs == 0) = 32;
% Resetting of index_in_y is needed because this may be a reconfiguration
index_in_y = [];
ny  = 0; % Will hold number of outputs = length(y)
ds = []; % Will hold descriptions for the fields of index_in_y
iy = []; % Indices to continuous outputs with smoothing between linearizations
zcount = 0; % Number of entries with zero-signals
usig = 0; % Number of user-defined signals
user_signal = []; % Descriptions of user-defined signals
for i = 1:size(outputs,1)
  str = deblank(outputs(i,:));
  if isfield(index_in_y,str)
    disp(['Warning in gs_configure: output ', str, ...
     ' appears more than once in outputs, only first entry used.'])
  elseif strcmp(str,'xs')
    ds.xs = 'state vector';
    index_in_y.xs = ny+(1:nxs);
    ny = ny+nxs;
  elseif strcmp(str,'xsdotB2')
    ds.xsdotB2 = 'time-derivative of state vector w.r.t. to vcd';
    index_in_y.xsdotB2 = ny+(1:nxs);
    ny = ny+nxs;
  elseif strcmp(str,'psizr')
    ds.psizr = 'Poloidal flux on the grid';
    index_in_y.psizr = ny+(1:ngg);
    iy = [iy ny+(1:ngg)];
    ny = ny+ngg;
  elseif strcmp(str,'ic')
    ds.ic = 'Coil currents, such that psizr_app = mpc*ic';
    index_in_y.ic = ny+(1:nc);
    iy = [iy ny+(1:nc)];
    ny = ny+nc;
  elseif strcmp(str,'icdot')
    ds.icdot = 'rate of change of coil currents [A/s]';
    index_in_y.icdot = ny+(1:nc);
    ny = ny+nc;
  elseif strcmp(str,'icdotA')
    ds.icdotA = 'rate of change of coil currents driven passively';
    index_in_y.icdotA = ny+(1:nc);
    ny = ny+nc;
  elseif strcmp(str,'icdotB')
    ds.icdotB = 'rate of change of coil currents driven by actuators';
    index_in_y.icdotB = ny+(1:nc);
    ny = ny+nc;
  elseif strcmp(str,'icdotB2')
    ds.icdotB2 = 'rate of change of coil currents driven by H&CD';
    index_in_y.icdotB2 = ny+(1:nc);
    ny = ny+nc;
  elseif strcmp(str,'dicdotdvps')
    ds.dicdotdvps = 'icdot response to circuit voltage from power supplies';
    index_in_y.dicdotdvps = ny+(1:nc)'*ones(1,nps)+ones(nc,1)*(0:nps-1)*nc;
    ny = ny+nc*nps;
  elseif strcmp(str,'dicdotdvcd')
    ds.dicdotdvcd = 'icdot response to circuit voltage from current-drive';
    index_in_y.dicdotdvcd = ny+(1:nc)'*ones(1,ncd)+ones(nc,1)*(0:ncd-1)*nc;
    ny = ny+nc*ncd;
  elseif strcmp(str,'iv')
    ds.iv = 'Vessel currents, such that psizr_app = mpv*iv';
    index_in_y.iv = ny+(1:nv);
    iy = [iy ny+(1:nv)];
    ny = ny+nv;
  elseif strcmp(str,'sp')
    ds.sp = ['x = linspace(0,1,nr)', 10, ...
             'pres = c0+c1*x+c2*x^2+c3*x^3', 10, ...
             'sp(1:nkn) are c3 for spline intervals', 10, ...
             'sp(nkn+1) is c2 for last interval', 10, ...
	     'sp(nkn+2) is c1 for last interval'];
    ds.sp = 'Spline vector for pres';
    index_in_y.sp = ny+(1:nkn+2);
    iy = [iy ny+(1:nkn+2)];
    ny = ny+nkn+2;
  elseif strcmp(str,'sf')
    ds.sf = ['x = linspace(0,1,nr)', 10, ...
             'fvac = rzero*bzero', 10, ...
             'f^2/2-fvac^2/2 = c0+c1*x+c2*x^2+c3*x^3', 10, ...
             'sf(1:nkn) are c3 for spline intervals', 10, ...
             'sf(nkn+1) is c2 for last interval', 10, ...
	     'sf(nkn+2) is c1 for last interval'];
    ds.sf = 'Spline vector for fpol^2/2-rzero^2*bzero^2/2';
    index_in_y.sf = ny+(1:nkn+2);
    iy = [iy ny+(1:nkn+2)];
    ny = ny+nkn+2;
  elseif strcmp(str,'er')
    ds.er = ['Remaining fraction of flux error. ', 10, ...
      'The flux error is found by gs_eq_analysis and ', ...
      'can be subtracted *after invoking gs_response* with:', 10, ...
      '  psizr(:) = psizr(:) - dpsizrdx(:,nx)*er'];
    index_in_y.er = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'fl')
    ds.fl = 'Flux loop signals as defined by config fields mlc, mlv, mpl';
    if isfield(config,'mlc') & isfield(config,'mlv') & isfield(config,'mpl')
      mlc = config.mlc;
      mlv = config.mlv;
      mpl = config.mpl;
      nfl = size(mlc,1);
      index_in_y.fl = ny+(1:nfl);
      iy = [iy ny+(1:nfl)];
      ny = ny+nfl;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include flux loops (fl) among outputs! ', ...
       'At least one of fields mlc, mlv, mpl missing from config.'])
    end
  elseif strcmp(str,'lv')
    ds.lv = 'Loop voltage signals as defined by config fields mhc, mhv, mph';
    if isfield(config,'mhc') & isfield(config,'mhv') & isfield(config,'mph')
      mhc = config.mhc;
      mhv = config.mhv;
      mph = config.mph;
      nlv = size(mhc,1);
      index_in_y.lv = ny+(1:nlv);
      iy = [iy ny+(1:nlv)];
      ny = ny+nlv;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include flux loops for loop voltage (lv) among outputs! ', ...
       'At least one of fields mhc, mhv, mph missing from config.'])
    end
  elseif strcmp(str,'bp')
    ds.bp = 'Magnetic probe signals as defined by config fields gbc, gbv, gpb';
    if isfield(config,'gbc') & isfield(config,'gbv') & isfield(config,'gpb')
      gbc = config.gbc;
      gbv = config.gbv;
      gpb = config.gpb;
      nbp = size(gbc,1);
      index_in_y.bp = ny+(1:nbp);
      iy = [iy ny+(1:nbp)];
      ny = ny+nbp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include magnetic probes (bp) in outputs! ', ...
       'At least one of fields gbc, gbv, gpb missing from config.'])
    end
  elseif strcmp(str,'rog')
    ds.rog = 'Rogowski signals as defined by config field rldata';
    if isfield(config,'rldata')
      if size(config.rldata,1) == nc+nv+1 % This is what it should be
        rldata = config.rldata'; % The transpose is used in gs codes
      elseif size(config.rldata,2) == nc+nv+1 % Transposed? We'll take it!
        rldata = config.rldata;
      else
        rldata = zeros(0,nc+nv+1);
	disp(['Warning gs_configure: ', ....
	 'Failed to include Rogowski loops (rog) in outputs! ', ...
	 'The field rldata in config must have dimensions [nc+nv+1,nrl].'])
      end
      nrog = size(rldata,1);
      index_in_y.rog = ny+(1:nrog);
      iy = [iy ny+(1:nrog)];
      ny = ny+nrog;
    else
      disp(['Warning gs_configure: ', ....
       'Failed to include Rogowski loops (rog) among outputs! ', ...
       'The field rldata is missing from config.'])
    end
  elseif strcmp(str,'psidp')
    ds.psidp = 'Flux at diagnostic points rdp, zdp';
    if ndp > 0
      index_in_y.psidp = ny+(1:ndp);
      iy = [iy ny+(1:ndp)];
      ny = ny+ndp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include flux at diagnostic points (psidp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing from config.'])
    end
  elseif strcmp(str,'brdp')
    ds.brdp = 'Radial magnetic field at diagnostic points rdp, zdp';
    if ndp > 0
      index_in_y.brdp = ny+(1:ndp);
      iy = [iy ny+(1:ndp)];
      ny = ny+ndp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include radial magnetic field at diagnostic points (brdp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing from config.'])
    end
  elseif strcmp(str,'bzdp')
    ds.bzdp = 'Vertical magnetic field at diagnostic points rdp, zdp';
    if ndp > 0
      index_in_y.bzdp = ny+(1:ndp);
      iy = [iy ny+(1:ndp)];
      ny = ny+ndp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include vertical magnetic field at diagnostic points (bzdp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing from config.'])
    end
  elseif strcmp(str,'btdp')
    ds.btdp = 'Toroidal magnetic field at diagnostic points rdp, zdp';
    if ndp > 0
      index_in_y.btdp = ny+(1:ndp);
      iy = [iy ny+(1:ndp)];
      ny = ny+ndp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include vertical magnetic field at diagnostic points (btdp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing from config.'])
    end
  elseif strcmp(str,'gapdp')
    ds.gapdp = 'Distance from rdp, zdp to boundary in direction toward axis';
    if ndp > 0
      index_in_y.gapdp = ny+(1:ndp);
      iy = [iy ny+(1:ndp)];
      ny = ny+ndp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include gaps from rdp, zdp to boundary (gapdp) among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing from config.'])
    end
  elseif strcmp(str,'gapdpdot')
    ds.gapdpdot = 'Boundary speed toward axis along gapdp vector';
    if ndp > 0
      index_in_y.gapdpdot = ny+(1:ndp);
      iy = [iy ny+(1:ndp)];
      ny = ny+ndp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include gapdpdot among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing from config.'])
    end
  elseif strcmp(str,'gapdpdotA')
    ds.gapdpdotA = 'Boundary speed toward axis when vps, vcd = 0';
    if ndp > 0
      index_in_y.gapdpdotA = ny+(1:ndp);
      iy = [iy ny+(1:ndp)];
      ny = ny+ndp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include gapdpdotA among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing from config.'])
    end
  elseif strcmp(str,'gapdpdotB')
    ds.gapdpdotB = 'Boundary speed toward axis caused by vps, vcd';
    if ndp > 0
      index_in_y.gapdpdotB = ny+(1:ndp);
      iy = [iy ny+(1:ndp)];
      ny = ny+ndp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include gapdpdotB among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing from config.'])
    end
  elseif strcmp(str,'gapdpdotB2')
    ds.gapdpdotB2 = 'Boundary speed toward axis caused by vcd';
    if ndp > 0
      index_in_y.gapdpdotB2 = ny+(1:ndp);
      iy = [iy ny+(1:ndp)];
      ny = ny+ndp;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include gapdpdotB2 among outputs! ', ...
       10, '   At least one of fields rdp, zdp missing from config.'])
    end
  elseif strcmp(str,'dgapdpdvps')
    ds.dgapdpdvps = 'd(gapdp)/dt = reshape(dgapdpdvps,ndp,nps)*vps';
    if ndp > 0 & nps > 0
      index_in_y.dgapdpdvps = ...
        ny+(1:ndp)'*ones(1,nps)+ones(ndp,1)*(0:nps-1)*ndp;
      %iy = [iy ny+(1:ndp*nps)];
      ny = ny+ndp*nps;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include dgapdpdvps among outputs! '])
      if ndp == 0
         disp('   At least one of fields rdp, zdp missing from config.')
      end
      if nps == 0
         disp('   Number of power supplies (nps) was set to zero.')
      end
    end
  elseif strcmp(str,'dgapdpdvcd')
    ds.dgapdpdvcd = 'd(gapdp)/dt = reshape(dgapdpdvcd,ndp,nps)*vcd';
    if ndp > 0 & ncd > 0
      index_in_y.dgapdpdvcd = ...
        ny+(1:ndp)'*ones(1,ncd)+ones(ndp,1)*(0:ncd-1)*ndp;
      %iy = [iy ny+(1:ndp*ncd)];
      ny = ny+ndp*ncd;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include dgapdpdvcd among outputs! '])
      if ndp == 0
         disp('   At least one of fields rdp, zdp missing from config.')
      end
      if ncd == 0
         disp('   Number of H/CD actuators (ncd) was set to zero.')
      end
    end
  elseif strcmp(str,'dgapdpdxs')
    ds.dgapdpdxs = 'd(gapdp)/dxs';
    if ndp > 0
      index_in_y.dgapdpdxs = ...
        ny+(1:ndp)'*ones(1,nxs)+ones(ndp,1)*(0:nxs-1)*ndp;
      ny = ny+ndp*nxs;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include dgapdpdxs among outputs! '])
      if ndp == 0
         disp('   At least one of fields rdp, zdp missing from config.')
      end
    end
  elseif strcmp(str,'gaps')
    ds.gaps = 'Gaps from r,z along gr,gz to separatrix, gapspec=[r,z,gr,gz]';
    if ngaps > 0
      index_in_y.gaps = ny+(1:ngaps);
      iy = [iy ny+(1:ngaps)];
      ny = ny+ngaps;
    else
      disp(['Warning gs_configure: ', ...
       'Failed to include gaps among outputs! ', ...
       'No gaps were specified in config.gapspec'])
    end
  elseif strcmp(str,'ys')
    ds.ys = 'Flux at all conductors';
    index_in_y.ys = ny+(1:nc+nv);
    iy = [iy ny+(1:nc+nv)];
    ny = ny+nc+nv;
  elseif strcmp(str,'rx')
    if     xpointorder == 1
      ds.rx = ['Radius for ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to normalized flux'];
    elseif xpointorder == 2
      ds.rx = ['Radius for ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to height'];
    elseif xpointorder == 3
      ds.rx = ['Radius for ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to normalized flux, then tracked'];
    elseif xpointorder == 4
      ds.rx = ['Radius for ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to height, then tracked'];
    end
    index_in_y.rx = ny+(1:nxpoints);
    iy = [iy ny+(1:nxpoints)];
    ny = ny+nxpoints;
  elseif strcmp(str,'zx')
    if     xpointorder == 1
      ds.rx = ['Height for ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to normalized flux'];
    elseif xpointorder == 2
      ds.rx = ['Height for ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to height'];
    elseif xpointorder == 3
      ds.rx = ['Height for ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to normalized flux, then tracked'];
    elseif xpointorder == 4
      ds.rx = ['Height for ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to height, then tracked'];
    end
    index_in_y.zx = ny+(1:nxpoints);
    iy = [iy ny+(1:nxpoints)];
    ny = ny+nxpoints;
  elseif strcmp(str,'psix')
    if     xpointorder == 1
      ds.rx = ['The flux at ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to normalized flux'];
    elseif xpointorder == 2
      ds.rx = ['The flux at ', num2str(nxpoints), ...
        ' x-points inside limiter, sorted according to height'];
    elseif xpointorder == 3
      ds.rx = ['The flux at ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to normalized flux, then tracked'];
    elseif xpointorder == 4
      ds.rx = ['The flux at ', num2str(nxpoints), ...
        ' x-points inside limiter, first sorted according to height, then tracked'];
    end
    index_in_y.psix = ny+(1:nxpoints);
    iy = [iy ny+(1:nxpoints)];
    ny = ny+nxpoints;
  elseif strcmp(str,'rcur')
    ds.rcur = 'Radius of current centroid';
    index_in_y.rcur = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'zcur')
    ds.zcur = 'Height of current centroid';
    index_in_y.zcur = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'zcurdot')
    ds.zcurdot = 'Vertical speed of current centroid';
    index_in_y.zcurdot = ny+1;
    ny = ny+1;
  elseif strcmp(str,'zcurdotA')
    ds.zcurdotA = 'Vertical speed of current centroid driven passively';
    index_in_y.zcurdotA = ny+1;
    ny = ny+1;
  elseif strcmp(str,'zcurdotB')
    ds.zcurdotB = 'Vertical speed of current centroid driven by actuators';
    index_in_y.zcurdotB = ny+1;
    ny = ny+1;
  elseif strcmp(str,'zcurdotB2')
    ds.zcurdotB2 = 'Vertical speed of current centroid driven by H&CD';
    index_in_y.zcurdotB2 = ny+1;
    ny = ny+1;
  elseif strcmp(str,'dzcurdotdvps')
    ds.dzcurdotdvps = 'Vertical speed response to circuit voltage from power supplies';
    index_in_y.dzcurdotdvps = ny+(1:nps);
    ny = ny+nps;
  elseif strcmp(str,'dzcurdotdvcd')
    ds.dzcurdotdvcd = 'Vertical speed response to heating & current-drive actuators';
    index_in_y.dzcurdotdvcd = ny+(1:ncd);
    ny = ny+ncd;
  elseif strcmp(str,'cpasma')
    ds.cpasma = 'Total toroidal plasma current';
    index_in_y.cpasma = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'cpasmadot')
    ds.cpasma = 'Time derivative of plasma current';
    index_in_y.cpasmadot = ny+1;
    ny = ny+1;
  elseif strcmp(str,'cpasmadotA')
    ds.cpasmadotA = 'Time derivative of plasma current without actuators';
    index_in_y.cpasmadotA = ny+1;
    ny = ny+1;
  elseif strcmp(str,'cpasmadotB')
    ds.cpasmadotB = 'Time derivative of plasma current by actuators';
    index_in_y.cpasmadotB = ny+1;
    ny = ny+1;
  elseif strcmp(str,'cpasmadotB2')
    ds.cpasmadotB2 = 'Time derivative of plasma current by H&CD';
    index_in_y.cpasmadotB2 = ny+1;
    ny = ny+1;
  elseif strcmp(str,'dcpasmadotdvps')
    ds.cpasma = 'cpasmadot response to circuit voltage from power supplies';
    index_in_y.dcpasmadotdvps = ny+(1:nps);
    ny = ny+nps;
  elseif strcmp(str,'dcpasmadotdvcd')
    ds.cpasma = 'cpasmadot response to heating & current-drive actuators';
    index_in_y.dcpasmadotdvcd = ny+(1:ncd);
    ny = ny+ncd;
  elseif strcmp(str,'aminor')
    ds.aminor = 'Horizontal minor radius of plasma';
    index_in_y.aminor = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'nbbbs')
    ds.nbbbs = 'number of valid points in rbbbs and zbbbs';
    index_in_y.nbbbs = ny+1;
    ny = ny+1;
  elseif strcmp(str,'rbbbs')
    ds.rbbbs = 'Radius of boundary points';
    index_in_y.rbbbs = ny+(1:nbbbs_max);
    ny = ny+nbbbs_max;
  elseif strcmp(str,'zbbbs')
    ds.zbbbs = 'Height of boundary points';
    index_in_y.zbbbs = ny+(1:nbbbs_max);
    ny = ny+nbbbs_max;
  elseif strcmp(str,'rmaxis')
    ds.rmaxis = 'Radius of magnetic axis';
    index_in_y.rmaxis = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'zmaxis')
    ds.zmaxis = 'Height of magnetic axis';
    index_in_y.zmaxis = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'psimag')
    ds.psimag = 'Flux at magnetic axis';
    index_in_y.psimag = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'rbdef')
    ds.rbdef = 'Radius of point that defines the boundary';
    index_in_y.rbdef = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'zbdef')
    ds.zbdef = 'Height of point that defines the boundary';
    index_in_y.zbdef = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'psibry')
    ds.psibry = 'Flux at boundary';
    index_in_y.psibry = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'li')
    ds.li = 'Normalized inductance';
    index_in_y.li = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'betap')
    ds.betap = 'Poloidal beta';
    index_in_y.betap = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'betan')
    ds.betan = 'Normalized toroidal beta';
    index_in_y.betan = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'Wth')
    ds.Wth = 'Thermal energy';
    index_in_y.Wth = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'psipla')
    ds.psipla = 'Plasma flux';
    index_in_y.psipla = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'bp2flx')
    ds.bp2flx = '(mu0*cpasma/Cl)^2';
    index_in_y.bp2flx = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'psiplaapp')
    ds.psiplaapp = 'Plasma flux from conductors';
    index_in_y.psiplaapp = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'vplas')
    ds.vplas = 'Plasma voltage from conductors';
    index_in_y.vplas = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'Lpla')
    ds.Lpla = 'Plasma self inductance';
    index_in_y.Lpla = ny+1;
    %iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'lconn')
    ds.lconn = 'Representative connection length';
    index_in_y.lconn = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'Ltot')
    ds.Ltot = 'Area-integral of 1/R';
    index_in_y.Ltot = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'Atot')
    ds.Atot = 'Area of plasma cross section';
    index_in_y.Atot = ny+1;
    ny = ny+1;
  elseif strcmp(str,'Vtot')
    ds.Vtot = 'Total plasma volume';
    index_in_y.Vtot = ny+1;
    iy = [iy ny+1];
    ny = ny+1;
  elseif strcmp(str,'rhot')
    ds.rhot = 'Square root of normalized toroidal flux';
    index_in_y.rhot = ny+(1:nr);
    iy = [iy ny+(1:nr)];
    ny = ny+nr;
    calculate_profiles = 1;
    rhot = zeros(nr,1);
    drhotdx = zeros(nr,nx);
  elseif strcmp(str,'jtav')
    ds.jtav = 'Flux-surface averaged toroidal current density profile';
    index_in_y.jtav = ny+(1:nr);
    iy = [iy ny+(1:nr)];
    ny = ny+nr;
    calculate_profiles = 1;
    jtav = zeros(nr,1);
    djtavdx = zeros(nr,nx);
  elseif strcmp(str,'vres')
    ds.vres = 'Resistive voltage profile';
    index_in_y.vres = ny+(1:ncont);
    iy = [iy ny+(1:ncont)];
    ny = ny+ncont;
    calculate_helical_voltage = 1;
    vres = zeros(ncont,1);
    dvresdx = zeros(ncont,nx);
  elseif strcmp(str,'vind')
    ds.vind = 'Inductive voltage profile';
    index_in_y.vind = ny+(1:ncont);
    iy = [iy ny+(1:ncont)];
    ny = ny+ncont;
    calculate_helical_voltage = 1;
    vind = zeros(ncont,1);
    dvindcdxdot = zeros(ncont,nx);
  elseif strcmp(str,'jpar')
    ds.jpar = 'Profile of parallel current';
    index_in_y.jpar = ny+(1:ncont);
    iy = [iy ny+(1:ncont)];
    ny = ny+ncont;
    calculate_helical_voltage = 1;
    jpar = zeros(ncont,1);
    djpardx = zeros(ncont,nx);
  elseif strcmp(str,'T')
    ds.T = 'Profile of toroidal flux';
    index_in_y.T = ny+(1:ncont);
    iy = [iy ny+(1:ncont)];
    ny = ny+ncont;
    calculate_helical_voltage = 1;
    Tc = zeros(ncont,1);
    dTcdx = zeros(ncont,nx);
  elseif strcmp(str,'qpsi')
    ds.qpsi = 'Safety factor q';
    index_in_y.qpsi = ny+(1:ncont);
    iy = [iy ny+(1:ncont)];
    ny = ny+ncont;
    calculate_helical_voltage = 1;
    qc = zeros(ncont,1);
    dqcdx = zeros(ncont,nx);
  elseif strcmp(str,'fluxerror')
    ds.fluxerror = 'Error in normalized flux, < 1e-2 for valid Grad-Shafranov solutions';
    index_in_y.fluxerror = ny+1;
    ny = ny+1;
  elseif strcmp(str,'time')
    ds.time = 'Simulation time (used for plots and time-dependent configuration data)';
    index_in_y.time = ny+1;
    ny = ny+1;
  elseif strcmp(str,'gamma')
    ds.gamma = 'Growth rate of most unstable mode';
    index_in_y.gamma = ny+1;
    ny = ny+1;
  elseif strcmp(str,'drcurdv')
    ds.drcurdv = 'Response of rcur to most unstable eigen vector';
    index_in_y.drcurdv = ny+1;
    ny = ny+1;
  elseif strcmp(str,'dzcurdv')
    ds.dzcurdv = 'Response of zcur to most unstable eigen vector';
    index_in_y.dzcurdv = ny+1;
    ny = ny+1;
  elseif strcmp(str,'response_count')
    ds.response_count = 'Number of calls to gs_response or gscp_analysis_response';
    index_in_y.response_count = ny+1;
    ny = ny+1;
  elseif strcmp(str,'frc')
    ds.frc = 'Radial forces on coils';
    index_in_y.frc = ny+(1:nc);
    iy = [iy ny+(1:nc)];
    ny = ny+nc;
    need_mrmz = 1;
    frc = zeros(nc,1);
    dfrcdx = zeros(nc,nx);
  elseif strcmp(str,'fzc')
    ds.fzc = 'Vertical forces on coils';
    index_in_y.fzc = ny+(1:nc);
    iy = [iy ny+(1:nc)];
    ny = ny+nc;
    need_mrmz = 1;
    fzc = zeros(nc,1);
    dfzcdx = zeros(nc,nx);
  elseif ~isempty(str2num(outputs(i,:))) % This is just a number
    n = str2num(outputs(i,:));
    zcount = zcount + 1;
    ds = setfield(ds, ['zeros', num2str(zcount)], [num2str(n), ' zero-signals']);
    index_in_y = setfield(index_in_y, ['zeros', num2str(zcount)], ny+(1:n));
    ny = ny+n;
  else % Check if this might be a user-defined signal
    n = 0;
    if isfield(config,[str '0'])
      n = max(n,size(getfield(config,[str '0']),1));
    end
    if isfield(config,['d' str 'dic'])
      n = max(n,size(getfield(config,['d' str 'dic']),1));
    end
    if isfield(config,['d' str 'div'])
      n = max(n,size(getfield(config,['d' str 'div']),1));
    end
    if isfield(config,['d' str 'dpcurrt'])
      n = max(n,size(getfield(config,['d' str 'dpcurrt']),1));
    end
    if n > 0
      usig = usig + 1;
      user_signal(usig).name = str;
      ds = setfield(ds, str, 'user-defined signal');
      user_signal(usig).y0 = zeros(n,1);
      if isfield(config,[str '0'])
	signal0 = getfield(config,[str '0']);
	for j = 1:min(n,length(signal0))
	  user_signal(usig).y0(j) = signal0(j);
	end
      end
      user_signal(usig).dydic = zeros(n,nc);
      if isfield(config,['d' str 'dic'])
	dsignaldic = getfield(config,['d' str 'dic']);
	for j = 1:min(n,size(dsignaldic,1))
	for k = 1:min(nc,size(dsignaldic,2))
	  user_signal(usig).dydic(j,k) = dsignaldic(j,k);
	end
	end
      end
      user_signal(usig).dydiv = zeros(n,nv);
      if isfield(config,['d' str 'div'])
	dsignaldiv = getfield(config,['d' str 'div']);
	for j = 1:min(n,size(dsignaldiv,1))
	for k = 1:min(nv,size(dsignaldiv,2))
	  user_signal(usig).dydiv(j,k) = dsignaldiv(j,k);
	end
	end
      end
      user_signal(usig).dydpcurrt = zeros(n,ngg);
      if isfield(config,['d' str 'dpcurrt'])
	dsignaldpcurrt = getfield(config,['d' str 'dpcurrt']);
	for j = 1:min(n,size(dsignaldpcurrt,1))
	for k = 1:min(ngg,size(dsignaldpcurrt,2))
	  user_signal(usig).dydpcurrt(j,k) = dsignaldpcurrt(j,k);
	end
	end
      end
      index_in_y = setfield(index_in_y, str, ny+(1:n));
      ny = ny+n;
    else
      disp(['Warning in gs_configure: output ' str ' is not recognized'])
    end
  end
end
if ny == 0
  ds.psizr = 'poloidal flux on the grid';
  index_in_y.psizr = ny+(1:ngg);
  iy = [iy ny+(1:nr)];
  ny = ny+ngg;
end
index_in_y.descriptions = ds;
lae.y = zeros(ny,1);
dydx = zeros(ny,nx);
if isfield(index_in_y,'ic')
  dydx(index_in_y.ic,indic) = eye(nc);
end
if isfield(index_in_y,'iv')
  dydx(index_in_y.iv,indiv) = eye(nv);
end
if isfield(index_in_y,'sp')
  dydx(index_in_y.sp,indsp) = eye(nkn+2);
end
if isfield(index_in_y,'sf')
  dydx(index_in_y.sf,indsf) = eye(nkn+2);
end
if isfield(index_in_y,'er')
  dydx(index_in_y.er,inder) = eye(1);
end
dgapdpdx = zeros(ndp,nx);


% For making y continuous when plasma response is updated
if isfield(config,'phaseouterrors')
  phaseouterrors  = config.phaseouterrors;
end
if ~exist('phaseouterrors','var') | isempty(phaseouterrors)
  phaseouterrors = 1;
end
err.y = zeros(ny,1);
err.ys = zeros(nc+nv,1);


% Circular plasma model
nrh = 2*nr-1;
rh = linspace(1,2,nrh); % This is just memory allocation
drh = 1/(nrh-1); % This is just memory allocation
gscp_configure

% For parsing efits with cc_efit_to_tok
if isfield(config,'imks')
  imks = config.imks;
end
if ~exist('imks','var') | isempty(imks)
  imks = 1;
end
if isfield(config,'iterminal')
  iterminal = config.iterminal;
end
if ~exist('iterminal','var') | isempty(iterminal)
  iterminal = 1;
end
if isfield(config,'ecnturn')
  ecnturn = config.ecnturn;
end
if ~exist('ecnturn','var')
  ecnturn = [];
end
if isfield(config,'fcnturn')
  fcnturn = config.fcnturn;
end
if ~exist('fcnturn','var')
  fcnturn = [];
end
if isfield(config,'idx_efit_to_tok')
  idx_efit_to_tok = config.idx_efit_to_tok;
end
if ~exist('idx_efit_to_tok','var')
  idx_efit_to_tok = [];
end


% Check for common issues
if constraints
  p0 = c0(iknotg,:)*sp0;
  p1 = c1(iknotg,:)*sp0;
  p2 = c2(iknotg,:)*sp0;
  p3 = c3(iknotg,:)*sp0;
  pres0 = p0 + p1.*psibar + p2.*psibar.^2 + p3.*psibar.^3;
  if any(pres0 < 0)
    if isfield(config,'pres0')
    elseif isfield(config,'pprime0')
    end
  end
end



if need_mrmz
  if isfield(config,'mrcc') & isfield(config,'mzcc') & ...
     isfield(config,'mrpc') & isfield(config,'mzpc')
    mrcc = config.mrcc;
    mzcc = config.mzcc;
    mrpc = config.mrpc;
    mzpc = config.mzpc;
  elseif isfield(config,'fcdata') & isfield(config,'fcnturn')
    fcdata = config.fcdata;
    fcnturn = config.fcnturn;
    mrcc = zeros(nc);
    mzcc = zeros(nc);
    mrpc = zeros(ngg,nc);
    mzpc = zeros(ngg,nc);
    nfc = size(fcdata,2);
    nec = nc-nfc;
    for i = 1:nfc
      fprintf(...
'Creating objects for coil force calculations: coil %d of %d...',i+nec,nc)
      ri = fcdata(2,i)+[-1 -1 1 1 -1]/2*fcdata(4,i);
      zi = fcdata(1,i)+[-1 1 1 -1 -1]/2*fcdata(3,i);
      for j = 1:nfc
	rj = fcdata(2,j)+[-1 -1 1 1 -1]/2*fcdata(4,j);
	zj = fcdata(1,j)+[-1 1 1 -1 -1]/2*fcdata(3,j);
	[ms, mrs, mzs] = mpolygon2polygon(ri,zi,rj,zj);
	%MCC(i,j) = fcnturn(i)*fcnturn(j)*ms;
	mrcc(i+nec,j+nec) = fcnturn(i)*fcnturn(j)*mrs;
	mzcc(i+nec,j+nec) = fcnturn(i)*fcnturn(j)*mzs;
      end
      % For validation
      %[ms2, mrs2, mzs2] = mpolygon2point(ri,zi,rgg(:),zgg(:),2);
      %MPC(:,i) = fcnturn(i)*ms2;
      %GBR2C(:,i) = -fcnturn(i)*mzs2./rgg(:)/2/pi;
      %GBZ2C(:,i) = +fcnturn(i)*mrs2./rgg(:)/2/pi;
      [ms, mrs, mzs] = mpolygon2point(ri,zi,rgg(:),zgg(:),1);
      mrpc(:,i+nec) = fcnturn(i)*mzs;
      mzpc(:,i+nec) = fcnturn(i)*mrs;
      fprintf(char(8+zeros(1,80)))
    end
    fprintf('                                                               ')
    fprintf(char(8+zeros(1,80)))
  else
    disp('Warning gs_configure: mrcc, mzcc, mrpc, mzpc can not be calculated')
    disp('The fields fcdata and fcnturn are missing from config')
  end
end
