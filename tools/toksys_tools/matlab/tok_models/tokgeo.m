function istok = tokgeo(tok,varname)
%
%  USAGE:  tokgeo          % Find and show tokamaks in the workspace
%          tokgeo(tok)     % Show the tokamak in tok
%          tokgeo(tok,opt) % Show tok with menu options in opt
%
%  INPUTS: tok, tokamak geometry data in TokSys format
%          opt, initial menu, type "tokgeo options" for a description
%
%  OUTPUTS: Figure with extra menus "View", "Lines', "Surfaces", "Names"

%  WRITTEN BY: Anders Welander ON 2018-06-27
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If nargin==0 then make recursive calls with nargin=2
if nargin == 0
  % Need to find a tok, check caller's workspace
  a = evalin('caller','whos');
  n = 0; % Will contain number of tokamaks that were found in caller's workspace
  for i = 1:numel(a)
    if strcmp(a(i).class,'struct')
      tok = evalin('caller',a(i).name);
      if tokgeo(tok,a(i).name) % true if tok does contain a tokamak
	n = n+1;
        if n == 1
	  %disp('Variable containing tokamak geometry:') 
	end
        %disp(a(i).name)
      end
    end
  end
  if n == 0
    disp('No TokSys tokamak variable was found in the workspace')
  end
  return
end


% Color selection
colorSelection = {'black','white','red','orange','yellow','green','cyan','blue','magenta',...
  'bronze','silver','gold','cocoa','peach','wheat','almond','pink','beige','azure','plum'};

ncolors = numel(colorSelection);

black     = [000 000 000]/255;
white     = [255 255 255]/255;
red       = [255 000 000]/255;
orange    = [255 128 000]/255;
yellow    = [255 255 000]/255;
green     = [000 255 000]/255;
cyan      = [000 255 255]/255;
blue      = [000 000 255]/255;
magenta   = [255 000 255]/255;
bronze    = [205 127 050]/255; % Original bronze color
bronze    = [136 102 054]/255; % Consistent with image
silver    = [192 192 192]/255; % Consistent with image
gold      = [207 181 059]/255; % Original gold color
gold      = [219 170 046]/255; % Consistent with image
cocoa     = [094 067 048]/255;
plum      = [221 160 221]/255;
raspberry = [233 000 103]/255;
hotpink   = [255 105 180]/255;
pink      = [255 192 203]/255;
beige     = [245 245 220]/255;
almond    = [238 217 196]/255;
wheat     = [245 222 179]/255;
peach     = [255 218 185]/255;
linen     = [250 240 230]/255;
ivory     = [255 255 240]/255;
azure     = [240 255 255]/255;

% Background images
imstat = zeros(1,ncolors); % 0 = didn't check for image, 1 = no image, 2 = image exists
im{ncolors} = [];          % Loaded images


% Allow the user to make the call with tokgeo tok instead of tokgeo(tok)
if ischar(tok)
  try
    tok = evalin('caller',tok);
  catch
    disp(' ')
    disp('opt can be used to initialize the extra menus')
    disp('__________________________________________________________________________________')
    disp(' ')
    disp('opt.View.is3D (boolean) sets the view to 3D (default=false)')
    disp('opt.View.Fly (double) fly into machine and back out opt.View.Fly times')
    disp('opt.View.Spin (double) gives an initial spin to the view opt.View.Spin turns')
    disp('opt.View.Grid (boolean) execute grid on')
    disp('opt.View.GridUnits (boolean) switch to grid units tok.rg, tok.zg (2D view only)')
    disp('opt.View.whitebg (double 1x3) rgb triplet for background color (default=[1 1 1])')
    disp('opt.View.whitebg can alternatively be one of the unabbreviated strings:')
    str = '';
    for i = 1:ncolors
      str = [str ' ' char(colorSelection(i))];
      if round(i/12) == i/12 & i ~= ncolors
        str(end+1) = 10;
      end
    end
    disp(str)
    disp('Lines and surfaces change to 1-normal color value for rgb channels that are zero')
    disp('opt.View.Flytime (double>0) duration of flight in and out of machine (default=10)')
    disp('opt.View.Spinrpm (double~=0) rotations per minute during spin (default=20)')
    disp('__________________________________________________________________________________')
    disp(' ')
    disp('opt.Lines.LIM (boolean) show contour of limiter cross section (default=true)')
    disp('opt.Lines.EC (boolean) show contour of E coil cross section (default=true)')
    disp('opt.Lines.FC (boolean) show contour of F coil cross sections (default=true)')
    disp('opt.Lines.VV (boolean) show contour of Vessel cross sections (default=true)')
    disp('opt.Lines.FL (boolean) show contour of Flux loop cross sections (default=false)')
    disp('opt.Lines.LV (boolean) show contour of Voltage loop cross sections (default=false)')
    disp('opt.Lines.BP (boolean) show B probes (default=false)')
    disp('opt.Lines.MSE (boolean) show MSE points (default=false)')
    disp('opt.Lines.TD (boolean) show 3D coils (default=false)')
    disp('opt.Lines.SL (boolean) show Saddle loops (default=false)')
    disp('opt.Lines.Width (double>0) line width for contour objects (default=1)')
    disp('__________________________________________________________________________________')
    disp(' ')
    disp('opt.Surfaces.LIM (boolean) show limiter surface in 3D view (default=true)')
    disp('opt.Surfaces.EC (boolean) show E coil surface in 3D view (default=true)')
    disp('opt.Surfaces.FC (boolean) show F coil surfaces in 3D view (default=true)')
    disp('opt.Surfaces.VV (boolean) show Vessel surfaces in 3D view (default=true)')
    disp('opt.Surfaces.FL (boolean) show Flux loop surfaces in 3D view (default=false)')
    disp('opt.Surfaces.LV (boolean) show Voltage loop surfaces in 3D view (default=false)')
    disp('opt.Surfaces.EdgeColor (double) 0=black, 5=color, 10=white (default=0)')
    disp('opt.Surfaces.FaceColor (double) 0=black, 5=color, 10=white (default=5)')
    disp('opt.Surfaces.Dphi (double) toroidal extent, range -360 to +360 (default=180)')
    disp('opt.Surfaces.dphi (double) angle between lines, range >0 to 45 (default=10)')
    disp('__________________________________________________________________________________')
    disp(' ')
    disp('opt.Names.EC (boolean) show names of E coils (default=false)')
    disp('opt.Names.FC (boolean) show names of F coils (default=false)')
    disp('opt.Names.VV (boolean) show names of Vessel elements (default=false)')
    disp('opt.Names.FL (boolean) show names of flux loops (default=false)')
    disp('opt.Names.LV (boolean) show names of voltage loops (default=false)')
    disp('opt.Names.BP (boolean) show names of B probes (default=false)')
    disp('opt.Names.MSE (boolean) show names of MSE points (default=false)')
    disp('opt.Names.TD (boolean) show names of 3D coils (default=false)')
    disp('opt.Names.SL (boolean) show names of saddle loops (default=false)')
    disp('opt.Names.FontSize (double>0) font size for names (default=12)')
    disp('__________________________________________________________________________________')
    disp(' ')
    disp('Initial 3D view is view(opt.az,opt.el) if both fields exist (default view(-30,45))')
    disp('Initial axis is axis(opt.axis) if the field axis exists')
    disp('Initial CameraViewAngle can be specified with opt.CameraViewAngle')
    disp('Initial CameraPosition can be specified with opt.CameraPosition')
    disp('Initial CameraTarget can be specified with opt.CameraTarget')
    disp('opt.imbg can be image or filename for image to use as background (e.g. goldbg.jpg)')
    disp('opt.imbgf is like opt.imbg but the background is not removed by changing whitebg')
    disp('__________________________________________________________________________________')
    disp(' ')
    disp('Fun tips for more 3D effects:')
    disp('Check out the camera toolbar under the regular View menu')
    disp('Create a hole in the vessel by panning the tokamak to partly outside the 3dbox')
    disp(' ')
    return
  end
end

% Decide if tok is a tokamak, a very low bar is used to allow tokamaks under development
if nargout > 0
  istok = isfield(tok,'fcdata');
  if ~istok
    return
  end
end

% If tokgeo is called without istok as output then always show whatever can be shown
figure

% Second input is either a variable name or opt
opt = [];
if nargin > 1
  if isstruct(varname)
    opt = varname;
  else
    set(gcf,'Name',varname)
  end
end
if ~isfield(opt,'View')
  opt.View = [];
end
if ~isfield(opt,'Lines')
  opt.Lines = [];
end
if ~isfield(opt,'Surfaces')
  opt.Surfaces = [];
end
if ~isfield(opt,'Names')
  opt.Names = [];
end

if isfield(tok,'tokamak') & ischar(tok.tokamak)
  tokamak = tok.tokamak;
else
  tokamak = '';
end


% Many variables are declared here just to make them available in all functions within this function

% Handle to axis
h = nan;

% Handle to background image
hbg = [];

% Handles used in name searches
d = [];     % Handle to dialog box
term = [];  % Handle to edit window with search string
texth = []; % True if first input to text can be an axis handle

% Create all data needed to draw lines and display names
Xline = {};
Yline = {};
Zline = {};
Rline = {};
Cline = []; % Colors
Xname = [];
Yname = [];
Zname = [];
Rname = [];
Sname = {}; % Names
N = 0; % Number of elements

id.LIM = [];
if isfield(tok,'limdata')
  [n1,n2] = size(tok.limdata);
  if n1 > n2
    r = tok.limdata(:,2);
    z = tok.limdata(:,1);
  else
    r = tok.limdata(2,:);
    z = tok.limdata(1,:);
  end
  [~,m] = max(z);
  Xline(end+1) = {r};
  Yline(end+1) = {r*0};
  Zline(end+1) = {z};
  Rline(end+1) = {r};
  Cline(end+1,1:3) = [0 0 0];
  Xname(end+1) = r(m);
  Yname(end+1) = 0;
  Zname(end+1) = z(m);
  Rname(end+1) = r(m);
  Sname(end+1) = {'Limiter'};
  N = N+1;
  id.LIM(end+1) = N;
end

id.EC = [];
if isfield(tok,'ecdata')
  if size(tok.ecdata,1) < 6
    tok.ecdata(6,:) = 0;
  end
  nec = max(tok.ecdata(5,:));
  if nec == 0
    colors = ones(99,1)*[0 1 0];
  else
    colors = jet(3*nec);
  end
  for k = unique(tok.ecdata(5,:))
    r = [];
    z = [];
    for j = find(tok.ecdata(5,:) == k)
      r = [r tok.ecdata(2,j)+[-1 -1 +1 +1 -1]*tok.ecdata(4,j)/2 NaN];
      z = [z tok.ecdata(1,j)+[-1 +1 +1 -1 -1]*tok.ecdata(3,j)/2 NaN];
    end
    [~,m] = max(z);
    Xline(end+1) = {r};
    Yline(end+1) = {r*0};
    Zline(end+1) = {z};
    Rline(end+1) = {r};
    Cline(end+1,1:3) = colors(nec+k,:);
    Xname(end+1) = r(m);
    Yname(end+1) = 0;
    Zname(end+1) = z(m);
    Rname(end+1) = r(m);
    if isfield(tok,'ecnames') & k > 0 & size(tok.ecnames,1) >= k
      Sname(end+1) = {deblank(tok.ecnames(k,:))};
    else
      Sname(end+1) = {''};
    end
    N = N+1;
    id.EC(end+1) = N;
  end
end

id.FC = [];
if isfield(tok,'fcdata')
  if size(tok.fcdata,1) < 6
    tok.fcdata(6,:) = 0;
  end
  [r,z] = data2poly(tok.fcdata);
  for j = 1:size(r,2)
    Xline(end+1) = {r(:,j)};
    Yline(end+1) = {r(:,j)*0};
    Zline(end+1) = {z(:,j)};
    Rline(end+1) = {r(:,j)};
    Cline(end+1,1:3) = [1 0 0];
    Xname(end+1) = tok.fcdata(2,j);
    Yname(end+1) = 0;
    Zname(end+1) = tok.fcdata(1,j);
    Rname(end+1) = tok.fcdata(2,j);
    if isfield(tok,'fcnames') & size(tok.fcnames,1) >= j
      Sname(end+1) = {deblank(tok.fcnames(j,:))};
    else
      Sname(end+1) = {''};
    end
    N = N+1;
    id.FC(end+1) = N;
  end
end

id.VV = [];
if isfield(tok,'vvdata')
  if size(tok.vvdata,1) < 6
    tok.vvdata(6,:) = 0;
  end
  [r,z] = data2poly(tok.vvdata);
  for j = 1:size(r,2)
    Xline(end+1) = {r(:,j)};
    Yline(end+1) = {r(:,j)*0};
    Zline(end+1) = {z(:,j)};
    Rline(end+1) = {r(:,j)};
    Cline(end+1,1:3) = [0 0 1];
    Xname(end+1) = tok.vvdata(2,j);
    Yname(end+1) = 0;
    Zname(end+1) = tok.vvdata(1,j);
    Rname(end+1) = tok.vvdata(2,j);
    if isfield(tok,'vvnames') & size(tok.vvnames,1) >= j
      Sname(end+1) = {deblank(tok.vvnames(j,:))};
    else
      Sname(end+1) = {''};
    end
    N = N+1;
    id.VV(end+1) = N;
  end
end

id.FL = [];
if isfield(tok,'fldata')
  if size(tok.fldata,1) < 6
    tok.fldata(6,:) = 0;
  end
  tok.fldata(3,tok.fldata(3,:)==0) = 0.01;
  tok.fldata(4,tok.fldata(4,:)==0) = 0.01;
  [r,z] = data2poly(tok.fldata);
  for j = 1:size(r,2)
    Xline(end+1) = {r(:,j)};
    Yline(end+1) = {r(:,j)*0};
    Zline(end+1) = {z(:,j)};
    Rline(end+1) = {r(:,j)};
    Cline(end+1,1:3) = [0 0 0];
    Xname(end+1) = tok.fldata(2,j);
    Yname(end+1) = 0;
    Zname(end+1) = tok.fldata(1,j);
    Rname(end+1) = tok.fldata(2,j);
    if isfield(tok,'flnames') & size(tok.flnames,1) >= j
      Sname(end+1) = {deblank(tok.flnames(j,:))};
    else
      Sname(end+1) = {''};
    end
    N = N+1;
    id.FL(end+1) = N;
  end
end

id.LV = [];
if isfield(tok,'lvdata')
  if size(tok.lvdata,1) < 6
    tok.lvdata(6,:) = 0;
  end
  tok.lvdata(3,tok.lvdata(3,:)==0) = 0.01;
  tok.lvdata(4,tok.lvdata(4,:)==0) = 0.01;
  [r,z] = data2poly(tok.lvdata);
  for j = 1:size(r,2)
    Xline(end+1) = {r(:,j)};
    Yline(end+1) = {r(:,j)*0};
    Zline(end+1) = {z(:,j)};
    Rline(end+1) = {r(:,j)};
    Cline(end+1,1:3) = [1 1 1]*0.25;
    Xname(end+1) = tok.lvdata(2,j);
    Yname(end+1) = 0;
    Zname(end+1) = tok.lvdata(1,j);
    Rname(end+1) = tok.lvdata(2,j);
    if isfield(tok,'lvnames') & size(tok.lvnames,1) >= j
      Sname(end+1) = {deblank(tok.lvnames(j,:))};
    else
      Sname(end+1) = {''};
    end
    N = N+1;
    id.LV(end+1) = N;
  end
end

id.BP = [];
if isfield(tok,'bpdata')
  if size(tok.bpdata,1) < 6
    tok.bpdata(6,:) = 0;
  end
  for j = 1:size(tok.bpdata,2)
    if tok.bpdata(4,j) > 0
      dl = tok.bpdata(4,j);
    else
      dl = 0.01;
    end
    r = tok.bpdata(2,j) + dl/2*[-1 1]*cos(tok.bpdata(3,j)*pi/180);
    z = tok.bpdata(1,j) + dl/2*[-1 1]*sin(tok.bpdata(3,j)*pi/180);    
    Xline(end+1) = {r*cos(tok.bpdata(5,j)*pi/180)};
    Yline(end+1) = {r*sin(tok.bpdata(5,j)*pi/180)};
    Zline(end+1) = {z};
    Rline(end+1) = {r};
    Cline(end+1,1:3) = [0 1 1];
    Xname(end+1) = tok.bpdata(2,j)*cos(tok.bpdata(5,j)*pi/180);
    Yname(end+1) = tok.bpdata(2,j)*sin(tok.bpdata(5,j)*pi/180);
    Zname(end+1) = tok.bpdata(1,j);
    Rname(end+1) = tok.bpdata(2,j);
    if isfield(tok,'bpnames') & size(tok.bpnames,1) >= j
      Sname(end+1) = {deblank(tok.bpnames(j,:))};
    else
      Sname(end+1) = {''};
    end
    N = N+1;
    id.BP(end+1) = N;
  end
end

id.MSE = [];
if isfield(tok,'msedata')
  if size(tok.msedata,1) < 6
    tok.msedata(6,:) = 0;
  end
  for j = 1:size(tok.msedata,2)
    if tok.msedata(4,j) > 0
      dl = tok.msedata(4,j);
    else
      dl = 0.01;
    end
    r = tok.msedata(2,j) + dl/2*[-1 1]*cos(tok.msedata(3,j)*pi/180);
    z = tok.msedata(1,j) + dl/2*[-1 1]*sin(tok.msedata(3,j)*pi/180);    
    Xline(end+1) = {r*cos(tok.msedata(5,j)*pi/180)};
    Yline(end+1) = {r*sin(tok.msedata(5,j)*pi/180)};
    Zline(end+1) = {z};
    Rline(end+1) = {r};
    Cline(end+1,1:3) = [1 0 1];
    Xname(end+1) = tok.msedata(2,j)*cos(tok.msedata(5,j)*pi/180);
    Yname(end+1) = tok.msedata(2,j)*sin(tok.msedata(5,j)*pi/180);
    Zname(end+1) = tok.msedata(1,j);
    Rname(end+1) = tok.msedata(2,j);
    if isfield(tok,'msenames') & size(tok.msenames,1) >= j
      Sname(end+1) = {deblank(tok.msenames(j,:))};
    else
      Sname(end+1) = {''};
    end
    N = N+1;
    id.MSE(end+1) = N;
  end
end

id.TD = [];
if isfield(tok,'tddata')
  if isempty(tok.tddata)
    ntd = 0;
  elseif iscell(tok.tddata)
    ntd = numel(tok.tddata);
  else
    ntd = 1;
  end
  for j = 1:ntd
    if iscell(tok.tddata)
      fil = tok.tddata{j};
    else
      fil = tok.tddata;
    end
    n = size(fil,1);
    x = nan(3*n,1);
    y = nan(3*n,1);
    z = nan(3*n,1);
    x(1:3:3*n-2,1) = fil(:,1);
    y(1:3:3*n-2,1) = fil(:,2);
    z(1:3:3*n-2,1) = fil(:,3);
    x(2:3:3*n-1,1) = fil(:,4);
    y(2:3:3*n-1,1) = fil(:,5);
    z(2:3:3*n-1,1) = fil(:,6);
    Xline(end+1) = {x};
    Yline(end+1) = {y};
    Zline(end+1) = {z};
    Rline(end+1) = {sqrt(x.^2+y.^2)};
    Cline(end+1,1:3) = [1 0 1];
    k = ~isnan(x);
    xm = mean(x(k));
    ym = mean(y(k));
    ph = angle(xm+1i*ym);
    rm = sqrt(mean(x(k).^2+y(k).^2));
    zm = (max(z)+min(z))/2;
    Xname(end+1) = rm*cos(ph);
    Yname(end+1) = rm*sin(ph);
    Zname(end+1) = zm;
    Rname(end+1) = rm;
    if isfield(tok,'tdnames') & size(tok.tdnames,1) >= j
      Sname(end+1) = {deblank(tok.tdnames(j,:))};
    else
      Sname(end+1) = {''};
    end
    N = N+1;
    id.TD(end+1) = N;
  end
end

id.SL = [];
if isfield(tok,'sldata')
  if isempty(tok.sldata)
    nsl = 0;
  elseif iscell(tok.sldata)
    nsl = numel(tok.sldata);
  else
    nsl = 1;
  end
  for j = 1:nsl
    if iscell(tok.sldata)
      fil = tok.sldata{j};
    else
      fil = tok.sldata;
    end
    n = size(fil,1);
    x = nan(3*n,1);
    y = nan(3*n,1);
    z = nan(3*n,1);
    x(1:3:3*n-2,1) = fil(:,1);
    y(1:3:3*n-2,1) = fil(:,2);
    z(1:3:3*n-2,1) = fil(:,3);
    x(2:3:3*n-1,1) = fil(:,4);
    y(2:3:3*n-1,1) = fil(:,5);
    z(2:3:3*n-1,1) = fil(:,6);
    Xline(end+1) = {x};
    Yline(end+1) = {y};
    Zline(end+1) = {z};
    Rline(end+1) = {sqrt(x.^2+y.^2)};
    Cline(end+1,1:3) = [0 1 0];
    k = ~isnan(x);
    xm = mean(x(k));
    ym = mean(y(k));
    ph = angle(xm+1i*ym);
    rm = sqrt(mean(x(k).^2+y(k).^2));
    zm = (max(z)+min(z))/2;
    Xname(end+1) = rm*cos(ph);
    Yname(end+1) = rm*sin(ph);
    Zname(end+1) = zm;
    Rname(end+1) = rm;
    if isfield(tok,'slnames') & size(tok.slnames,1) >= j
      Sname(end+1) = {deblank(tok.slnames(j,:))};
    else
      Sname(end+1) = {''};
    end
    N = N+1;
    id.SL(end+1) = N;
  end
end

% Default axis for 2D
Rmin = +inf;
Rmax = -inf;
Zmin = +inf;
Zmax = -inf;
for i = [id.LIM id.EC id.FC id.VV]
  Rmin = min(Rmin,min(Rline{i}));
  Rmax = max(Rmax,max(Rline{i}));
  Zmin = min(Zmin,min(Zline{i}));
  Zmax = max(Zmax,max(Zline{i}));
end
DR = Rmax-Rmin;
DZ = Zmax-Zmin;
a2D = [round2(Rmin-DR/20) round2(Rmax+DR/20) round2(Zmin-DZ/20) round2(Zmax+DZ/20)];
for i = [id.FL id.LV id.BP id.MSE id.TD id.SL]
  Rmin = min(Rmin,min(Rline{i}));
  Rmax = max(Rmax,max(Rline{i}));
  Zmin = min(Zmin,min(Zline{i}));
  Zmax = max(Zmax,max(Zline{i}));
end
DR = Rmax-Rmin;
DZ = Zmax-Zmin;
a2Dx = [round2(Rmin-DR/20) round2(Rmax+DR/20) round2(Zmin-DZ/20) round2(Zmax+DZ/20)];

% Surface parameters
phi0 = 0;
Xsurf = Xline;
Ysurf = Yline;
Zsurf = Zline;
Csurf = Cline;
Esurf = Cline;

% Handles for lines
Hline = nan(1,N);

% Handles for surfaces
Hsurf = nan(1,N);

% Handles for names
Hname = nan(1,N);

% Handles for objects found in searches
sHline = Hline;
sHname = Hname;

% Visibility flags
Fline = false(1,N);
Fsurf = false(1,N);
Fname = false(1,N);

% Flag for view
is2D = true;

% Grid
isgrid = false;    % True if grid on
gridunits = false; % Can be on in 2D view only

% In grid units ir = (r-rgmin)/dr+1
if isfield(tok,'rg')
  rgmin = tok.rg(1);
  dr = (tok.rg(end)-tok.rg(1))/(length(tok.rg)-1);
else
  rgmin = 1;
  dr = 1;
end
if isfield(tok,'zg')
  zgmin = tok.zg(1);
  dz = (tok.zg(end)-tok.zg(1))/(length(tok.zg)-1);
else
  zgmin = 1;
  dz = 1;
end

a2Dgu = [(a2D(1)-rgmin)/dr+1 (a2D(2)-rgmin)/dr+1 (a2D(3)-zgmin)/dz+1 (a2D(4)-zgmin)/dz+1];

% View
mView = uimenu('Label','View');
mView2D = uimenu(mView,'Label','2D','Callback',@view2D);
mView3D = uimenu(mView,'Label','3D','Callback',@view3D);
mViewFly = uimenu(mView,'Label','Fly','Callback',@Fly);
mViewSpin = uimenu(mView,'Label','Spin','Callback',@Spin);
mViewGrid = uimenu(mView,'Label','Grid','Callback',@ToggleGrid);
if isfield(opt.View,'Grid') & opt.View.Grid
  ToggleGrid
end
if isfield(tok,'rg') & isfield(tok,'zg')
  mViewGU = uimenu(mView,'Label','Grid units','Callback',@ToggleGridUnits);
  if isfield(opt.View,'GridUnits') & opt.View.GridUnits
    ToggleGridUnits
  end
end
mViewNavigate = uimenu(mView,'Label','Navigate...','Callback',@Navigate);

try
  bgrgb = eval(opt.View.whitebg);
catch
  try
    bgrgb = opt.View.whitebg(1:3);
  catch
    bgrgb = [1 1 1];
  end
end
ibg = 1; % This index is used in change_bg
mViewBG = uimenu(mView,'Label','whitebg','Separator','on');
for i = 1:ncolors
  str = char(colorSelection(i));
  mViewBGs(i) = uimenu(mViewBG,'Label',str,'userdata',i,'Callback',@change_bg);
  bgrgbs(i,:) = eval(str);
  if bgrgb == bgrgbs(i,:)
    ibg = i;
    set(mViewBGs(i),'Check','on')
  end
  if strcmp(char(colorSelection(i)),'bronze')
    set(mViewBGs(i),'Separator','on')
  end
end
[A,B] = color_complement(bgrgb);

flytimes = [5 10 30 60 300 600];
if isfield(opt.View,'Flytime') & ~isempty(opt.View.Flytime) & abs(opt.View.Flytime(1)) > 0
  flytime = opt.View.Flytime(1);
else
  flytime = 10; % default
end
iflytime = find(flytime == flytimes);
mViewFT = uimenu(mView,'Label','Fly time');
for i = 1:numel(flytimes)
  mViewFTs(i) = uimenu(mViewFT,'Label',num2str(flytimes(i)),'userdata',i,'Callback',@change_flytime);
  if i == iflytime
    set(mViewFTs(i),'Check','on')
  end
end

spinspeeds = [200 100 50 20 10 5 2 1 -1 -2 -5 -10 -20 -50 -100 -200];
if isfield(opt.View,'Spinrpm') & ~isempty(opt.View.Spinrpm) & abs(opt.View.Spinrpm(1)) > 0
  spinspeed = opt.View.Spinrpm(1);
else
  spinspeed = 20; % default
end
ispin = find(spinspeed == spinspeeds);
mViewSS = uimenu(mView,'Label','Spin rpm');
for i = 1:numel(spinspeeds)
  mViewSSs(i) = uimenu(mViewSS,'Label',num2str(spinspeeds(i)),'userdata',i,'Callback',@change_spinrpm);
  if i == ispin
    set(mViewSSs(i),'Check','on')
  end
end



% Lines
mLin = uimenu('Label','Lines');

if ~isempty(id.LIM)
  mLinLIM = uimenu(mLin,'Label','Limiter','Callback',@ToggleLines,'userdata',id.LIM);
  if ~isfield(opt.Lines,'LIM') | opt.Lines.LIM
    set(mLinLIM,'Check','on')
    Fline(id.LIM) = true;
  end
end
if ~isempty(id.EC)
  mLinEC = uimenu(mLin,'Label','E coil','Callback',@ToggleLines,'userdata',id.EC);
  if ~isfield(opt.Lines,'EC') | opt.Lines.EC
    set(mLinEC,'Check','on')
    Fline(id.EC) = true;
  end
end
if ~isempty(id.FC)
  mLinFC = uimenu(mLin,'Label','F coils','Callback',@ToggleLines,'userdata',id.FC);
  if ~isfield(opt.Lines,'FC') | opt.Lines.FC
    set(mLinFC,'Check','on')
    Fline(id.FC) = true;
  end
end
if ~isempty(id.VV)
  mLinVV = uimenu(mLin,'Label','Vessel','Callback',@ToggleLines,'userdata',id.VV);
  if ~isfield(opt.Lines,'VV') | opt.Lines.VV
    set(mLinVV,'Check','on')
    Fline(id.VV) = true;
  end
end
if ~isempty(id.FL)
  mLinFL = uimenu(mLin,'Label','Flux loops','Callback',@ToggleLines,'userdata',id.FL);
  if isfield(opt.Lines,'FL') & opt.Lines.FL
    set(mLinFL,'Check','on')
    Fline(id.FL) = true;
  end
end
if ~isempty(id.LV)
  mLinLV = uimenu(mLin,'Label','Voltage loops','Callback',@ToggleLines,'userdata',id.LV);
  if isfield(opt.Lines,'LV') & opt.Lines.LV
    set(mLinLV,'Check','on')
    Fline(id.LV) = true;
  end
end
if ~isempty(id.BP)
  mLinBP = uimenu(mLin,'Label','B probes','Callback',@ToggleLines,'userdata',id.BP);
  if isfield(opt.Lines,'BP') & opt.Lines.BP
    set(mLinBP,'Check','on')
    Fline(id.BP) = true;
  end
end
if ~isempty(id.MSE)
  mLinMSE = uimenu(mLin,'Label','MSE points','Callback',@ToggleLines,'userdata',id.MSE);
  if isfield(opt.Lines,'MSE') & opt.Lines.MSE
    set(mLinMSE,'Check','on')
    Fline(id.MSE) = true;
  end
end
if ~isempty(id.TD)
  mLinTD = uimenu(mLin,'Label','3D coils','Callback',@ToggleLines,'userdata',id.TD);
  if isfield(opt.Lines,'TD') & opt.Lines.TD
    set(mLinTD,'Check','on')
    Fline(id.TD) = true;
  end
end
if ~isempty(id.SL)
  mLinSL = uimenu(mLin,'Label','Saddle loops','Callback',@ToggleLines,'userdata',id.SL);
  if isfield(opt.Lines,'SL') & opt.Lines.SL
    set(mLinSL,'Check','on')
    Fline(id.SL) = true;
  end
end

% Line width
linewidths = 1:10;
if isfield(opt.Lines,'Width') & ~isempty(opt.Lines.Width) & opt.Lines.Width(1) > 0
  linewidth = opt.Lines.Width(1);
else
  linewidth = 1; % default
end
iline = find(linewidth == linewidths);
mLinLW = uimenu(mLin,'Label','Line width','Separator','on');
for i = 1:numel(linewidths)
  mLinLWs(i) = uimenu(mLinLW,'Label',num2str(linewidths(i)),'userdata',i,'Callback',@change_linewidth);
  if i == iline
    set(mLinLWs(i),'Check','on')
  end
end



% Surfaces
mSurf = uimenu('Label','Surfaces');

if ~isempty(id.LIM)
  mSurfLIM = uimenu(mSurf,'Label','Limiter','userdata',id.LIM,'Callback',@ToggleSurfaces);
  if ~isfield(opt.Surfaces,'LIM') | opt.Surfaces.LIM
    set(mSurfLIM,'Check','on')
    Fsurf(id.LIM) = true;
  end
end
if ~isempty(id.EC)
  mSurfEC = uimenu(mSurf,'Label','E coil','userdata',id.EC,'Callback',@ToggleSurfaces);
  if ~isfield(opt.Surfaces,'EC') | opt.Surfaces.EC
    set(mSurfEC,'Check','on')
    Fsurf(id.EC) = true;
  end
end
if ~isempty(id.FC)
  mSurfFC = uimenu(mSurf,'Label','F coils','userdata',id.FC,'Callback',@ToggleSurfaces);
  if ~isfield(opt.Surfaces,'FC') | opt.Surfaces.FC
    set(mSurfFC,'Check','on')
    Fsurf(id.FC) = true;
  end
end
if ~isempty(id.VV)
  mSurfVV = uimenu(mSurf,'Label','Vessel','userdata',id.VV,'Callback',@ToggleSurfaces);
  if ~isfield(opt.Surfaces,'VV') | opt.Surfaces.VV
    set(mSurfVV,'Check','on')
    Fsurf(id.VV) = true;
  end
end
if ~isempty(id.FL)
  mSurfFL = uimenu(mSurf,'Label','Flux loops','userdata',id.FL,'Callback',@ToggleSurfaces);
  if isfield(opt.Surfaces,'FL') & opt.Surfaces.FL
    set(mSurfFL,'Check','on')
    Fsurf(id.FL) = true;
  end
end
if ~isempty(id.LV)
  mSurfLV = uimenu(mSurf,'Label','Voltage loops','userdata',id.LV,'Callback',@ToggleSurfaces);
  if isfield(opt.Surfaces,'LV') & opt.Surfaces.LV
    set(mSurfLV,'Check','on')
    Fsurf(id.LV) = true;
  end
end

EdgeColors = 0:10;
if isfield(opt.Surfaces,'EdgeColor') & ~isempty(opt.Surfaces.EdgeColor)
  EdgeColor = max(0,min(10,opt.Surfaces.EdgeColor(1)));
else
  EdgeColor = 0; % default
end
iEdgeColor = find(EdgeColor == EdgeColors);
mSurfLineC = uimenu(mSurf,'Label','Edge color','Separator','on');
for i = 1:numel(EdgeColors)
  mSurfLineCs(i) = uimenu(mSurfLineC,'Label',num2str(EdgeColors(i)),...
    'userdata',i,'Callback',@change_EdgeColor);
  if i == iEdgeColor
    set(mSurfLineCs(i),'Check','on')
  end
end

FaceColors = 0:10;
if isfield(opt.Surfaces,'FaceColor') & ~isempty(opt.Surfaces.FaceColor)
  FaceColor = opt.Surfaces.FaceColor(1);
else
  FaceColor = 5; % default
end
iFaceColor = find(FaceColor == FaceColors);
mSurfBright = uimenu(mSurf,'Label','Face color');
for i = 1:numel(FaceColors)
  mSurfBrights(i) = uimenu(mSurfBright,'Label',num2str(FaceColors(i)),...
    'userdata',i,'Callback',@change_FaceColor);
  if i == iFaceColor
    set(mSurfBrights(i),'Check','on')
  end
end

Dphis = [360 270 180 90 -90 -180 -270];
if isfield(opt.Surfaces,'Dphi') & ~isempty(opt.Surfaces.Dphi)
  Dphi = opt.Surfaces.Dphi(1);
else
  Dphi = 180; % default
end
iDphi = find(Dphi == Dphis);
mSurfDphi = uimenu(mSurf,'Label',char([916 966]));
for i = 1:numel(Dphis)
  mSurfDphis(i) = uimenu(mSurfDphi,'Label',num2str(Dphis(i)),'userdata',i,'Callback',@change_Dphi);
  if i == iDphi
    set(mSurfDphis(i),'Check','on')
  end
end

dphis = [1 2 5 10 22.5];
if isfield(opt.Surfaces,'dphi') & ~isempty(opt.Surfaces.dphi)
  dphi = opt.Surfaces.dphi(1);
else
  dphi = 10; % default
end
idphi = find(dphi == dphis);
mSurfdphi = uimenu(mSurf,'Label',char([948 966]));
for i = 1:numel(dphis)
  mSurfdphis(i) = uimenu(mSurfdphi,'Label',num2str(dphis(i)),'userdata',i,'Callback',@change_dphi);
  if i == idphi
    set(mSurfdphis(i),'Check','on')
  end
end



% Names
mNames = uimenu('Label','Names');
if isfield(tok,'ecnames')
  mNamesEC = uimenu(mNames,'Label','ecnames','Callback',@ToggleNames,'userdata',id.EC);
  if isfield(opt.Names,'EC') & opt.Names.EC
    set(mNamesEC,'Check','on')
    Fname(id.EC) = true;
  end
end
if isfield(tok,'fcnames')
  mNamesFC = uimenu(mNames,'Label','fcnames','Callback',@ToggleNames,'userdata',id.FC);
  if isfield(opt.Names,'FC') & opt.Names.FC
    set(mNamesFC,'Check','on')
    Fname(id.FC) = true;
  end
end
if isfield(tok,'vvnames')
  mNamesVV = uimenu(mNames,'Label','vvnames','Callback',@ToggleNames,'userdata',id.VV);
  if isfield(opt.Names,'VV') & opt.Names.VV
    set(mNamesVV,'Check','on')
    Fname(id.VV) = true;
  end
end
if isfield(tok,'flnames')
  mNamesFL = uimenu(mNames,'Label','flnames','Callback',@ToggleNames,'userdata',id.FL);
  if isfield(opt.Names,'FL') & opt.Names.FL
    set(mNamesFL,'Check','on')
    Fname(id.FL) = true;
  end
end
if isfield(tok,'lvnames')
  mNamesLV = uimenu(mNames,'Label','lvnames','Callback',@ToggleNames,'userdata',id.LV);
  if isfield(opt.Names,'LV') & opt.Names.LV
    set(mNamesLV,'Check','on')
    Fname(id.LV) = true;
  end
end
if isfield(tok,'bpnames')
  mNamesBP = uimenu(mNames,'Label','bpnames','Callback',@ToggleNames,'userdata',id.BP);
  if isfield(opt.Names,'BP') & opt.Names.BP
    set(mNamesBP,'Check','on')
    Fname(id.BP) = true;
  end
end
if isfield(tok,'msenames')
  mNamesMSE = uimenu(mNames,'Label','msenames','Callback',@ToggleNames,'userdata',id.MSE);
  if isfield(opt.Names,'MSE') & opt.Names.MSE
    set(mNamesMSE,'Check','on')
    Fname(id.MSE) = true;
  end
end
if isfield(tok,'tdnames')
  mNamesTD = uimenu(mNames,'Label','tdnames','Callback',@ToggleNames,'userdata',id.TD);
  if isfield(opt.Names,'TD') & opt.Names.TD
    set(mNamesTD,'Check','on')
    Fname(id.TD) = true;
  end
end
if isfield(tok,'slnames')
  mNamesSL = uimenu(mNames,'Label','slnames','Callback',@ToggleNames,'userdata',id.SL);
  if isfield(opt.Names,'SL') & opt.Names.SL
    set(mNamesSL,'Check','on')
    Fname(id.SL) = true;
  end
end

% Font size
fontsizes = [4 6 8 10 12 16 20 24 28 32 40 48];
if isfield(opt.Names,'FontSize') & ~isempty(opt.Names.FontSize) & opt.Names.FontSize(1) > 0
  fontsize = opt.Names.FontSize(1);
else
  fontsize = 12; % default
end
ifont = find(fontsize == fontsizes);
mNamesFS = uimenu(mNames,'Label','Font size','Separator','on');
for i = 1:numel(fontsizes)
  mNamesFSs(i) = uimenu(mNamesFS,'Label',num2str(fontsizes(i)),'userdata',i,'Callback',@change_fontsize);
  if i == ifont
    set(mNamesFSs(i),'Check','on')
  end
end
mNamesSearch = uimenu(mNames,'Label','Search...','Separator','on','Callback',@SearchWindow);

% Default camera
defcam = [];

% Update surfaces
UpdateSurfaces

% The background color
whitebg(gcf,bgrgb)

% Begin by displaying the 2D view
view2D
  
if isfield(opt.View,'is3D') & opt.View.is3D
  view3D
end

% Initialize the view to any user specification
if isfield(opt,'az') & isfield(opt,'el')
  view(opt.az,opt.el)
end

% Initialize the axis to any user specification
if isfield(opt,'axis')
  a = axis;
  n = min(length(a),length(opt.axis));
  a(1:n) = opt.axis(1:n);
  axis(a)
end

% When bgfixed is true the background image can't be changed in the GUI
bgfixed = false;

% Display user-defined background image
if isfield(opt,'imbg')
  if ischar(opt.imbg)
    try
      imbg = imread(opt.imbg);
    catch
      imbg = [];
    end
  else
    imbg = opt.imbg;
  end
  if ishandle(hbg)
    delete(hbg)
  end
  hbg = axes('units','normalized', 'position',[0 0 1 1]);
  uistack(hbg,'bottom');
  try
    image(imbg)
  end
  set(hbg,'handlevisibility','off','visible','off')
end
if isfield(opt,'imbgf')
  if ischar(opt.imbgf)
    try
      imbg = imread(opt.imbgf);
    catch
      imbg = [];
    end
  else
    imbg = opt.imbgf;
  end
  if ishandle(hbg)
    delete(hbg)
  end
  hbg = axes('units','normalized', 'position',[0 0 1 1]);
  uistack(hbg,'bottom');
  try
    image(imbg)
    bgfixed = true;
  end
  set(hbg,'handlevisibility','off','visible','off')
end

% Fly if requested
if isfield(opt.View,'Fly') & opt.View.Fly(1) ~= 0
  nfly = opt.View.Fly(1);
  for ii = 1:round(nfly)
    Fly
  end
end

% Finally, perform a spin if requested in opt
if isfield(opt.View,'Spin') & opt.View.Spin(1) ~= 0
  nspin = opt.View.Spin(1);
  Spin
end

% Number of spins
nspin = 1;

% Everything done, the contents of opt no longer matters, only user interaction matters from now
% Now it's up to the user to cause more code execution by clicking in the menu, etc
% The following functions are inside the function tokgeo and therefore share variables

% Callback functions
function view2D(source,event)
  
  % Remember handles to searched objects so they can be restored for the new axis
  sFname = ishandle(sHname);
  
  % If the 3D view is rotated to almost the 2D view before executing view(0,90) then double clicking zoom
  % rotates the view away from view(2). This is a bug. The workaround is to erase the axis and redraw everything.
  if ishandle(h)
    delete(h)
  end
  
  % Create axis
  h = axes;
  hold on
  if isgrid
    grid on
  end

  % Update check marks
  set(mView2D,'Check','on')
  set(mView3D,'Check','off')
  is2D = true;
  
  DrawLines
  
  DisplayNames
  
  % Restore search hits
  for i = find(sFname)
    if gridunits
      sHname(i) = text(h,(Rname(i)-rgmin)/dr+1,(Zname(i)-zgmin)/dz+1,Sname{i});
    else
      sHname(i) = text(h,Rname(i),Zname(i),Sname{i});
    end
    set(sHname(i),'FontSize',fontsize,'Color',A+B.*[0.8 0 0],'HorizontalAlignment','Center');
  end
  
  axis image
  
  if gridunits
    axis(a2Dgu)
  else
    axis(a2D)
  end
  
  if isempty(texth)
    try
      dum = text(h,a2D(1),a2D(3),'TEST');
      delete(dum)
      texth = true;
    catch
      texth = false; % Can't use handle as first input to text
    end
  end
  
  zoom on
  
  title(tokamak)
      
end

function view3D(source,event)

  % Update check marks
  set(mView2D,'Check','off')
  set(mView3D,'Check','on')
  is2D = false;
  
  % R must be zero where it really is zero, hence only real units in 3D view
  if gridunits
    ToggleGridUnits
  end
  
  axis auto
  
  DrawLines
  
  DrawSurfaces
  
  DisplayNames

  view(-30,45)
  rotate3d on

  axis equal
  
  % Store default camera settings in defcam
  defcam = getcam;
    
end

% Fly through vessel.
% It takes 1/3 of the time to get to the center of the limiter cross section at phi=0
% It takes 1/3 of the time to go Dphi degrees through the vessel at constant speed
% It takes 1/3 of the time to return to the starting point
% No discontinuities in the velocity or angle being viewed
% Acceleration is constant within the three intervals
function Fly(source,event)
  
  if is2D
    return
  end
  
  % Time when flight begins [sec]
  t0 = now*24*3600;

  % Camera target
  trg = get(h,'CameraTarget');

  % Starting and final point
  xyz0 = get(h,'CameraPosition');

  % Camera view angle
  vang = get(h,'CameraViewAngle');

  % Change the view angle to 20 for a good viewing experience inside the vessel
  set(h,'CameraViewAngle',20)

  % Move the camera so that the initial image doesn't change much
  xyz = trg + vang/20*(xyz0-trg);
  set(h,'CameraPosition',xyz)

  % Entry into vessel is always at phi = 0
  rlmin = min(Rline{id.LIM});
  rlmax = max(Rline{id.LIM});
  zlmin = min(Zline{id.LIM});
  zlmax = max(Zline{id.LIM});
  r0 = (rlmax+rlmin)/2;
  z0 = (zlmax+zlmin)/2;

  % Dphi in radians
  dp = Dphi*pi/180;

  % The flight takes 3 time units
  % At t=0, the position is xyz and velocity is zero and R = R0+R1*t+R2*t^2+R3*t^3
  % At t=1, the position is [r0,0,z0] and velocity is r0*dp*[0,1,0] and phi = dp*(t-1)
  % At t=2, the position is [r0*cos(dp),r0*sin(dp),z0] and velocity r0*dp*[-sin(dp),cos(dp),0]

  % Polynomial coefficients for correct positions and velocities
  xA = [0 0 0 1; 1 1 1 1; 0 0 1 0; 3 2 1 0]\[xyz(1); r0; 0; 0];
  yA = [0 0 0 1; 1 1 1 1; 0 0 1 0; 3 2 1 0]\[xyz(2); 0; 0; r0*dp];
  zA = [0 0 0 1; 1 1 1 1; 0 0 1 0; 3 2 1 0]\[xyz(3); z0; 0; 0];

  % Polynomial coefficients for start, final position, and acceleration=0 at t=1 at right speed
  xA = [0 0 0 1; 1 1 1 1; 6 2 0 0; 3 2 1 0]\[xyz(1); r0; 0; 0];
  yA = [0 0 0 1; 1 1 1 1; 6 2 0 0; 3 2 1 0]\[xyz(2); 0; 0; r0*dp];
  zA = [0 0 0 1; 1 1 1 1; 6 2 0 0; 3 2 1 0]\[xyz(3); z0; 0; 0];

  xC = [8 4 2 1; 27 9 3 1; 12 4 1 0; 27 6 1 0]\[r0*cos(dp); xyz(1); -r0*dp*sin(dp); 0];
  yC = [8 4 2 1; 27 9 3 1; 12 4 1 0; 27 6 1 0]\[r0*sin(dp); xyz(2); +r0*dp*cos(dp); 0];
  zC = [8 4 2 1; 27 9 3 1; 12 4 1 0; 27 6 1 0]\[z0; xyz(3); 0; 0];

  % When at phi, look at the inner limiter surface at phi+alpha
  alpha = acos(rlmin/r0);

  % phi+alpha is the angle where we will be in dt time units
  dt = abs(alpha/dp);

  while now*24*3600 < t0+flytime
    t = (now*24*3600-t0)/flytime*3;
    if t < 1
      x = polyval(xA,t);
      y = polyval(yA,t);
      z = polyval(zA,t);
      vx = polyval([3;2;1].*xA(1:3),t);
      vy = polyval([3;2;1].*yA(1:3),t);
      vz = polyval([3;2;1].*zA(1:3),t);
    elseif t < 2
      x = r0*cos(dp*(t-1));
      y = r0*sin(dp*(t-1));
      z = z0;
      vx = -r0*sin(dp*(t-1))*dp;
      vy = +r0*cos(dp*(t-1))*dp;
      vz = 0;
    else
      x = polyval(xC,t);
      y = polyval(yC,t);
      z = polyval(zC,t);
      vx = polyval([3;2;1].*xC(1:3),t);
      vy = polyval([3;2;1].*yC(1:3),t);
      vz = polyval([3;2;1].*zC(1:3),t);
    end
    set(h,'CameraPosition',[x y z])

    % Camera target
    if t < 1-dt
      f = t/(1-dt);
      x = f*r0 + (1-f)*trg(1);
      y =        (1-f)*trg(2);
      z = f*z0 + (1-f)*trg(3);
    elseif t < 1
      f = (t+dt-1)/dt;
      x = (f*rlmin + (1-f)*r0)*cos(dp*(t-1+dt));
      y = (f*rlmin + (1-f)*r0)*sin(dp*(t-1+dt));
      z = f*zlmax + (1-f)*z0;
    elseif t < 2-dt
      f = (t-1)/(1-dt);
      x = rlmin*cos(dp*(t-1+dt));
      y = rlmin*sin(dp*(t-1+dt));
      z = f*zlmin + (1-f)*zlmax;
    elseif t < 2
      x = rlmin*cos(dp);
      y = rlmin*sin(dp);
      z = zlmin;
    else
      f = (t-2);
      x = f*trg(1) + (1-f)*rlmin*cos(dp);
      y = f*trg(2) + (1-f)*rlmin*sin(dp);
      z = f*trg(3) + (1-f)*zlmin;
    end

    set(h,'CameraTarget',[x y z])
    set(h,'CameraViewAngle',20)
    drawnow
    if is2D | ~ishandle(h)
      return % User selected View 2D and function view2D executed, or figure has been deleted
    end
  end  
  set(h,'CameraPosition',xyz)
  set(h,'CameraViewAngleMode','auto')
end

% During the spin, the command window will not echo characters
% It doesn't help to pause in while loops and I have no other ideas
% An impatient user will need to press control-c to see the echo right away
function Spin(source,event)
  
  if is2D
    return
  end
  
  % To prevent size changes during spin, CameraViewAngleMode is set to manual
  % Alas this also overrides the "stretch-to-fill behavior" in auto mode
  % Therefore size changes can occur by switching between manual and auto
  set(h,'CameraViewAngleMode','manual')
  
  spintime = nspin*60/abs(spinspeed); % spinspeed is in rpm
  t0 = now*24*3600;
  trg = get(h,'CameraTarget');
  xyz = get(h,'CameraPosition');
  ang = angle(xyz(1)+1i*xyz(2));
  rho = sqrt((xyz(1)-trg(1))^2+(xyz(2)-trg(2))^2);
  xfinal = rho*cos(nspin*2*pi+ang);
  yfinal = rho*sin(nspin*2*pi+ang);
  while now*24*3600 < t0+spintime
    phi = (now*24*3600-t0)*spinspeed/60*2*pi; % radians
    xyz = get(h,'CameraPosition');
    set(h,'CameraPosition',[trg(1)+rho*cos(phi+ang) trg(2)+rho*sin(phi+ang) xyz(3)])
    drawnow
    if is2D | ~ishandle(h)
      return % User selected View 2D and function view2D executed, or figure has been deleted
    end
  end  
  set(h,'CameraPosition',[xfinal yfinal xyz(3)])
  set(h,'CameraViewAngleMode','auto')
end

function ToggleGrid(source,event)
  isgrid = ~isgrid;
  if isgrid
    set(mViewGrid,'Check','on')
  else
    set(mViewGrid,'Check','off')
  end
  if ishandle(h)
    if isgrid
      grid on
    else
      grid off
    end
  end
end

function ToggleGridUnits(source,event)
  gridunits = ~gridunits & is2D;
  if gridunits
    set(mViewGU,'Check','on')
  else
    set(mViewGU,'Check','off')
  end
  if is2D
    view2D % Delete h and redraw, otherwise double-clicking zoom won't always work
  end
end

function Navigate(source,event)
  
  if is2D
    return
  end

  cam0 = getcam;
  
  W = 450;
  H = 450;
  scr = get(0,'ScreenSize');
  pos = get(gcf,'Position');
  xpos = max(scr(3)/3,min(pos(1)+pos(3)+10, scr(3)*2/3-W));
  ypos = pos(2)+pos(4)-H;
  d = dialog('Position',[xpos ypos W H],'Name','Go to position and look at target');

  % Initialize
  vmin = [0      -180 -4*Rmax -4*Rmax Zmin-3*DZ];
  vmax = [4*Rmax +180 +4*Rmax +4*Rmax Zmax+3*DZ];
  vang = get(h,'CameraViewAngle');
  v = cam0.CameraPosition;
  v1 = [sqrt(v(1)^2+v(2)^2), angle(v(1)+1i*v(2))*180/pi, v];
  v = cam0.CameraTarget;
  v2 = [sqrt(v(1)^2+v(2)^2), angle(v(1)+1i*v(2))*180/pi, v];
  V = [v1(:) v2(:)];
  
  % The text "CameraViewAngle"
  ta = uicontrol(d,'Style','text','Position',[35 420 150 22],...
	'HorizontalAlignment','Left','FontWeight','Bold','String','CameraViewAngle');

  % The text description
  titO = uicontrol(d,'Style','text','Position',[5 400 25 22],...
    'HorizontalAlignment','Center','FontWeight','Bold','String',char(937));

  % The edit box with the value
  valO = uicontrol(d,'Style','edit','Position',[35 400 100 22],...
    'HorizontalAlignment','Left','String',num2str(vang),'userdata',[0 0],'Callback',@Camera);

  % The slider with the same value as in the edit box
  sliO = uicontrol(d,'Style','Slider','Position',[140 400 275 22],'Value',vang,...
    'SliderStep',[5e-3 1e-2],...
    'Min',1,'Max',45,...
    'userdata',[0 0],'Callback',@Camera);

  % The text "CameraPosition"
  tp = uicontrol(d,'Style','text','Position',[35 350 200 22],...
	'HorizontalAlignment','Left','FontWeight','Bold','String','CameraPosition');

  % The text "CameraTarget"
  tt = uicontrol(d,'Style','text','Position',[35 175 200 22],...
	'HorizontalAlignment','Left','FontWeight','Bold','String','CameraTarget');

  tits = ['R' char(966) 'xyz'];
  for i = 1:5
    for j = 1:2
    
      % The value
      V(i,j) = min(vmax(i),max(vmin(i),V(i,j)));
    
      % Position of GUI elements
      y = 525-25*i-175*j;
      
      % The text description for each line
      tit(i,j) = uicontrol(d,'Style','text','Position',[5 y 25 22],...
	'HorizontalAlignment','Center','FontWeight','Bold','String',tits(i));
      
      % The edit box with the value
      val(i,j) = uicontrol(d,'Style','edit','Position',[35 y 100 22],...
	'HorizontalAlignment','Left','String',num2str(V(i,j)),'userdata',[i j],'Callback',@Camera);

      % The slider with the same value as in the edit box
      sli(i,j) = uicontrol(d,'Style','Slider','Position',[140 y 275 22],'Value',V(i,j),...
        'SliderStep',[5e-3 1e-2],...
        'Min',vmin(i),'Max',vmax(i),...
        'userdata',[i j],'Callback',@Camera);
      
    end
  end

  ok  = uicontrol(d,'Position',[035 10 100 25],'String','Ok',     'FontWeight','Bold','Callback',@Ok);
  def = uicontrol(d,'Position',[175 10 100 25],'String','Default','FontWeight','Bold','Callback',@Defaults);
  can = uicontrol(d,'Position',[315 10 100 25],'String','Cancel', 'FontWeight','Bold','Callback',@Cancel);
  
  % Focus on Ok
  uicontrol(ok)

  function Camera(source,event)
    if strcmp(get(source,'Style'),'edit')
      str = get(source,'String');
      num = str2num(str);
    else
      num = get(source,'Value');
    end
    if isempty(num)
      set(source,'String','')
      return
    end
    ij = get(source,'userdata');
    i = ij(1);
    j = ij(2);
    if i == 0 & j == 0
      num = min(get(sliO,'Max'),max(get(sliO,'Min'),num));
      set(valO,'String',num2str(num))
      set(sliO,'Value',num)
      set(h,'CameraViewAngle',num)
      return
    end
    if strcmp(get(source,'Style'),'edit')
      R = str2num(get(val(1,j),'String'));
      p = str2num(get(val(2,j),'String'));
      x = str2num(get(val(3,j),'String'));
      y = str2num(get(val(4,j),'String'));
      z = str2num(get(val(5,j),'String'));
    else % Slider
      R = get(sli(1,j),'Value');
      p = get(sli(2,j),'Value');
      x = get(sli(3,j),'Value');
      y = get(sli(4,j),'Value');
      z = get(sli(5,j),'Value');
    end
    rmax = get(sli(i,j),'Max');
    if i == 1
      x = R*cos(p*pi/180);
      y = R*sin(p*pi/180);
    elseif i == 2
      x = R*cos(p*pi/180);
      y = R*sin(p*pi/180);
      % Wrap around logic for slider
      if p == 180
        set(sli(i,j),'Max',+180.001) % Make it possible to go 1 step more
      elseif p > 180
        set(sli(i,j),'Max',+180)     % Restore upper limit
	p = -180;                    % Wrap around
        set(sli(i,j),'Min',-180.001) % Make it possible to turn back
      elseif p == -180
        set(sli(i,j),'Min',-180.001) % Make it possible to go 1 step more
      elseif p < -180
        set(sli(i,j),'Min',-180)     % Restore lower limit
	p = 180;                     % Wrap around
        set(sli(i,j),'Max',+180.001) % Make it possible to turn back
      else
        set(sli(i,j),'Min',-180)
        set(sli(i,j),'Max',+180)
      end
    elseif i == 3
      xmax = sqrt(rmax^2-y^2);
      % Bounce back the x slider if it makes R too big
      x = max(-xmax,min(xmax,x));
      R = sqrt(x^2+y^2);
      p = angle(x+1i*y)*180/pi;
    elseif i == 4
      ymax = sqrt(rmax^2-x^2);
      % Bounce back the y slider if it makes R too big
      y = max(-ymax,min(ymax,y));
      R = sqrt(x^2+y^2);
      p = angle(x+1i*y)*180/pi;
    end
    if abs(x) < 1e-13
      x = 0;
    end
    if abs(y) < 1e-13
      y = 0;
    end
    if j == 1
      set(h,'CameraPosition',[x y z])
    elseif j == 2
      set(h,'CameraTarget',[x y z])
    end
    pp = nan(1,2);
    pp(j) = p;
    UpateNavControls(getcam,j,[R p x y z])
  end
  
  function UpateNavControls(cam,k,vk)
    vO = cam.CameraViewAngle;
    vO = min(get(sliO,'Max'),max(get(sliO,'Min'),vO));
    set(valO,'String',num2str(vO))
    set(sliO,'Value',vO)
    v(3:5,1) = cam.CameraPosition;
    v(3:5,2) = cam.CameraTarget;
    v(1,:) = sqrt(v(3,:).^2+v(4,:).^2);
    v(2,:) = angle(v(3,:) + 1i*v(4,:))*180/pi;
    if nargin > 2
      v(:,k) = vk;
    end
    for i = 1:5
      for j = 1:2
        v(i,j) = min(get(sli(i,j),'Max'),max(get(sli(i,j),'Min'),v(i,j)));
        set(val(i,j),'String',num2str(v(i,j)))
        set(sli(i,j),'Value',v(i,j))
      end
    end
  end
  
  function Ok(source,event)
    delete(d)
  end
  
  function Defaults(source,event)
    setcam(defcam)
    UpateNavControls(getcam)
  end
  
  function Cancel(source,event)
    setcam(cam0)
    set(h,'CameraViewAngleMode',cam0.CameraViewAngleMode)
    delete(d)
  end
  
end

function change_bg(source,event)
  set(mViewBGs(ibg),'Check','off')
  ibg = get(source,'userdata');
  bgrgb = bgrgbs(ibg,:);
  set(mViewBGs(ibg),'Check','on')
  whitebg(gcf,bgrgb)
  [A,B] = color_complement(bgrgb);
  for i = 1:N
    if ishandle(Hline(i))
      set(Hline(i),'Color',A+B.*Cline(i,:))
    end
    if ishandle(Hsurf(i))
      set(Hsurf(i),'FaceColor',A+B.*Csurf(i,:),'EdgeColor',A+B.*Esurf(i,:))
    end
  end
  if bgfixed
    return
  end
  if ishandle(hbg)
    delete(hbg)
  end
  if imstat(ibg) == 0
    try
      im(ibg) = {imread([char(colorSelection(ibg)) 'bg.jpg'])};
      if isa(im{ibg},'uint8') & size(im{ibg},3) == 3
        imstat(ibg) = 2;
      else
        imstat(ibg) = 1;
      end
    catch
      imstat(ibg) = 1;
    end
  end
  if imstat(ibg) == 2
    hbg = axes('units','normalized', 'position',[0 0 1 1]);
    uistack(hbg,'bottom');
    image(im{ibg})
    set(hbg,'handlevisibility','off','visible','off')
  end
end

function change_flytime(source,event)
  set(mViewFTs(iflytime),'Check','off')
  iflytime = get(source,'userdata');
  flytime = flytimes(iflytime);
  set(mViewFTs(iflytime),'Check','on')
end

function change_spinrpm(source,event)
  set(mViewSSs(ispin),'Check','off')
  ispin = get(source,'userdata');
  spinspeed = spinspeeds(ispin);
  set(mViewSSs(ispin),'Check','on')
end

function ToggleLines(source,event)
  % Get indices
  k = get(source,'userdata');
  % Toggle the check mark
  if strcmp(get(source,'Check'),'off')
    set(source,'Check','on')
    Fline(k) = true;
  else
    set(source,'Check','off')
    Fline(k) = false;
  end
  DrawLines
end

function DrawLines
  for i = 1:N
    if ishandle(Hline(i))
      delete(Hline(i))
    end
    if Fline(i)
      % Draw object for 2D or 3D viewing
      if is2D
        if gridunits
          Hline(i) = plot((Rline{i}-rgmin)/dr+1,(Zline{i}-zgmin)/dz+1);
	else
          Hline(i) = plot(Rline{i},Zline{i});
	end
      else
        Hline(i) = plot3(Xline{i},Yline{i},Zline{i});
      end
      set(Hline(i),'LineWidth',linewidth,'Color',A+B.*Cline(i,:))
    end
  end
end

function change_linewidth(source,event)
  set(mLinLWs(iline),'Check','off')
  iline = get(source,'userdata');
  linewidth = linewidths(iline);
  set(mLinLWs(iline),'Check','on')
  for i = 1:N
    if ishandle(Hline(i))
      set(Hline(i),'LineWidth',linewidth)
    end
  end
end

function ToggleSurfaces(source,event)
  k = get(source,'userdata');
  if strcmp(get(source,'Check'),'on')
    set(source,'Check','off')
    Fsurf(k) = false;
  else
    set(source,'Check','on')
    Fsurf(k) = true;
  end
  DrawSurfaces
end

function DrawSurfaces
  for i = 1:N
    if ishandle(Hsurf(i))
      delete(Hsurf(i))
    end
    if Fsurf(i) & ~is2D
      Hsurf(i) = surf(Xsurf{i},Ysurf{i},Zsurf{i},'FaceColor',A+B.*Csurf(i,:),'EdgeColor',A+B.*Esurf(i,:));
    end
  end
end

function change_EdgeColor(source,event)
  set(mSurfLineCs(iEdgeColor),'Check','off')
  iEdgeColor = get(source,'userdata');
  EdgeColor = EdgeColors(iEdgeColor);
  set(mSurfLineCs(iEdgeColor),'Check','on')
  UpdateSurfaces
  for i = 1:N
    if ishandle(Hsurf(i))
      set(Hsurf(i),'EdgeColor',A+B.*Esurf(i,:))
    end
  end  
end

function change_FaceColor(source,event)
  set(mSurfBrights(iFaceColor),'Check','off')
  iFaceColor = get(source,'userdata');
  FaceColor = FaceColors(iFaceColor);
  set(mSurfBrights(iFaceColor),'Check','on')
  UpdateSurfaces
  for i = 1:N
    if ishandle(Hsurf(i))
      set(Hsurf(i),'FaceColor',A+B.*Csurf(i,:))
    end
  end  
end

function change_Dphi(source,event)
  set(mSurfDphis(iDphi),'Check','off')
  iDphi = get(source,'userdata');
  Dphi = Dphis(iDphi);
  set(mSurfDphis(iDphi),'Check','on')
  UpdateSurfaces
  for i = 1:N
    if ishandle(Hsurf(i))
      set(Hsurf(i),'XData',Xsurf{i})
      set(Hsurf(i),'YData',Ysurf{i})
      set(Hsurf(i),'ZData',Zsurf{i})
      set(Hsurf(i),'CData',Xsurf{i}*0)
    end
  end  
end

function change_dphi(source,event)
  set(mSurfdphis(idphi),'Check','off')
  idphi = get(source,'userdata');
  dphi = dphis(idphi);
  set(mSurfdphis(idphi),'Check','on')
  UpdateSurfaces
  for i = 1:N
    if ishandle(Hsurf(i))
      set(Hsurf(i),'XData',Xsurf{i})
      set(Hsurf(i),'YData',Ysurf{i})
      set(Hsurf(i),'ZData',Zsurf{i})
      set(Hsurf(i),'CData',Xsurf{i}*0)
    end
  end  
end

function ToggleNames(source,event)
  % Get indices
  k = get(source,'userdata');
  % Toggle the check mark
  if strcmp(get(source,'Check'),'off')
    set(source,'Check','on')
    Fname(k) = true;
  else
    set(source,'Check','off')
    Fname(k) = false;
  end
  DisplayNames
end

function DisplayNames
  for i = 1:N
    if ishandle(Hname(i))
      delete(Hname(i))
    end
    if Fname(i)
      % Print names for 2D or 3D viewing
      if is2D
        if gridunits
          Hname(i) = text((Rname(i)-rgmin)/dr+1,(Zname(i)-zgmin)/dz+1,Sname{i});
	else
          Hname(i) = text(Rname(i),Zname(i),Sname{i});
	end
      else
        Hname(i) = text(Xname(i),Yname(i),Zname(i),Sname{i});
      end
      set(Hname(i),'Color',A,'FontSize',fontsize,'HorizontalAlignment','Center');
    end
  end
end

function change_fontsize(source,event)
  set(mNamesFSs(ifont),'Check','off')
  ifont = get(source,'userdata');
  fontsize = fontsizes(ifont);
  set(mNamesFSs(ifont),'Check','on')
  for i = 1:N
    if ishandle(Hname(i))
      set(Hname(i),'FontSize',fontsize)
    end
    if ishandle(sHname(i))
      set(sHname(i),'FontSize',fontsize)
    end
  end
end

function SearchWindow(source,event)
  scr = get(0,'ScreenSize');
  pos = get(gcf,'Position');
  xpos = max(scr(3)/3,min(pos(1)+pos(3)+10, scr(3)*2/3-240));
  ypos = pos(2)+pos(4)-140;
  d = dialog('Position',[xpos ypos 240 140],'Name','Search names');
  txt = uicontrol(d,'Style','text','Position',[20 110 210 22],...
    'HorizontalAlignment','Left','String','String to look for:');
  term = uicontrol(d,'Style','edit','Position',[20 80 210 22],...
    'HorizontalAlignment','Left','TooltipString','Search string','Callback',@Search);
  doit = uicontrol(d,'Position',[35 50 70 25],'String','Search','Callback',@Search);
  endit = uicontrol(d,'Position',[135 50 70 25],'String','Close','Callback','delete(gcf)');
  forget = uicontrol(d,'Position',[35 20 170 25],'String','Clear search hits','Callback',@clear_hits);
  % Focus on edit box
  uicontrol(term)
end

function Search(source,event)
  str = upper(get(term,'String'));
  set(term,'String','')
  foundnow = false(1,N);
  for i = 1:N
    foundnow(i) = any(strfind(upper(Sname{i}),str));
  end
  if numel(str) >= 3 & strcmp(str(1:3),'VES')
    foundnow(id.VV) = true;
  end
  % Animate vibration of found objects
  XLim = get(h,'XLim');
  YLim = get(h,'YLim');
  ZLim = get(h,'ZLim');
  for dzvib = [1 0 -1 0]*DZ/50
    for i = find(foundnow)
      if ishandle(Hline(i))
        if is2D
	  if gridunits
	    set(Hline(i),'YData',(Zline{i}+dzvib-zgmin)/dz+1)
	  else
	    set(Hline(i),'YData',Zline{i}+dzvib)
	  end
	else
	  set(Hline(i),'ZData',Zline{i}+dzvib)
	end
      end
      if ishandle(Hsurf(i))
	set(Hsurf(i),'ZData',Zsurf{i}+dzvib)
      end
    end
    set(h,'XLim',XLim)
    set(h,'YLim',YLim)
    set(h,'ZLim',ZLim)
    drawnow
  end
  
  % Using h as first input to text if possible because this is better at keeping focus on the edit box  
  if ~texth % This is an old version that doesn't accept handle as first input
    axes(h) % This works back to at least version 2013a where text(h,x,y,z,string) did not work
  end
  for i = find(foundnow)
    if ishandle(sHname(i))
      delete(sHname(i))
    end
    if is2D
      if gridunits
        if texth
          sHname(i) = text(h,(Rname(i)-rgmin)/dr+1,(Zname(i)-zgmin)/dz+1,Sname{i});
	else
          sHname(i) = text((Rname(i)-rgmin)/dr+1,(Zname(i)-zgmin)/dz+1,Sname{i});
	end
      else
        if texth
          sHname(i) = text(h,Rname(i),Zname(i),Sname{i});
	else
          sHname(i) = text(Rname(i),Zname(i),Sname{i});
	end
      end
    else
      sHname(i) = text(h,Xname(i),Yname(i),Zname(i),Sname{i});
    end
    set(sHname(i),'Color',A+B.*[0.8 0 0],'FontSize',fontsize,'HorizontalAlignment','Center')
  end
  
  % Focus back on the edit box called term
  if ~texth
    figure(d)   % Returning focus to term for older matlab version
    pause(0.25) % This works with 2014b on saturn
  end
  uicontrol(term)

end

function clear_hits(source,event)
  for i = 1:N
    if ishandle(sHname(i))
      delete(sHname(i))
    end
  end
  % Focus back on edit box
  uicontrol(term)
end

function [r,z] = data2poly(data)
  if size(data,1) < 6
    data(6,:) = 0;
  end
  dta = data(3,:)./tan(data(6,:)*pi/180);
  dta(data(6,:)==0) = 0;
  dr5 = [-1; -1; +1; +1; -1]*data(4,:)/2 + [-1; +1; +1; -1; -1]*dta/2;
  dta = data(4,:).*tan(data(5,:)*pi/180);
  dz5 = [-1; +1; +1; -1; -1]*data(3,:)/2 + [-1; -1; +1; +1; -1]*dta/2;
  r = ones(5,1)*data(2,:) + dr5;
  z = ones(5,1)*data(1,:) + dz5;
end

function UpdateSurfaces
  for i = [id.LIM id.EC id.FC id.VV id.FL id.LV]
    phi1 = phi0+Dphi;
    r = Rline{i};
    z = Zline{i};
    r = r(:)';
    z = z(:)';
    x = [];
    y = [];
    k = 0;
    for phi = linspace(phi0,phi1,1+ceil(abs(phi1-phi0)/dphi))
      k = k+1; 
      x(k,:) = r*cos(phi*pi/180);
      y(k,:) = r*sin(phi*pi/180);
      z(k,:) = z(1,:);
    end
    Xsurf(i) = {x};
    Ysurf(i) = {y};
    Zsurf(i) = {z};
    if FaceColor > 5
      f = 2*(1-FaceColor/10);
      Csurf(i,:) = f*Cline(i,:) + (1-f);
    else
      f = 2*FaceColor/10;
      Csurf(i,:) = f*Cline(i,:);
    end
    if EdgeColor > 5
      f = 2*(1-EdgeColor/10);
      Esurf(i,:) = f*Cline(i,:) + (1-f);
    else
      f = 2*EdgeColor/10;
      Esurf(i,:) = f*Cline(i,:);
    end
  end
  for i = id.LIM
    Csurf(i,:) = FaceColor/10;
    Esurf(i,:) = EdgeColor/10;
  end
end

function cam = getcam
  cam.CameraPosition      = get(h,'CameraPosition');
  cam.CameraPositionMode  = get(h,'CameraPositionMode');
  cam.CameraTarget        = get(h,'CameraTarget');
  cam.CameraTargetMode    = get(h,'CameraTargetMode');
  cam.CameraUpVector      = get(h,'CameraUpVector');
  cam.CameraUpVectorMode  = get(h,'CameraUpVectorMode');
  cam.CameraViewAngle     = get(h,'CameraViewAngle');
  cam.CameraViewAngleMode = get(h,'CameraViewAngleMode');
end

function setcam(cam)
  set(h,'CameraPosition',cam.CameraPosition)
  %set(h,'CameraPositionMode',cam.CameraPositionMode)
  set(h,'CameraTarget',cam.CameraTarget)
  %set(h,'CameraTargetMode',cam.CameraTargetMode)
  %set(h,'CameraUpVector',cam.CameraUpVector)
  %set(h,'CameraUpVectorMode',cam.CameraUpVectorMode)
  set(h,'CameraViewAngle',cam.CameraViewAngle)
  %set(h,'CameraViewAngleMode',cam.CameraViewAngleMode)
end

% Use this instead of new round function so that old versions of matlab can run
% Same as round(x,2,'significant')
function y = round2(x)
  j = floor(log10(abs(x)));
  y = round(x*10^(1-j))/10^(1-j);
end

% To calculate color to show: A+B.*C where C is Cline or Csurf or Esurf
function [A,B] = color_complement(bgrgb);
  A = [0 0 0];
  B = [1 1 1];
  A(bgrgb==0) = +1;
  B(bgrgb==0) = -1;
end

% End of tokgeo
end



