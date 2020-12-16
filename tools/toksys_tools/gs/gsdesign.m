function eq = gsdesign(spec, init, config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   eq = gsdesign(spec, init, config)
%
%  PURPOSE: Design a 2-D (Grad-Shafranov) equilibrium by minimizing
%           a cost function: sum(weights.param*(param-targets.param))
%
%  INPUTS: spec, specification, a structure with:
%                targets, weights, limits, locks
%
%          init, initial equilibrium
%
%        config, Toksys description of the tokamak, and options
%
%  OUTPUTS: eq, equilibrium 
%
%For detailed information, run gsdesign without inputs
%List of scripts that demonstrate use of gsdesign:
%gsdesign_demo_d3d_DN   - Design of EAST-like DN plasma for DIII-D
%gsdesign_demo_d3d_DSNF - Design of a double snowflake plasma
%See also gsdesign_recreate_init, gsdesign_iss
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Anders Welander  ON	2/15/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How to calculate outputs
persistent outputs index_in_y iy ny user_signal phaseouterrors
persistent fl dfldx nfl mlc mlv mpl
persistent lv dlvdx nlv mhc mhv mph
persistent bp dbpdx nbp gbc gbv gpb 
persistent rog drogdx nrog rldata

% Extra diagnostics
persistent rdp zdp ndp mdc mdv mpd grdc grdv grpd gzdc gzdv gzpd
persistent psipd brdp bzdp btdp dpsipddx dbrdpdx dbzdpdx dbtdpdx
persistent nxpoints kxpoints xpointorder

% The tokamak
persistent tokamak
persistent nr nz ngg rg zg dr dz rgg zgg mpc mpv mpp mcc mcv mvv nc nv ress rss mss
persistent piccc nbbbs_max circle
persistent nlim Rlim Zlim nl rl zl drl dzl dl wl wld wlb wlt wlr wlz iil
persistent ilimgg klimgg rvmin rvmax zvmin zvmax
persistent mrcc mzcc mrpc mzpc

% Options
persistent constraints nkn psikn plotit default_plotit 
persistent converge default_converge max_iterations
persistent evolve_option plot_settings

% For time-dependent & conditional plotting
persistent plot_times plot_counter plot_if_new_response default_time ha

% Remember helping variablesfigure_handle
persistent mu0 twopi R13 mx neighbors ir4 iz4 ir8 iz8 ns nx nxc
persistent sp0 sf0 sg0 c0 c1 c2 c3 d0 d1 d2 d3 d4 
persistent iknotg psibar a b c d aa bb q1 q2 x1 x2 x3
persistent id4 id8 dpsid4 dpsid8 dpsib4 dpsib8 dpsib16 dCldrbbbs dCldzbbbs
persistent M kk initial_plasma plasma turnin concavel halo
persistent indic indiv indis indsp indsf inder AR Ag RA 
persistent fdx fer dzbbbs_max plasma_size ibbbs gbbbs xx wbbbs wds wbs wrs wzs
persistent rhobbbs dpsibbbsdr dpsibbbsdz dpsizrpladpcurrt
persistent drb3dpsia drb3dpsib drb3dpsip dzb3dpsia dzb3dpsib dzb3dpsip
persistent redge zedge fedge gedge iedge xedge thedge rhoedge

% For small-plasma model
persistent rhogmin rhogmax ep

% Equilibrium for which response was calculated
persistent ref lae lar

% For morphing in new cpasma, li, betap, ys in gs_update_eq after gs_eq_analysis
persistent err felip

% Last analyzed equilibrium
persistent iia wa iib wb cpasma maxis rmaxis zmaxis psimag rbdef zbdef psibry
persistent x1exists x1inside ix1 rx1 zx1 iix1 wx1
persistent x2exists x2inside ix2 rx2 zx2 iix2 wx2
persistent itl iitl wtl rbreak zbreak
persistent rzero bzero preszr pprimezr ffprimzr
persistent iused iplasma rbbbs zbbbs zbot ztop psipla li betap pcurrt nbbbs 
persistent fpol ffprim pres pprime psix1 psix2 thbbbs Wth
persistent psizr_err Bp2zr Cl Bp2V Vtot bp2flx fluxerror
persistent p2 p3 f2 f3 r1 r2 r3 z1 z2 z3
persistent rcell zcell Acell ncell
persistent V W A I T jtav rhot qpsi

% Remember present equilibrium
persistent ic iv sp sf er psizr psibarzr

% Plasma response
persistent dydx dysdx dpsipladx
persistent dpsizrdx dAcelldx dpcurrtdx dcpasmadx dlidx dbetapdx lstar
persistent dxdxc dpsizr dpsizrx dpsizr_err
persistent drbbbsdx dzbbbsdx dpcoredpsi dpsizrdpsizrapp
persistent drcurdx dzcurdx drmaxisdx dzmaxisdx drbdefdx dzbdefdx
persistent response_count

% Variables that are used to construct sp0, sf0, sg0
persistent psibry0 psimag0 pprime0 ffprim0
persistent no_edge_current no_edge_gradient

% Parameters used by gs_dynamics to evolve x
persistent Amat Bmat plares rxx Rpla Rpla_previous gamma drcurdv1 dzcurdv1
persistent dxdx3 dx3dx1
persistent eta_vs_rhot sb0
persistent calculate_profiles calculate_profile_responses

% Used by gs_update_eq 
persistent dpsibar_limit_for_linear

% Circuits
persistent icci picci ci nci dxdxs Rhat Vhat nps ncd
persistent Pxx Pxxi Rextra Lextra Rckt_extra Lckt_extra

% For EFITs and cc_efit_to_tok
persistent ecid fcid vvid fcturn ecturn turnfc fcnturn ecnturn vvfrac
persistent idx_efit_to_tok imks iterminal


if nargin == 0
  s = which('gsdesign_help.txt');
  unix(['more ' s]);
  if nargout == 0
    clear eq
  else
    eq = ['This variable will be a structure with designed equilibrium', ...
          10, 'if a specification is given as input to gsdesign.'];
  end
  return
end

% Read or write to persistent variables from caller's ws
if ischar(spec)
  if nargin == 1
    if exist(spec,'var')
      eq = eval(spec);
    elseif strcmp(spec,'eq') % Return a structure eq
      gs_create_struct_eq
    elseif strcmp(upper(spec),'FUN2WS') % All persistent to caller's ws
      vars = who;
      for i = 1:size(vars,1)
	assignin('caller',char(vars(i)),eval(char(vars(i))))
      end
    elseif strcmp(upper(spec),'WS2FUN') % From caller's ws to persistent
      vars = who;
      for i = 1:size(vars,1)
	if evalin('caller',['exist(''' char(vars(i)) ''',''var'')'])
	  eval([char(vars(i)) ' = evalin(''caller'',''' char(vars(i)) ''');'])
	end
      end
    else
      eq = [spec ' is not a persistent variable in gsdesign'];
    end
    return
  elseif nargin == 2
    if exist(spec,'var')
      eval([spec ' = init;']);
      eq = init;
    elseif strcmp(upper(spec),'HELP')
      gshelp(init,'gsdesign')
    else
      disp([spec ' is not a persistent variable in gsevolve']);
      eq = [spec ' is not a persistent variable in gsevolve'];
    end
    if nargout == 0
      clear eq
    end
    return
  end
end

xc = [];

% Configure
if exist('config','var') & isstruct(config)
  % Define diagnostics as outputs when designing
  if isfield(config,'outputs')
    outputs = char('',config.outputs);
  else
    outputs = '';
  end
  if isfield(config,'mlc') & isfield(config,'mlv') & ...
     isfield(config,'mpl') & ~isempty(config.mlc)
    outputs = char(outputs,'fl');
  end
  if isfield(config,'mhc') & isfield(config,'mhv') & ...
     isfield(config,'mph') & ~isempty(config.mhc)
    outputs = char(outputs,'lv');
  end
  if isfield(config,'gbc') & isfield(config,'gbv') & ...
     isfield(config,'gpb') & ~isempty(config.gbc)
    outputs = char(outputs,'bp');
  end
  if isfield(config,'rldata') & ~isempty(config.rldata)
    outputs = char(outputs,'rog');
  end
  if isfield(config,'rdp') & isfield(config,'zdp')
    outputs = char(outputs,'psidp');
    outputs = char(outputs,'brdp');
    outputs = char(outputs,'bzdp');
    outputs = char(outputs,'btdp');
  end
  if size(outputs,1) > 1
    config.outputs = outputs(2:end,:);
  end
  need_mrmz = ...
    isfield(spec,'targets') && ...
    (isfield(spec.targets,'frc') | isfield(spec.targets,'fzc')) | ...
    isfield(spec,'locks') && ...
    (isfield(spec.locks,'frc') | isfield(spec.locks,'fzc')) | ...
    isfield(spec,'limits') && ...
    (isfield(spec.limits,'frc') | isfield(spec.limits,'fzc')) && ...
    (isempty(mrcc) | isempty(mzcc) | isempty(mrpc) | isempty(mzpc));
  gs_configure
  if ~isfield(config,'max_iterations')
    max_iterations = 99;
  end
end


if isfield(spec,'time')
  time = spec.time;
else
  time = nan;
end

% Initialize
if exist('init','var') & isstruct(init)  
  gs_initialize
end

if isempty(spec) & nargout == 0
  return
end

plotit = 2;

if isfield(spec,'plot_settings')
  plot_settings = spec.plot_settings; % Overwrite result of gs_configure
end

if isfield(spec,'fig') & ~isempty(spec.fig)
  if isa(spec.fig,'matlab.ui.Figure')
    fig = spec.fig.Number;
  else
    fig = spec.fig(1);
  end
  if isnan(fig(1))
    fig = 0;
  end;
  fig = round(fig(1));
else
  fig = 0;
end

if isfield(spec,'showgrid') & ~isempty(spec.showgrid)
  showgrid = spec.showgrid;
else
  showgrid = 0;
end

if isfield(spec,'max_iterations') & ~isnan(spec.max_iterations)
  max_iterations = spec.max_iterations;
end

if isfield(spec,'targets');
  targets = spec.targets;
else
  targets = [];
end
if isfield(spec,'weights');
  weights = spec.weights;
else
  weights = [];
end
if isfield(spec,'limits');
  limits = spec.limits;
else
  limits = [];
end
if isfield(spec,'locks');
  locks = spec.locks;
else
  locks = [];
end

if ~isfield(locks,'iv')
  if isfield(targets,'iv')
    locks.iv = nan(nv,1);
  else
    locks.iv = zeros(nv,1); % By default vessel currents are locked to 0
  end
end

if isfield(limits,'rx') & isfield(limits,'zx')
  % Check for all kinds of user errors, only warn if discovered
  n1 = min(size(limits.rx,1),size(limits.zx,1));
  if size(limits.rx,1) ~= size(limits.zx,1) & n1 > 0
    disp('Warning gsdesign: size(limits.rx,1) ~= size(limits.zx,1)')
    disp(['gsdesign will use the first ' num2str(n1) ' rows.'])
  end
  limits.rx = limits.rx(1:n1,:);
  limits.zx = limits.zx(1:n1,:);
  if size(limits.rx,2) == 1
    limits.rx(:,2) = inf;
  end
  if size(limits.zx,2) == 1
    limits.zx(:,2) = inf;
  end
  n2 = min(size(limits.rx,2),size(limits.zx,2));
  if size(limits.rx,2) ~= size(limits.zx,2) & n2 > 0
    disp('Warning gsdesign: size(limits.rx,2) ~= size(limits.zx,2)')
    disp(['gsdesign will use the first ' num2str(n2) ' columns.'])
  end
  limits.rx = limits.rx(:,1:n2);
  limits.zx = limits.zx(:,1:n2);
  k = ~all(isnan(limits.rx'+limits.zx') | isinf(limits.rx'+limits.zx'));
  limits.rx = limits.rx(k,:);
  limits.zx = limits.zx(k,:);
  n1 = size(limits.rx,1);
  % Interpret the numbers
  limits.zx(limits.rx <= 0) = nan;
  limits.zx(isnan(limits.rx)) = nan;
  limits.rx(isnan(limits.zx)) = nan;
  limits.rx(1:n1,n2+1:n2+4) = nan;
  limits.zx(1:n1,n2+1:n2+4) = nan;
  limits.zx(limits.zx < zg(1)+dz) = zg(1)+dz;
  limits.zx(limits.zx > zg(nz)-dz) = zg(nz)-dz;
  limits.rx(limits.rx > rg(nr)-dr) = rg(nr)-dr;
  % Convert all rows of limits.rx, limits.zx to polygons
  dum = min(dr,dz);
  for i = 1:n1
    k = find(~isnan(limits.rx(i,:)));
    n = length(k);
    limits.rx(i,1:n) = limits.rx(i,k);
    limits.rx(i,n+1:end) = nan;
    limits.zx(i,1:n) = limits.zx(i,k);
    limits.zx(i,n+1:end) = nan;
    if n == 2
      rmi = min(limits.rx(i,:));
      rma = max(limits.rx(i,:));
      zmi = min(limits.zx(i,:));
      zma = max(limits.zx(i,:));
      limits.rx(i,:) = nan;
      limits.zx(i,:) = nan;
      limits.rx(i,[1 2 5]) = rmi;
      limits.rx(i,[3 4  ]) = rma;
      limits.zx(i,[2 3  ]) = zmi;
      limits.zx(i,[1 4 5]) = zma;
      n = 5;
    end
    if n > 0
      if limits.rx(i,1) ~= limits.rx(i,n) | limits.zx(i,1) ~= limits.zx(i,n)
	n = n+1;
	limits.rx(i,n) = limits.rx(i,1);
	limits.zx(i,n) = limits.zx(i,1);
      end
      % Create an additional polygon with more coordinates
      limits.rxx(i,1) = limits.rx(i,1);
      limits.zxx(i,1) = limits.zx(i,1);
      limits.cxx(i,1) = 1; % Flag corner coordinates
      k = 1;
      for j = 1:n-1
	r = limits.rx(i,j+1)-limits.rx(i,j);
	z = limits.zx(i,j+1)-limits.zx(i,j);
	m = 1+ceil(sqrt(r^2+z^2)/dum);
	limits.rxx(i,k+(1:m)) = limits.rx(i,j)+r*(1:m)/m;
	limits.zxx(i,k+(1:m)) = limits.zx(i,j)+z*(1:m)/m;
	k = k+m;
	limits.cxx(i,k) = 1;
      end
      limits.nxx(i) = k;
    else
      limits.rxx(i,:) = nan;
      limits.zxx(i,:) = nan;
      limits.nxx(i) = 0;
    end
  end
end
if isfield(spec,'findallxpoints')
  findallxpoints = spec.findallxpoints;
else
  findallxpoints = 0;
end
if isfield(limits,'rx') & isfield(limits,'zx') & ~isempty(limits.rx)
  findallxpoints = 1;
  showallxpoints = 1;
else
  showallxpoints = 0;
end
if isfield(spec,'showallxpoints')
  showallxpoints = spec.showallxpoints;
end
if showallxpoints
  findallxpoints = 1;
end

if ~isfield(limits,'ic')
  if strcmp(upper(tokamak),'D3D')
     % Ampere limits for E&F coils
    limits.ic = [2e3 2e3 ...
      59 59 59 59 54 124 124 59 89 ...
      59 59 59 59 54 124 124 59 89]'*[-100 100];
  end
  if strcmp(upper(tokamak),'NSTXU')
    % Limits from script by Patrick J. Vail
    limits.ic = [-20 0 0 -8 0 0 -13 -13 -13 ...
    -13 0 0 0 0 0 0 -24 -24 -24 -24 -13 -13 -13 -13 0 0 -8 0 0]'*1000;
    limits.ic(:,2) = [20 15 0 13.5 15 15 8 8 8 8 ...
    13 13 13 13 13 13 0 0 0 0 8 8 8 8 15 15 13.5 0 15]'*1000;
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
% Check user-input of limits
if isfield(limits,'ic')
  [s1,s2] = size(limits.ic);
  if s2 > 2 % Something is wrong with the dimensions
    disp('Warning gsdesign: size(limits.ic,2) should be <= 2')
    if s1 < 3 % The dimensions seem mixed up
      disp('limits.ic appears to be transposed from the right format')
      limits.ic = limits.ic';
      disp(['limits.ic was transposed to size ' num2str(s2) 'x' num2str(s1) '!'])
    else
      disp('Remaing columns will be ignored!')
    end
  end
  [s1,s2] = size(limits.ic);
  if s2 >= 2
    if any(limits.ic(:,1) > limits.ic(:,2))
      disp('Warning gsdesign: Lower limits.ic(:,1) higher than upper limits.ic(:,2)')
      limits.ic(:,1:2) = sort(limits.ic(:,1:2),2);
      disp('The two columns of limits.ic were sorted to fix this!')
    end
  end
  if s2 == 1
    limits.ic = abs(limits.ic)*[-1 1];
  end
end

% icci makes TokSys terminal currents ic from circuit currents ci
if isfield(spec,'ciic')
  if size(spec.ciic,2) > nc
    disp(['Warning gsdesign: Projection matrix ciic has more coils ' ...
      'than available and will be cut after coil ' num2str(nc)])
  end
  if isfield(spec,'icci')
    disp(['Warning gsdesign: Circuits were specified with ciic ' ...
      'and icci, ciic will be used, icci will be ignored'])
  end
  if isfield(spec,'cccirc')
    disp(['Warning gsdesign: Circuits were specified with ciic ' ...
      'and cccirc, ciic will be used, cccirc will be ignored'])
  end
  if isfield(spec,'buscode')
    disp(['Warning gsdesign: Circuits were specified with ciic ' ...
      'and buscode, ciic will be used, buscode will be ignored'])
  end
  nci = size(spec.ciic,1);
  n = min(nc,size(spec.ciic,2));
  icci = [pinv(spec.ciic(:,1:n)); zeros(nc-n,nci)];
elseif isfield(spec,'icci')
  if size(spec.icci,1) > nc
    disp(['Warning gsdesign: Projection matrix icci has more coils ' ...
      'than available and will be cut after coil ' num2str(nc)])
  end
  if isfield(spec,'cccirc')
    disp(['Warning gsdesign: Circuits were specified with icci ' ...
      'and cccirc, icci will be used, cccirc will be ignored'])
  end
  if isfield(spec,'buscode')
    disp(['Warning gsdesign: Circuits were specified with icci ' ...
      'and buscode, icci will be used, buscode will be ignored'])
  end
  nci = size(spec.icci,2);
  n = min(nc,size(spec.icci,1));
  icci = [spec.icci(1:n,:); zeros(nc-n,nci)];
elseif isfield(spec,'cccirc')
  nci = max(abs(spec.cccirc));
  icci = zeros(nc,nci);
  for j = 1:nc
    k = abs(spec.cccirc(j));
    s = sign(spec.cccirc(j));
    icci(j,k) = s;
  end
  if isfield(spec,'buscode')
    disp(['Warning gsdesign: Circuits were specified with both ' ...
      'cccirc and buscode, cccirc will be used, buscode will be ignored'])
  end
elseif isfield(spec,'buscode')
  disp('warning in gsdesign: current limits for first coil on bus is ignored')
  icci = eye(nc);
  k = find(spec.buscode);
  icci(min(k),k) = -1;
  icci = icci(:,[1:k(1)-1 k(1)+1:nc]);
  % Add as an extra constraint for D3D that ECOILA = ECOILB
  if strcmp(upper(tokamak),'D3D')
    icci = icci(:,2:end);
    icci(1,1) = 1;
    disp('Since this is D3D, ECOILA and ECOILB were put into 1 circuit')
  end
else
  icci = eye(nc);
end
nci = size(icci,2);

picci = pinv(icci);

limits.ci = ones(nci,1)*[-inf inf];
if isfield(limits,'ic')
  for i=1:nc
    for j = 1:nci
      if icci(i,j) > 0
	if icci(i,j)*limits.ci(j,1) < limits.ic(i,1)
          limits.ci(j,1) = limits.ic(i,1)/icci(i,j);
	end
	if icci(i,j)*limits.ci(j,2) > limits.ic(i,2)
          limits.ci(j,2) = limits.ic(i,2)/icci(i,j);
	end
      elseif icci(i,j) < 0
	if icci(i,j)*limits.ci(j,2) < limits.ic(i,1)
          limits.ci(j,2) = limits.ic(i,1)/icci(i,j);
	end
	if icci(i,j)*limits.ci(j,1) > limits.ic(i,2)
          limits.ci(j,1) = limits.ic(i,2)/icci(i,j);
	end
      end
    end
  end
else
  disp('Warning gsdesign: No limits have been set for coil currents.')
end

if isfield(locks,'ci')
  if isfield(locks,'ic')
    disp(['Warning gsdesign: Locks were specified for both circuits (ci)' ...
    ' and individual coils (ic). Only locks.ci will be honored.'])    
  end
  locks.ci = [locks.ci(:); nan(nci-length(locks.ci(:)),1)];
  locks.ci = locks.ci(1:nci,1); % Just in case user input had too many
else
  locks.ci = nan(nci,1);
  if isfield(locks,'ic')
    locks.ic(length(locks.ic)+1:nc) = nan; % Make sure it's filled
    % Calculate locks.ci
    for i = 1:nc      
      for j = 1:nci
	if picci(j,i)
	  if isnan(locks.ci(j))
            locks.ci(j) = picci(j,i)*locks.ic(i);
	  else
            locks.ci(j) = locks.ci(j) + picci(j,i)*locks.ic(i);
	  end
	end
      end
    end
  end
end

Ipi = nan; % The initial Ip
if isnan(Ipi) & isfield(locks,'cpasma') & any(~isnan(locks.cpasma))
  Ipi = mean(locks.cpasma(~isnan(locks.cpasma)));
end
if isnan(Ipi) & isfield(targets,'cpasma') & any(~isnan(targets.cpasma))
  Ipi = mean(targets.cpasma(~isnan(targets.cpasma)));
end
if isnan(Ipi) & isfield(locks,'rog') & any(~isnan(locks.rog))
  Ipi = mean(locks.rog(~isnan(locks.rog)));
end
if isnan(Ipi) & isfield(targets,'rog') & any(~isnan(targets.rog))
  Ipi = mean(targets.rog(~isnan(targets.rog)));
end

% Initialize automatically if the user didn't supply an init
if max(psizr(:)) == min(psizr(:))
  if isfield(targets,'rsep')
    rsep = targets.rsep;
  else
    rsep = [];
  end
  if isfield(targets,'zsep')
    zsep = targets.zsep;
  else
    zsep = [];
  end
  % Make sure they are the same size
  rsep = rsep+0*zsep;
  zsep = zsep+0*rsep;
  ksep = isinpoly(rsep,zsep);
  rsep = rsep(ksep);
  zsep = zsep(ksep);
  if sum(ksep) > 6
    aminor = (max(rsep)-min(rsep))/2;
    bminor = (max(zsep)-min(zsep))/2;
    rmaxis = min(rsep)+aminor;
    zmaxis = min(zsep)+bminor;
    thsep = angle(rsep-rmaxis+1i*zsep-1i*zmaxis);
    [thsep,ksep] = sort(thsep);
    rsep = rsep(ksep);
    zsep = zsep(ksep);
  else
    aminor = (max(Rlim)-min(Rlim))/3;
    bminor = (max(Zlim)-min(Zlim))/3;
    rmaxis = min(Rlim)+aminor;
    zmaxis = (max(Zlim)+min(Zlim))/2;
    v = linspace(0,2*pi,nbbbs_max)';
    rsep = rmaxis+aminor*cos(v);
    zsep = zmaxis+bminor*sin(v);
  end
  nsep = length(rsep);
  maxis = 1;
  init.pprime = ones(nr,1);
  init.ffprim = ones(nr,1)*mu0;
  init.rg = rg;
  init.zg = zg;
  izr = inpolygon(rgg,zgg,rsep,zsep);
  pcurrt = zeros(nz,nr);
  pcurrt(izr) = RA(izr) + AR(izr);
  if ~isnan(Ipi)
    dum = Ipi/sum(pcurrt(:));
    pcurrt = pcurrt*dum;
    init.pprime = init.pprime*dum;
    init.ffprim = init.ffprim*dum;
  end
  init.psizr_pla = reshape(dpsizrpladpcurrt*pcurrt(:),nz,nr);
  dpsisepdic = zeros(nsep,nc);
  for i = 1:nc
    dpsisepdic(:,i) = gs_interp2(rg,zg,mpc(:,i),rsep,zsep);
  end
  dpsisepdci = dpsisepdic*icci;
  ci = locks.ci;
  iopen = isnan(ci);
  psiseppla = gs_interp2(rg,zg,init.psizr_pla,rsep(:),zsep(:));
  H = 2*dpsisepdci(:,iopen)'*dpsisepdci(:,iopen); % Hessian for devdpv
  h = (2*(psiseppla+dpsisepdci(:,~iopen)*ci(~iopen))'*dpsisepdci(:,iopen))';             % Hessian r.h.s.
  n = sum(iopen);
  Ai = [-eye(n); eye(n)]; % Inequility constraints matrix
  bi = [-limits.ci(iopen,1); limits.ci(iopen,2)]; % Inequility r.h.s.
  ii = ~isinf(bi);
  if any(iopen)
    if license('checkout','optimization_toolbox')
      opts = optimset(...
	'Algorithm','interior-point-convex',...
	'Display','off',...
	'TolX',1e-12,...
	'TolFun',1e-12);
      xqp1 = quadprog((H+H')/2,h,Ai(ii,:),bi(ii),[],[],[],[],[],opts);
      ci(iopen) = xqp1;
    else
      rep = pdipmqpneq2(H,h,Ai(ii,:),bi(ii),100,1e-9,0.97);
      ci(iopen) = rep.x;
    end
  else
    % All coil currents locked, only way to match is adjust cpasma
    % The equation is psisepapp + x1*psiseppla = x2
    % Or let's say: x1*psiseppla + x2 = psisepapp
    psisepapp = dpsisepdci*ci;
    dum2 = -pinv([psiseppla ones(nsep,1)])*psisepapp;
    init.psizr_pla = dum2(1)*init.psizr_pla;
  end
  init.ic = icci*ci;
  init.psizr_app = reshape(mpc*init.ic,nz,nr);
  init.psizr = init.psizr_app+init.psizr_pla;
  init.psimag = gs_interp2(rg,zg,init.psizr,rmaxis,zmaxis);
  init.psibry = 0;
  init.rmaxis = rmaxis;
  init.zmaxis = zmaxis;
  init.circle = 0;
  gs_initialize
  % Reduce bi if ci would exceed limits due to Ipi/cpasma scaling
  % Appears to improve convergence more often than making it worse
  if ~isnan(Ipi) & Ipi ~= 0
    gs_eq_analysis
    k1 = iopen & limits.ci(:,1) < 0;
    k2 = iopen & limits.ci(:,2) > 0;
    dum = Ipi/cpasma*max([ci(k1)./limits.ci(k1,1); ci(k2)./limits.ci(k2,2)]);
    if dum > 1.05
      bi = bi/dum;
      if license('checkout','optimization_toolbox')
	xqp1 = quadprog((H+H')/2,h,Ai(ii,:),bi(ii),[],[],[],[],[],opts);
	ci(iopen) = xqp1;
      else
	rep = pdipmqpneq2(H,h,Ai(ii,:),bi(ii),100,1e-9,0.97);
	ci(iopen) = rep.x;
      end
      init.ic = icci*ci;
      init.psizr_app = reshape(mpc*init.ic,nz,nr);
      init.psizr = init.psizr_app+init.psizr_pla;
      init.psimag = gs_interp2(rg,zg,init.psizr,rmaxis,zmaxis);
      init.psibry = 0;
      init.rmaxis = rmaxis;
      init.zmaxis = zmaxis;
      init.circle = 0;
      gs_initialize
    end
  end
end

ci = picci*ic; % currents in circuits
ic = icci*ci; % ic currents that obey circuit constraints

if constraints == 0 % The coefficients in sp and sf all vary independently
  dxdxcirc = [icci zeros(nc,nx-nc); zeros(nx-nc,nci) eye(nx-nc)];
elseif constraints == 1 % sp=a*sp0 and sf=b*sf0+c*sg0, where now a,b,c vary
   dxdxcirc = [icci zeros(nc,nv+4); ... % Coil currents
               zeros(nv,nci) eye(nv) zeros(nv,4); ...
	       sp0*[zeros(1,nci+nv) 1 0 0 0]; ...
	       sf0*[zeros(1,nci+nv) 0 1 0 0]+ ...
	       sg0*[zeros(1,nci+nv) 0 0 1 0]; ...
	       zeros(1,nci+nv+3) 1];
end
nxcirc = size(dxdxcirc,2);

if ~isnan(Ipi)
  gs_eq_analysis
  if Ipi == 0
    psizr = psizr_app;
    sp(:) = 0;
    sf(:) = 0;
  else
    if Ipi/cpasma < 0.8 | Ipi/cpasma > 1.2
      psizr = psizr*Ipi/cpasma;
      ci    =    ci*Ipi/cpasma;
      ic = icci*ci;
      iv    =    iv*Ipi/cpasma;
      sp    =    sp*Ipi^2/cpasma^2;
      sf    =    sf*Ipi^2/cpasma^2;
      cpasma = Ipi;
    end
  end
end

% Solve for the equilibrium
gs_eq_analysis
if findallxpoints
  gs_nulls_xlim
end
gs_response

iteration_counter = 0;
extra_iterations = 0;
best.costfun = inf;
done = false;
drawoften = true;
costs = [];

while ~done

  gs_find_design_dx
  
  costfun = norm(evc); % Estimate of present error vector
  iteration_counter = iteration_counter+1;
  costs(iteration_counter) = costfun;

  % Estimate impact of the flux error
  maxflerr = max(abs(...
    (dpsizrdx(iused,nx)-...
    (1-psibarzr(iused))*dpsimagdx(nx)-...
    psibarzr(iused)*dpsibrydx(nx))/...
    (psibry-psimag)));
  
  % Only trust the costfun estimate if the flux error is reasonably small
  if maxflerr > 0.1
    % *If* error is coming down, this iteration is best
    best.iteration = iteration_counter;  
  elseif costfun < best.costfun
    % Error is low and iteration worth saving
    best.iteration = iteration_counter;
    best.costfun = costfun;
    best.psizr = psizr;
    best.ci = ci;
    best.iv = iv;
    best.sp = sp;
    best.sf = sf;
    if iteration_counter == 1
      % The first equilibrium can have the best cost function but 
      % not satisfy locks and limits, therefore never revert to it
      best.costfun = inf; % Ensure logic never reverts to this iteration
    end
  end
  
  % Decide what fraction of the change, dx to perform
  
  dpsizr(:) = dpsizrdx*dx;
  dpsimag = dpsimagdx*dx;
  dpsibry = dpsibrydx*dx;
  dpsibarzr = (dpsizr-(1-psibarzr)*dpsimag-psibarzr*dpsibry)/(psibry-psimag);
  drbbbs = drbbbsdx*dx;
  dzbbbs = dzbbbsdx*dx;
  
  maxdr = max(abs(drbbbs(1:nbbbs)));
  maxdz = max(abs(dzbbbs(1:nbbbs)));
  fr = maxdr/aminor;
  fz = maxdz/bminor;
  g = max(fr, fz);
  maxfl = max(abs(dpsibarzr(iused)));
  
  % No bbbs element ever gets to move more than 40% of cell size
  % No flux element gets to change more than 20% in normalized flux
  f = min([1,  0.2/maxfl,  0.4/g]);
  
  % No bbbs element ever gets to move more than 100% of cell size
  % No flux element gets to change more than 50% in normalized flux
  f = min([1,  0.5/maxfl,  1.0/g]);

  % No bbbs element ever gets to move more than 100% of cell size
  % No flux element gets to change more than 10% in normalized flux
  f = min([1,  0.1/maxfl,  1.0/g]);
  
  if iteration_counter > 25
  %  f = min([1,  0.05/maxfl,  0.5/g]);
  end
  if iteration_counter > 50
  %  f = min([1,  0.025/maxfl,  0.25/g]);
  end
  
  %f = min([1,  1.0/maxfl,  2.0/g]); % Too much
  if ~plasma
    f = 1;
  end
  
  lge.psibry = psibry;
  lge.psizr = psizr;
  lge.ci = ci;
  lge.iv = iv;
  lge.sp = sp;
  lge.sf = sf;
  lge.er = er;
  lge.rbdef = rbdef;
  lge.zbdef = zbdef;
  if bdef.target.exists
    lge.dum = ((bdef.target.r-rbdef)/dr)^2+((bdef.target.z-zbdef)/dz)^2;
  end

  f_low_enough = false;
  while ~f_low_enough
    psizr = lge.psizr + f*dpsizr;
    ci    = lge.ci    + f*dxcirc(1:nci);
    iv    = lge.iv    + f*dx(indiv);
    sp    = lge.sp    + f*dx(indsp);
    sf    = lge.sf    + f*dx(indsf);
    er    = lge.er    + f*dx(inder);
    ic = icci*ci;
    psibry_prediction = lge.psibry + f*dpsibrydx*dx;
    gs_eq_analysis
    if bdef.target.exists
      dum = ((bdef.target.r-rbdef)/dr)^2+((bdef.target.z-zbdef)/dz)^2;
    end
    if abs((psibry-psibry_prediction)/(psibry-psimag)) > 0.1 | ...
      bdef.target.exists & lge.dum < 4 & dum > 4
      f = f/2;
    else
      f_low_enough = true;
    end    
  end
  if findallxpoints
    gs_nulls_xlim
  end

  gs_response
 
  % Plot convergence progress if desired
  if plotit > 1
    nhires = 1;
    figure_message = ['Solving... iteration ', num2str(iteration_counter)];
    if isfield(spec,'time')
      figure_message = [figure_message ', time = ' num2str(spec.time)];
    end
    try
      gsdesign_plot_progress
    catch
      fig = ha.figure;
      close(fig)
      gsdesign_plot_progress
    end
  end
      
  done = ~plasma | iteration_counter >= max_iterations | ...
         iteration_counter >= best.iteration+5+extra_iterations;
  
  if exist('ha','var') & ...
     ishandle(ha.button.stop) & ...
     strcmp(get(ha.button.stop,'State'),'on')
    done = true;
  end
  
  % If changes are down to machine precision then stop immediately
  if psimag ~= psibry
    if max(abs(dpsizr(:)/(psimag-psibry))) < 1e-8
      done = true;
    end
  end
  
end

% Check if this last one happens to be the best or if we should restore
gs_find_design_dx
costfun = norm(evc); % Estimate of present error vector
costs(iteration_counter+1) = costfun;

% Skip restore based on best cost function
if 0 & costfun > best.costfun
  costfun = best.costfun;
  psizr = best.psizr;
  ci = best.ci;
  iv = best.iv;
  sp = best.sp;
  sf = best.sf;
  ic = icci*ci;
  gs_eq_analysis
  if findallxpoints
    gs_nulls_xlim
  end
  gs_find_design_dx
end

if plotit > 0
  if plasma
    nhires = 10;
  end
  figure_message = ['Design took ', num2str(iteration_counter) ' iterations'];
  if isfield(spec,'time')
    figure_message = [figure_message ', time = ' num2str(spec.time)];
  end
  try
    gsdesign_plot_progress
  catch
    fig = ha.figure;%........................
% Define PF coil currents

spec.locks.ic = nan*ones(29,1);
spec.targets.ic = nan*ones(29,1);

%coils to be locked at 0
locks_index = [3 4 11 12 13 14 15 16 27 28];

for coil = 1:length(spec.locks.ic)
    if ismember(coil, locks_index)
        spec.locks.ic(coil) = 0;
    else
        spec.targets.ic(coil) = equil_I.cc0t(coil);
    end
end

%current modifications, e.g.
%spec.targets.ic(15) = spec.targets.ic(15) + 5000;
%spec.targets.ic(16) = spec.targets.ic(16) + 5000;
    close(fig)
    gsdesign_plot_progress
  end
end

gs_trace_contours
gs_contour_profiles
V = Vc;
A = Ac;
L = Lc;
T = Tc;
W = Wc;
I = Ic;
B = Bc;
if ~exist('jtavc','var') % TEMP
jtavc = 2*pi/(psibry-psimag)/2e-6*...
  (spline(psibar,Ic,psibar+1e-6)-spline(psibar,Ic,psibar-1e-6))./Acp;
jtavc(1) = rmaxis*pprime(1)+ffprim(1)/mu0/rmaxis;
end
jtav = jtavc;
qpsi = qc;
gs_contour_response
gs_profile_response


% Turn on archiving of spec fields and ev in gs_create_struct_eq
archive_design_data = true;

gs_create_struct_eq

% Avoid false archiving, which might be caused by unconventional use
archive_design_data = false;
