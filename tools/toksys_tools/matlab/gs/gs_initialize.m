%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_initialize
%
%  PURPOSE: Initialize gseq, load an equilibrium
%
%  INPUTS: init, structure with initial equilibrium quantities.
%                These are used if available:
%                rg, zg, psizr, psimag, rmaxis, zmaxis, psibry, 
%                pprime, ffprim, rzero, bzero
%                The grid in init does not need to match the grid in config
%
%  OUTPUTS: Initialized equilibrium
%
%  RESTRICTIONS: The initial equilibrium in gs will only approximate init
%                if any of the following occurs:
%                1. the grids in init and config don't match
%                2. spline knots for pres, fpol in init don't match gs spec
%                3. edge currents exist (gs considers partial cell coverage)
%                4. init isn't converged
%                Normally none of these will produce a serious discrepancy
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  METHOD: 
	
%  NOTES:  
	
%
%  WRITTEN BY:  Anders Welander  ON	4/2/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For parsing efits with cc_efit_to_tok and find what the coil currents are
vac_objs.imks = imks;
vac_objs.iterminal = iterminal;
vac_objs.idx_efit_to_tok = idx_efit_to_tok;
vac_objs.mcc = mcc;
vac_objs.ecnturn = ecnturn;
vac_objs.fcnturn = fcnturn;
vac_objs.nv = nv;
if isfield(init,'ecid')
  vac_objs.def_connect.ecid = init.ecid;
else
  vac_objs.def_connect.ecid = [];
end
if isfield(init,'fcid')
  vac_objs.def_connect.fcid = init.fcid;
else
  vac_objs.def_connect.fcid = [];
end
if isfield(init,'vvid')
  vac_objs.def_connect.vvid = init.vvid;
else
  vac_objs.def_connect.vvid = [];
end
if isfield(init,'ecturn')
  vac_objs.def_connect.ecturn = init.ecturn;
else
  vac_objs.def_connect.ecturn = [];
end
if isfield(init,'fcturn')
  vac_objs.def_connect.fcturn = init.fcturn;
end
if isfield(init,'turnfc')
  vac_objs.def_connect.turnfc = init.turnfc;
end
if isfield(init,'vvfrac')
  vac_objs.def_connect.vvfrac = init.vvfrac;
else
  vac_objs.def_connect.vvfrac = [];
end
% Done creating vac_objs


% Now reading new initial values for the equilibrium states

if isfield(init,'ic')
  ic = init.ic(:);
elseif isfield(init,'cc')
  equil_I = cc_efit_to_tok(vac_objs,init);
  ic = equil_I.cc0t(:);
end
if ~exist('ic','var') | isempty(ic)
  ic = zeros(nc,1);
end
if length(ic) ~= nc
  error('The length(init.ic) must equal size(config.mpc,2)')
end

if isfield(init,'iv')
  iv = init.iv(:);
elseif isfield(init,'vc')
  equil_I = cc_efit_to_tok(vac_objs,init);
  iv = equil_I.vc0t(:);
end
if ~exist('iv','var') | isempty(iv)
  iv = zeros(nv,1);
end
if length(iv) ~= nv
  error('The length(init.iv) must equal size(config.mpv,2)')
end

% psimag and psibry are needed if profiles are given by pprime, ffprim
if isfield(init,'psibry')
  psibry = init.psibry;
end
if ~exist('psibry','var') | isempty(psibry)
  psibry = 0;
end
if isfield(init,'psimag')
  psimag = init.psimag;
end
if ~exist('psimag','var') | isempty(psimag)
  psimag = psibry+1;
end

% rzero and bzero are needed
if isfield(init,'rzero')
  rzero = init.rzero;
end
if ~exist('rzero','var') | isempty(rzero)
  if strcmp(upper(tokamak),'ITER')
    rzero = 6.2;
  elseif strcmp(upper(tokamak),'EAST')
    rzero = 1.965;
  elseif strcmp(upper(tokamak),'D3D')
    rzero = 1.6955;
  else
    rzero = (max(rl)+min(rl))/2;
  end
end
if isfield(init,'bzero')
  bzero = init.bzero;
end
if ~exist('bzero','var') | isempty(bzero)
  if strcmp(upper(tokamak),'ITER')
    bzero = -5.3;
  elseif strcmp(upper(tokamak),'EAST')
    bzero = 2;
  elseif strcmp(upper(tokamak),'D3D')
    bzero = -2;
  else
    bzero = rzero;
  end
end

% calculate coefficients for template polynomials in all nkn+2 intervals
p0 = c0*sp0;
p1 = c1*sp0;
p2 = c2*sp0;
p3 = c3*sp0;
f0 = c0*sf0;
f1 = c1*sf0;
f2 = c2*sf0;
f3 = c3*sf0;
g0 = c0*sg0;
g1 = c1*sg0;
g2 = c2*sg0;
g3 = c3*sg0;
pbnkn = psibar*ones(1,nkn+2);
% Sometimes init contains incorrect values, warn and correct
if ~isfield(init,'sp')
  if isfield(init,'pres') & all(init.pres==0) & ...
     isfield(init,'pprime') & ~all(init.pprime==0)
    disp('Warning gs_initialize: init.pres is all zeros but init.pprime isn''t')
    disp('                       init.pprime will be used for initialization')
    init.pres(:) = nan;
  end
end
if isfield(init,'sp') & length(init.sp(:)) == nkn+2
  sp = init.sp(:);
elseif isfield(init,'pres') & ~any(isnan(init.pres))
  pres = spline(linspace(0,1,length(init.pres))',init.pres,psibar);
  pp0 = (p0(iknotg) + p1(iknotg).*psibar + ...
    p2(iknotg).*psibar.^2 + p3(iknotg).*psibar.^3);
  if constraints
    sp = sp0*pinv(pp0)*pres;
  else
    sp = pinv(c0(iknotg,:) + c1(iknotg,:).*pbnkn + ...
      c2(iknotg,:).*pbnkn.^2 + c3(iknotg,:).*pbnkn.^3)*pres;
  end
elseif isfield(init,'pprime')
  pprime = spline(linspace(0,1,length(init.pprime))',init.pprime,psibar);
  pp0 = (p1(iknotg) + p2(iknotg)*2.*psibar + ...
    p3(iknotg)*3.*psibar.^2)*twopi/(psibry-psimag);
  if constraints
    sp = sp0*pinv(pp0)*pprime;
  else
    sp = pinv(c1(iknotg,:) + c2(iknotg,:)*2.*pbnkn + ...
      c3(iknotg,:)*3.*pbnkn.^2)*pprime/twopi*(psibry-psimag);
  end
end
if isfield(init,'sf') & length(init.sf(:)) == nkn+2
  sf = init.sf(:);
elseif isfield(init,'fpol')
  fpol = spline(linspace(0,1,length(init.fpol))',init.fpol,psibar);
  fvac = init.fpol(end);
  bzero = fvac/rzero;
  ff0 = (f0(iknotg) + f1(iknotg).*psibar + f2(iknotg).*psibar.^2 + ...
	  f3(iknotg).*psibar.^3);
  gg0 = (g0(iknotg) + g1(iknotg).*psibar + g2(iknotg).*psibar.^2 + ...
	  g3(iknotg).*psibar.^3);
  if constraints
    sf = [sf0 sg0]*pinv([ff0 gg0])*(fpol.^2/2-fvac^2/2);
  else
    sf = pinv(c0(iknotg,:) + c1(iknotg,:).*pbnkn + ...
      c2(iknotg,:).*pbnkn.^2 + c3(iknotg,:).*pbnkn.^3)*(fpol.^2/2-fvac^2/2);
  end
elseif isfield(init,'ffprim')
  ffprim = spline(linspace(0,1,length(init.ffprim))',init.ffprim,psibar);
  ff0 = (f1(iknotg) + f2(iknotg)*2.*psibar + ...
	  f3(iknotg)*3.*psibar.^2)*twopi/(psibry-psimag);
  gg0 = (g1(iknotg) + g2(iknotg)*2.*psibar + ...
	  g3(iknotg)*3.*psibar.^2)*twopi/(psibry-psimag);
  if constraints
    sf = [sf0 sg0]*pinv([ff0 gg0])*ffprim;
  else
    sf = pinv(c1(iknotg,:) + c2(iknotg,:)*2.*pbnkn + ...
      c3(iknotg,:)*3.*pbnkn.^2)*ffprim/twopi*(psibry-psimag);
  end
end
plasma = any(sp) | any(sf);

% Allow grids in init and config to be different
if isfield(init,'psizr') & ...
  isfield(init,'rg') & isfield(init,'zg')
  psizr = gs_interp2(init.rg,init.zg,init.psizr,rgg,zgg);
end

psizr_app = reshape(mpc*ic + mpv*iv, nz, nr);
if ~plasma % No plasma means we know psizr
  psizr(:) = psizr_app;
end

if isfield(init,'circle')
  circle = init.circle;
end
if ~exist('circle','var') | isempty(circle)
  circle = 0;
end

if isfield(init,'cpasma')
  cpasma = init.cpasma;
end
if isfield(init,'rmaxis')
  rmaxis = init.rmaxis;
end
if isfield(init,'zmaxis')
  zmaxis = init.zmaxis;
end

if ~exist('er','var') | isempty(er)
  er = 1;
end

% Initialize errors that are morphed away in gs_update_eq
err.ys     = zeros(nc+nv,1);
err.li     = 0;
err.betap  = 0;
err.cpasma = 0;

% Plotting
rtplot_handles.refresh = true;
plot_settings.refresh = true;
% For time-dependent & conditional plotting
% Things that change during the simulation can be initialized here
plot_counter = 1;
time_for_next_plot = -inf;
titstrlen = 0;

% EFIT stuff
if isfield(init,'fcturn')
  fcturn = init.fcturn;
end
if isfield(init,'turnfc')
  turnfc = init.turnfc;
end
if isfield(init,'fcid')
  fcid = init.fcid;
end
if isfield(init,'ecturn')
  ecturn = init.ecturn;
end
if isfield(init,'ecid')
  ecid = init.ecid;
end
if isfield(init,'idx_efit_to_tok')
  idx_efit_to_tok = init.idx_efit_to_tok;
end

idiot = init;

% iccc makes a TokSys ic from an efit cc
if isfield(idiot,'cc')
  ncc = length(init.cc);
else
  ncc = nc;
end
iccc = eye(nc,ncc);
try
  for j = 1:ncc
    idiot.cc = zeros(ncc,1);
    idiot.cc(j) = 1;
    equil_I = cc_efit_to_tok(vac_objs,idiot);
    iccc(:,j) = equil_I.cc0t;
  end    
end
piccc = pinv(iccc);

% ivvc makes a TokSys iv from an efit vc
if isfield(idiot,'vc')
  nvc = length(init.vc);
else
  nvc = nv;
end
ivvc = eye(nv,nvc);
try
  for j = 1:nvc
    idiot.vc(:) = zeros(nvc,1);
    idiot.vc(j) = 1;
    equil_I = cc_efit_to_tok(vac_objs,idiot);
    ivvc(:,j) = equil_I.vc0t;
  end    
end
pivvc = pinv(ivvc);
