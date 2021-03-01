function [y, eq] = gseq(xc, init, config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   [y, eq] = gseq(xc, init, config)
%
%  PURPOSE: Solve the 2-D (free-boundary) Grad-Shafranov equation
%           CONFIGURE & INITIALIZE: xc0 = gseq([], init, config)
%           SOLVE for equilibrium: [y, eq] = gseq(xc)
%
%  INPUTS:     xc, equilibrium states
%            init, initial equilibrium, for further instructions type:
%                  edit gseq_init.m
%          config, tokamak description and options, to set these type:
%                  edit gseq_config.m
%
%  OUTPUTS: y, outputs defined in config.outputs, 
%              gseq('index_in_y') returns indices for named outputs
%          eq, equilibrium in TokSys format
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	8/14/15
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Recall variables needed for fast loop
persistent plasma cpasma lae evolve_option constraints Pxx dxdxc converge
persistent stabilize stabilized halo circle_model plasma_is_tiny lim
persistent ic iv sp sf er psizr psibarzr psibry psimag dpsimagdx dpsibrydx
persistent dpsizrdx dpsizrx dpsizr dpsibar_limit_for_linear err phaseouterrors
persistent dpsihbardx dr0dx dz0dx da0dx dlidx dbetapdx dgapdpdx dzcurdx gamma
persistent wa wb iia iib iplasma maxis drmaxisdx dzmaxisdx dcpasmadx
persistent indic indiv indis indsp indsf inder index_in_y
persistent dydx dysdx dpdx ulim vlim daminordx drsurfdx
persistent Amat Bmat sb0 vresc
persistent nc nci nv nkn nx nxc nxs nps ncd nrh np a0 rdp zdp ndp
persistent rtplot_handles plot_settings default_time titstrlen title_time_format
persistent nbbbs_max nbbbs rbbbs zbbbs drbdx dzbdx rbdef zbdef drbdefdx dzbdefdx
persistent ecdata fcdata vvdata nl rl zl shape nr nz ngg rg zg dr dz rgg zgg 
persistent plot_times plot_counter plot_if_new_response dtplot time_for_next_plot
persistent xpix zpix nxpix nypix fpix Pcg ilimgg
persistent V W A I T jtav rhot qpsi calculate_helical_voltage dvindcdxdot
persistent calculate_profiles calculate_profile_responses drhotdx djtavdx

if isfield(xc,'time')
  time = xc.time;
else
  time = default_time;
end

% FAST LOOP
if nargin == 1 & ~ischar(xc) & ~converge    
  gs_update_fast % Use linear response
  if equilibrium_update_is_complete
    % PLOT
    if (plot_counter <= length(plot_times) && ...
       time >= plot_times(plot_counter)) || ...
       (dtplot > 0 && time >= time_for_next_plot) || ...
       (isfield(xc,'rtplotit') && xc.rtplotit)
      gs_rtplot_prepare
      gs_rtplot
      if dtplot > 0 & time >= time_for_next_plot  
        time_for_next_plot = (1+round(time/dtplot))*dtplot-dtplot/1e3;
      end
    end    
    if nargout == 0
      clear y
    end
    return  
  end
else
  equilibrium_update_is_complete = false;
end

% How to calculate outputs
persistent outputs iy ny user_signal
persistent fl dfldx nfl mlc mlv mpl
persistent lv dlvdx nlv mhc mhv mph
persistent bp dbpdx nbp gbc gbv gpb 
persistent rog drogdx nrog rldata

% Extra diagnostics
persistent mdc mdv mpd grdc grdv grpd gzdc gzdv gzpd
persistent psipd brdp bzdp btdp dpsipddx dbrdpdx dbzdpdx dbtdpdx
persistent nxpoints kxpoints xpointorder

% The tokamak
persistent tokamak mpc mpv mpp mcc mcv mvv ress rss mss
persistent nlim Rlim Zlim drl dzl dl wl wld wlb wlt wlr wlz iil
persistent piccc klimgg rvmin rvmax zvmin zvmax concavel turnin

% Options
persistent psikn plotit default_plotit use_circle_model
persistent default_converge min_iterations max_iterations
persistent ncp2gr ngr2cp

% Remember helping variables
persistent mu0 twopi R13 mx neighbors ir4 iz4 ir8 iz8 ns
persistent sp0 sf0 sg0 c0 c1 c2 c3 d0 d1 d2 d3 d4 M kk AR Ag RA 
persistent iknotg psibar a b c d aa bb q1 q2 x1 x2 x3
persistent id4 id8 dpsid4 dpsid8 dpsib4 dpsib8 dpsib16 dCldrbbbs dCldzbbbs
persistent fer dzbbbs_max plasma_size ibbbs gbbbs xx wbbbs wds wbs wrs wzs
persistent rhobbbs dpsibbbsdr dpsibbbsdz dpsizrpladpcurrt
persistent drb3dpsia drb3dpsib drb3dpsip dzb3dpsia dzb3dpsib dzb3dpsip
persistent redge zedge fedge gedge iedge xedge thedge rhoedge

% For open-field-line plasma model
persistent geohzr ihzr ihalo 
persistent smallest_aminor open_field_line_volume open_field_line_rmaxis

% Last analyzed equilibrium
persistent rmaxis zmaxis itl iitl wtl
persistent rzero bzero preszr pprimezr ffprimzr
persistent iused zbot ztop psipla li betap betan pcurrt 
persistent fpol ffprim pres pprime psix1 psix2 thbbbs Wth
persistent psizr_err Bp2zr Cl Bp2V Vtot bp2flx fluxerror
persistent p2 p3 f2 f3 r1 r2 r3 z1 z2 z3 rcell zcell Acell ncell

% Plasma response
persistent dpsipladx dAcelldx dpcurrtdx dWthdx lstar
persistent dpsibrydp r0 z0 psih
persistent drbbbsdx dzbbbsdx dpcoredpsi dpsizrdpsizrapp dpsizr_err
persistent drcurdx dbetandx
persistent response_count

% Variables that are used to construct sp0, sf0, sg0
persistent psibry0 psimag0 pprime0 ffprim0 no_edge_current no_edge_gradient

% Parameters used by gs_dynamics to evolve x
persistent plares rxx Rpla Rpla_previous drcurdv dzcurdv
persistent dxdx3 eta_vs_rhot

% Circuits
persistent icci picci ci Rhat Vhat Pxxi Rextra Lextra

% For EFITs and cc_efit_to_tok
persistent ecid fcid vvid fcturn ecturn turnfc fcnturn ecnturn vvfrac
persistent idx_efit_to_tok imks iterminal

% Read or write to persistent variables from caller's ws
if ischar(xc)
  if nargin == 1
    if exist(xc,'var')
      y = eval(xc);
    elseif strcmp(xc,'eq') % Return a structure eq
      gs_eq_analysis
      gs_create_struct_eq
      y = eval(xc);
    elseif strcmp(upper(xc),'FUN2WS') % All persistent to caller's ws
      vars = who;
      for i = 1:size(vars,1)
	assignin('caller',char(vars(i)),eval(char(vars(i))))
      end
    elseif strcmp(upper(xc),'WS2FUN') % From caller's ws to persistent
      vars = who;
      for i = 1:size(vars,1)
	if evalin('caller',['exist(''' char(vars(i)) ''',''var'')'])
	  eval([char(vars(i)) ' = evalin(''caller'',''' char(vars(i)) ''');'])
	end
      end
    else
      y = [xc ' is not a persistent variable in gseq'];
    end
    return
  elseif nargin == 2
    if exist(xc,'var')
      eval([xc ' = init;']);
      y = init;
    elseif strcmp(upper(xc),'HELP')
      gshelp(u,'gseq')
    else
      disp([xc ' is not a persistent variable in gseq']);
      y = [xc ' is not a persistent variable in gseq'];
    end
    if nargout == 0
      clear y
    end
    return
  end
end

% CONFIGURE (process input config)
if exist('config','var')
  gs_configure
end

new_response_was_calculated = false;

% INITIALIZE (process input init)
if exist('init','var')
  gs_initialize
  if circle_model
    gscp_analysis_response
  else
    gs_eq_analysis
    gs_response
  end
end

% UPDATING EQUILIBRIUM (process input xc to make output y)
if converge
  gs_converge % Find precise solution
else
  gs_update_fast % Use linear response
  if ~equilibrium_update_is_complete    
    gs_update_slow % Update linear response
  end
end

if isempty(xc) % Return xc0
  xc0.ic = ic;
  xc0.iv = iv;
  if constraints == 0
    xc0.sp = sp;
    xc0.sf = sf;
  elseif constraints == 1
    xc0.ip = cpasma;
    xc0.li = li;
    xc0.betap = betap;
  elseif constraints == 2
    xc0.ip = cpasma;
    xc0.li = li;
    xc0.Wth = Wth;
  end
  y = xc0;
end

% PLOT
if ~isempty(xc)
  if new_response_was_calculated & plot_if_new_response
    gs_rtplot
  elseif (plot_counter <= length(plot_times) && ...
         time >= plot_times(plot_counter)) || ...
         (dtplot > 0 && time >= time_for_next_plot) || ...
         (isfield(xc,'rtplotit') && xc.rtplotit)
    gs_rtplot_prepare
    gs_rtplot    
    if dtplot > 0 & time >= time_for_next_plot
      time_for_next_plot = (1+round(time/dtplot))*dtplot-dtplot/1e3;
    end
  end
end

if nargout == 0
  clear y
end

if nargout > 1
  if ~exist('psizr_app','var')
    gs_eq_analysis
  end
  gs_create_struct_eq
end
