function [tok_system,dbg_objs] = build_tokamak_system(build_inputs)
%
%  USAGE: tok_system= build_tokamak_system(build_inputs)
%
%  PURPOSE: Generic script to build tokamak axisymmetric system from 
%	pre-calculated objects. Includes R,Z,Ip response in system.
%
%  INPUTS: (in build_inputs structure)
%   tokamak = name of tokamak system to construct (required)
%   vacuum_objs = structure containing device model objects (required)
%   cccirc = connection vector for CC's 
%   netlist = network connection description (Only one of cccirc or netlist can
%             be input.  If neither is input, default connections are used.)
%   efit_gfile = equilibrium gfile to read (not required if vacuum_model=1)
%   vacuum_model = flag, = 1 get vacuum model, =0 includes plasma (default)
%   ichooseq= equilibrium file type:
%	1 = efit_gfile
%	2 = corsica generated flat files
%	3 = saved corsica equilibrium
%	4 = TokSys equilibrium data structure
%   corsica_inputs = name of directory containing corsica input (flat) files
%			(required only if ichooseq=2)
%   corsica_savefile = name of saved corsica equilibrium 
%			(required only if ichooseq=3)
%   equil_data = name of TokSys format equilibrium data structure 
%			(required only if ichooseq=4)
%   irzresp_dynamic, irzresp_output = flags determining plasma contribution to
%	      dynamic (state) and output equations (default = 3 for both):
%        0 = motionless plasma
%	 1 = rigid R motion
%	 2 = rigid Z
%	 3 = both R and Z rigid motions
%	 4 = nonrigid plasma response based on vst in corsica
%	 5 = gspert, (perturbed grad-shafronv equation plasma response)
%
%  OPTIONAL INPUTS:
%   vvgroup = vacuum vessel grouping vector: length = all vessel elements
%		set element value = k to belong to group k, 0 to not include
%		(default= [1:total # elements])
%   vvcirc  = circuit connection vector for vac. vessel(default=all independent)
%   iplcirc = flag to select inclusion of plasma circuit in state space(1)
%		or not(0) (optional, default = 1 if vacuum_model=0)
%   Rp       = (optional) user override of computed plasma resistance
%   Rext     = extra circuit resistance to add to each coil and vessel element
%		   can be either length=#coils or #coils + #vessel (default=0)
%   Lext     = extra circuit inductance to add to each coil/vessel
%		   can be either length=#coils or #coils + #vessel (default=0)
%   replace_Rext= if 1, instead of adding Rext, replace the computed coil/vessel
%			         resistance with Rext (default=0)
%   cc_file = file containing coil currents (optional, usually given by efit)
%   idx_efit_to_tok = optional map of efit indices to toksys indices, s.t. if
%       I=currents in efit order, then I(idx_efit_to_tok)=currents in toksys order
%   scldzdis = To increase gamma when irzresp_dynamic=irzresp_output=4, 
%			make scldzdis<1.0. (optional, default=1)
%   tok_geo_plot_fn = one of 'generic', 'plot_east_geo', 'plot_kstar_geo', etc.
%		(optional, default = '' or 'none' => no plots)
%   Te_res = plasma electron temp for plasma resistance calc [eV]
%		      (optional, default = 4000eV)
%   li_res = plasma internal inductance for resistance calculation
%		       (optional,default = 0.5)
%   Zeff_res = Zeff for resistance calc (optional, default = 1.5)
%   scale_R_scpf = scale factor to multiply resistance of superconducting
%		   PF coils, scalar or vector of length = number of coils
%	 	             (optional, default = 1 applied to ALL coils)
%   netlist_currents = array of strings defining branches for which to output
%                      currents (optional, default = no current outputs, 
%                      not used if cccirc used for circuit connectivity)
%   netlist_voltages = n x 2 matrix, each row containing node numbers [N1 N2],
%                      with voltage output defined as V(N1)-V(N2) 
%                      (optional, default = no voltage outputs, 
%                       not used if cccirc used for circuit connectivity)
%   gspert_options = options for gspert calculation (only used if irzresp values=5)
%   iwait = (optional) flag: 1=wait when error messages, 0(default)=no wait
%   verbose  = level of screen output. Allowed values are:
%         0 = no messages or uncritical warnings (default)
%         1 = display warnings
%         2 = display warnings and other messages
%
%  OUTPUTS:
%     tok_system = structure containing model of tokamak system
%
%  RESTRICTIONS:
%    	Must have specified objects files in specified path.
%	Dimensions of equilibrium objects much match electromagnetic system
%	dimensions (eg grid rg,zg must match crj dimensions).
%
%  METHOD:  
%       Loads and defines objects needed for control analysis/design.
%	cccirc uses Corsica connection vector convention: vector of indices 
%	identifying which circuit each coil belongs to ("index" is negative to
%	denote antiseries...). So eg [1 2 3 1 -2 -3] means coil 1 is series'd
%	with coil 4, coil 2 is antiseries'd with coil 5, coil 3 is antiseriesd
%	with coil 6, and there are 3 total circuits in the final connected
%	system.

%  WRITTEN BY:  Mike Walker	ON 	2/3/05
%   	based on build_east_sys.m by Dave Humphreys
%
%  MODIFICATIONS: 
%     2006-03-31  Robert Deranian   Added explicit coil current calculation
%     2008-02-07  NWE Added rogs & loop voltages to cmat,dmat,hmat
%     2008-10-02  ASW Added option to use vst response
%     2010-06-07  ASW Added option to use gspert plasma response   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OBSOLETE:
%   Rckt     = extra resistance to add to each coil circuit (default=0)
%   Lckt     = extra inductance to add to each coil circuit (default=0)
%   replace_Rckt= if 1, replace computed coil resistance with Rckt (default=0)
tok_system=[];	% define in case we must return on error
tol_symmetry = 0.01;	% 1 percent error in symmetry of xmatx is OK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse inputs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

struct_to_ws(build_inputs);

if(exist('verbose','var')==0)
   verbose = 0;
end

% Always required:
if(~exist('tokamak'))
   wait('ERROR build_tokamak_system: input "tokamak" is REQUIRED')
   return;
end
if(~exist('vacuum_objs'))
   wait('ERROR build_tokamak_system: input "vacuum_objs" is REQUIRED')
   return;
end
if(exist('netlist'))	% build circuit using netlist input
   if(exist('cccirc'))
      wait('ERROR build_tokamak_system: only one of cccirc or netlist allowed')
      return;
   end
   use_netlist= 1;
   use_cccirc = 0;
   if(~exist('netlist_currents'))
      netlist_currents = [];
   end
   if(~exist('netlist_voltages'))
      netlist_voltages = [];
   end
else 				% else process with cccirc, which has defaults.
   use_netlist= 0;
   use_cccirc = 1;
   if(~exist('cccirc'))
      cccirc = [1:vacuum_objs.nc];
      ncx = vacuum_objs.nc; 
   else
      if(length(cccirc)~=vacuum_objs.nc)
        fprintf('ERROR build_tokamak_system: length(cccirc) must match num. coils= %d\n',vacuum_objs.nc)
        wait;
        return;
      end
      ncx = max(abs(cccirc)); 
   end

% Check for consistency of specified vvgroup and vvcirc
   if(exist('vvgroup')==1) 
      if(isempty(vvgroup))
         temp = [];
      else
         temp = vvgroup(:,1);
      end
      temp = setdiff(unique(temp),0);
      if(length(temp) ~= max(temp))
         wait('ERROR build_tokamak_system: vvgroup must define contiguous numbered groups')
         return;
      end
      if(exist('vvcirc')==1) 
         temp1 = setdiff(union(temp,abs(vvcirc)),0);
         if(length(temp1)~=length(temp))
            wait('ERROR build_tokamak_system: vvcirc and vvgroup inconsistent')
            return;
         end
      end
   end

   resv = vacuum_objs.resv;
   if(~exist('vvgroup')) 	% vvgroup processing only with cccirc
      vvgroup = [1:vacuum_objs.nv]';
      nvx= vacuum_objs.nv;
   elseif(~isempty(vvgroup))
      if(size(vvgroup,2)>2)
        wait('ERROR build_tokamak_system: vvgroup must be 1 or 2 column matrix')
        return;
      end
      if(size(vvgroup,1)~=length(resv))
         wait(['ERROR build_tokamak_system: vvgroup col. size must match '...
		'size of resv = ' int2str(length(resv))])
         return;
      end
      nvx = max(abs(vvgroup(:,1)));	 %Find reduced size of VV
   else				% vvgroup exists but is empty
      vvgroup = [1:vacuum_objs.nv]';
      nvx= vacuum_objs.nv;
   end

% vvgroup now guaranteed to exist for operation on by vvcirc logic.

   if(~exist('vvcirc') | isempty(vvcirc)) 	
      vvcirc = unique(vvgroup(:,1));
      vvcirc = setdiff(vvcirc,0);	% remove any 0 elements
   end

% and now vvcirc is guaranteed to exist for all that follows.

   nvx = length(setdiff(unique(abs(vvcirc)),0));
end
if(~exist('vacuum_model') | isempty(vacuum_model))
   vacuum_model = 0;
end
if(~exist('ichooseq'))
   ichooseq = 1;	% default is efit_gfile
end
if(ichooseq==1 & ~exist('efit_gfile'))
   if(~vacuum_model)
      wait(['ERROR build_tokamak_system:' ...
	' input "efit_gfile" is REQUIRED if vacuum_model~=0'])
      return;
   end
end
if(ichooseq==1 & ~exist('efit_gfile'))
  wait('ERROR build_tokamak_system: efit_gfile must be specified when ichooseq=1')
  return;
end
if(ichooseq==2 & ~exist('corsica_inputs'))
   if(~vacuum_model)
      wait(['ERROR build_tokamak_system:' ...
	' input "corsica_inputs" is REQUIRED if vacuum_model~=0'])
      return;
   end
end
if(ichooseq==4 & ~exist('equil_data'))
   if(~vacuum_model)
      wait(['ERROR build_tokamak_system:' ...
	' input "equil_data" is REQUIRED if vacuum_model~=0'])
      return;
   end
end
if(exist('names')==1) 
   wait('WARNING build_tokamak_system: use of input "names" is obsolete and will be ignored. \n')
end
names = struct('dummy','');
if(~exist('scale_R_scpf') | isempty(scale_R_scpf))
   scale_R_scpf = 1;
end
if(~exist('iplcirc') | isempty(iplcirc))
   if(vacuum_model)
      iplcirc = 0;
   else
      iplcirc = 1;
   end
end
if(~exist('Te_res') | isempty(Te_res))
   Te_res = 4000;
end
if(~exist('li_res') | isempty(li_res))
   li_res = 0.5;
end
if(~exist('Zeff_res') | isempty(Zeff_res))
   Zeff_res = 1.5;
end
if(~exist('tok_geo_plot_fn') | isempty(tok_geo_plot_fn))
   tok_geo_plot_fn = '';
elseif(strncmp(tok_geo_plot_fn,'none',4))
   tok_geo_plot_fn = '';
end
if(~exist('scldzdis') | isempty(scldzdis))
   scldzdis=1;
end
if(~exist('Rext'))
   Rext = 0;
elseif(size(Rext,2)>1)  % make into column vector
   Rext = Rext';
end
if(exist('Rckt'))
  if verbose > 0
    disp('WARNING build_tokamak_system: Rckt is obsolete. Replace with Rext')
  end
else
   Rckt = 0;
end
if(~exist('Lext'))
   Lext = 0;
elseif(size(Lext,2)>1)  % make into column vector
   Lext = Lext';
end
if(exist('Lckt'))
  if verbose > 0
    disp('WARNING build_tokamak_system: Lckt is obsolete. Replace with Lext')
  end
else
   Lckt = 0;
end
if(exist('replace_Rext'))
   if(exist('replace_Rckt'))
      wait(['ERROR build_tokamak_system: ' ...
                'cannot specify both replace_Rext and replace_Rckt'])
      return;
   end
else
   replace_Rext=0;
end
if(~exist('replace_Rckt'))
   replace_Rckt=0;
end
%if(~exist('replace_L'))
%   replace_L=0;
%end
if(~exist('verbose'))
   verbose = 0;
end
if(~exist('cc_file'))
   cc_file = '';
end
if(~exist('IpRog_conds'))
   IpRog_conds=[];
end
if(~exist('irzresp_dynamic'))
   irzresp_dynamic=3;
end
if(~exist('irzresp_output'))
   irzresp_output=3;
end

if(vacuum_model & iplcirc)
  wait('ERROR build_tokamak_system: vacuum_model and iplcirc are INCONSISTENT')
  return;
end
if(~isfield(build_inputs,'iwait'))
   iwait = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data needed to construct models.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prelims and Constants:
mu0 = 0.4*pi;
twopi = 2*pi;

% Tokamak data:
ncc     = vacuum_objs.nc;
nvv     = vacuum_objs.nv;
nfl     = vacuum_objs.nfl;
nbp     = vacuum_objs.nbp;
if(isfield(vacuum_objs,'ecdata'))
   ecdata  = vacuum_objs.ecdata;
else
   ecdata = [];
end
fcdata  = vacuum_objs.fcdata;
vvdata  = vacuum_objs.vvdata;
limdata = vacuum_objs.limdata;
fldata  = vacuum_objs.fldata;
bpdata  = vacuum_objs.bpdata;
ccnturn = vacuum_objs.ccnturn;

mlc     = vacuum_objs.mlc;
mlv     = vacuum_objs.mlv;
mpl     = vacuum_objs.mpl;

if(isfield(vacuum_objs,'rldata')) % overrides Iprog_conds if exists
   rldata = vacuum_objs.rldata;
   nrl = vacuum_objs.nrl;
elseif(~isempty(IpRog_conds))     % backwards compatbility
   rldata = IpRog_conds';
   rldata(end+1,:) = 1;   % add Ip support by default
   nrl = size(rldata,2);
else
   rldata = [];
   nrl = 0;
end


if(isfield(vacuum_objs,'lvdata'))
   lvdata  = vacuum_objs.lvdata;
   nlv     = vacuum_objs.nlv;
   mhc     = vacuum_objs.mhc;
   mhv     = vacuum_objs.mhv;
   mph     = vacuum_objs.mph;
else
   nlv = 0;
   mhc = [];
   mhv = [];
   mph = [];
end

if(isfield(vacuum_objs,'msedata'))
   nmse    = vacuum_objs.nmse;
   msedata = vacuum_objs.msedata;
   gmsebrc = vacuum_objs.gmsebrc;
   gmsebzc = vacuum_objs.gmsebzc;
   gmsebrv = vacuum_objs.gmsebrv;
   gmsebzv = vacuum_objs.gmsebzv;
   gmsebrp = vacuum_objs.gmsebrp;
   gmsebzp = vacuum_objs.gmsebzp;
else
   nmse = 0;
   gmsebrc = [];
   gmsebzc = [];
   gmsebrv = [];
   gmsebzv = [];
end

if(~isempty(vacuum_objs.gpb))
   gbc  = vacuum_objs.gbc;
   gbv  = vacuum_objs.gbv;
   gpb  = vacuum_objs.gpb;
else
   gbc  = [];
   gbv  = [];
   gpb  = [];
end

if(isempty(ecdata))
   necoils=0;
   ecnturn = [];
else
   ecnturn = vacuum_objs.ecnturn;
   necoils = length(ecnturn);
end

% Scale pf for SC if needed:

if(length(scale_R_scpf)==1)
   scaleresc = [scale_R_scpf*ones(ncc,1)];
elseif(size(scale_R_scpf,1)==1)
   scaleresc = scale_R_scpf';
else
   scaleresc = scale_R_scpf;
end
if(length(scaleresc)~=ncc)
     wait('ERROR build_tokamak_system: length(scale_R_scpf) must = 1 or ncc')
     return;
end
resc = scaleresc.*vacuum_objs.resc;

%%%%%%%%%%%%%%%%%%%
% Equilibrium data
%%%%%%%%%%%%%%%%%%%

if(~vacuum_model)
  switch(ichooseq)
    case 1,
      if verbose > 1
        disp('Reading EFIT data...')
      end
      nesum = 1;
      equil_data = read_gfile_tok(efit_gfile,tokamak,[],[],[],[],cc_file);

    case 2,
      equil_data = read_corsica_data(corsica_inputs,vacuum_objs);

    case 3,


    case 4,


    otherwise,

  end
  if(isfield(build_inputs,'idx_efit_to_tok'))
    equil_data.idx_efit_to_tok = idx_efit_to_tok;
  end
  if(exist('efit_options','var')) 
    if(isfield(efit_options,'idxvv'))
      equil_data.vvid = equil_data.vvid(efit_options.idxvv);
      equil_data.vvfrac = equil_data.vvfrac(efit_options.idxvv);
      equil_data.vc = equil_data.vc(unique(equil_data.vvid));
      equil_data.idxvv_efit_to_tok = equil_data.idxvv_efit_to_tok(efit_options.idxvv); 
    end
  end

% jphi is in MA (kA?) here, while pcur and cpasma are in Amps

% Check that cc data is valid - otherwise model cannot be built.
   if(~isfield(equil_data,'cc') | length(equil_data.cc)==0)
      wait('ERROR build_tokamak_system: No cc data in equilibrium. Cannot build model.')
      return;
   end

% Check that equil_data and vacuum grids are consistent before going on.
   if(vacuum_objs.nr~=length(equil_data.rg) | vacuum_objs.nz~=length(equil_data.zg))
      fprintf('ERROR build_tokamak_system: equilibrium grid size %d x %d does not match vacuum size %d x %d\n', ...
		length(equil_data.rg), length(equil_data.zg), vacuum_objs.nr, vacuum_objs.nz);
      wait
      return;
   end

   drvac = mean(diff(vacuum_objs.rg));
   dzvac = mean(diff(vacuum_objs.zg));
   pct_diff = abs(drvac-equil_data.dr)/drvac*100;
   if(pct_diff>1e-3)
      disp(['WARNING build_tokamak_system: vacuum and equilibrium ' ...
			 'grid dr differ by ' num2str(pct_diff) '%'])
      if iwait, pause, end
   end
   pct_diff = abs(dzvac-equil_data.dz)/dzvac*100;
   if(pct_diff>1e-3)
      disp(['WARNING build_tokamak_system: vacuum and equilibrium ' ...
			 'grid dz differ by ' num2str(pct_diff) '%'])
      if iwait, pause, end
   end
end

%%%%%%%%%%%%%%%%%%%%%%
% plasma response data
%%%%%%%%%%%%%%%%%%%%%%

if(vacuum_model)
   dzdis = zeros(1,ncc+nvv);
   vc0 = [];
   rg = vacuum_objs.rg;
   zg = vacuum_objs.zg;
else
   idoplots=0;
   idoncal = 0;

% Generic radial/vertical stab calc: only needs cc,jphi,+environment,
% but need many things from this execution, eg dfsdr, cphi, etc...
% that are assumed present after this call.
% rzrig is done early so that we can use cc0 in vacuum calculations below.
   if verbose > 1
     disp('Entering rzrig calculation...')
   end
   [rzrig_data,cc0,vc0,dbg_rzrig] = rzrig(equil_data,tokamak,vacuum_objs,idoplots,idoncal,[],iwait,verbose);
   ip0 = rzrig_data.ip0;
   dfsdr = rzrig_data.dfsdr;
   dfsdz = rzrig_data.dfsdz;
%    vc0 = equil_data.iv;
%    ip0 = equil_data.cpasma;
%    cc0 = equil_data.cc;

    idxpl = find(equil_data.jphi(:)~=0);
    a0 = (max(vacuum_objs.rgg(idxpl))-min(vacuum_objs.rgg(idxpl)))/2;
    b0 = (max(vacuum_objs.zgg(idxpl))-min(vacuum_objs.zgg(idxpl)))/2;
    if(a0==b0)   % handle "point" distributions of current
       kap0=1;
    else
       kap0 =  b0/a0;
    end
    zcur = equil_data.jphi(:)'*vacuum_objs.zgg(:)/sum(equil_data.jphi(:));  %Zcentroid
    rcur = equil_data.jphi(:)'*vacuum_objs.rgg(:)/sum(equil_data.jphi(:));  %centroid major radius
    psivac = vacuum_objs.mpc*cc0;   %vacuum flux
    cphi = rzrig_data.cphi;
    cphi = equil_data.jphi*equil_data.dz*equil_data.dr;
    

   if(irzresp_dynamic < 4 | irzresp_output < 4)
      dzdis = scldzdis*rzrig_data.dzdis;
      dfpdr = sum(rzrig_data.dcdr(:).*psivac)/rzrig_data.ip0;
      dfpdz = sum(rzrig_data.dcdz(:).*psivac)/rzrig_data.ip0;
   end
   
   if(irzresp_dynamic == 4 | irzresp_output == 4)
      if verbose > 1
        disp('Entering corsica and vst calculation...')
      end
      if(ichooseq == 1)
        vst_data = plasma_response_vst(efit_gfile,tokamak,vacuum_objs,ichooseq);
      elseif(ichooseq == 3)
        vst_data = plasma_response_vst(corsica_savefile,tokamak,vacuum_objs,ichooseq);
      end
   end

   if(irzresp_dynamic == 5 | irzresp_output == 5)
      if verbose > 1
        disp('Entering gspert calculation...')
      end
      gspert_options.iconstraints = 1;
      if(~isfield(gspert_options,'idoplot'))
         gspert_options.idoplot = 0;
      end
      response = gspert(equil_data,vacuum_objs,gspert_options,gspert_options.idoplot);
      response.gspert_options = gspert_options;
      dzdis = response.dzdis;
   end
   
% Check that equilibrium cc0, vc0 are consistent with desired connectivity.
% HOW TO DO FOR vc0????
% HOW TO DO FOR netlist input???
    
   if(exist('cccirc'))
      test_idx=0;
      mean_cc0 = mean(abs(cc0));
      for ii=1:ncx
         idx = find(cccirc==ii);
         test = diff(cc0(idx));
         if(~isempty(test) & test~=0 & any(abs(test)/mean_cc0>1e-2))
            test_idx(ii,1:length(idx)) = idx;
         end
      end
      if(any(any(test_idx)))
            fprintf('WARNING build_tokamak_system:')
            fprintf('  Equilibrium not consistent with specified connectivity:\n');
            for ii=1:size(test_idx,1)
               if(any(test_idx(ii,:))~=0)
                  fprintf('Currents in coils '); fprintf(' %d ',test_idx(ii,:));
                  fprintf(' should be identical\n');
                  fprintf('Currents = '); fprintf('%f ', cc0(test_idx(ii,:))); fprintf('\n');
               end
            end
            wait
      end
   end

   psivac = vacuum_objs.mpc*cc0;   %vacuum flux 

   rg = equil_data.rg;
   zg = equil_data.zg;
   psizr = equil_data.psizr;
%   psibnd = equil_data.psibnd;	% psibry? Which is always defined? What's the difference?
   nz = equil_data.nh;
   nr = equil_data.nw;
end

if(~isempty(tok_geo_plot_fn))
   figure(1),clf,hold off
   if(strcmp(tok_geo_plot_fn,'generic'))
      options.iplteq=1; options.ipltflx=1;
      options.ipltbp=1; options.ipltfl=1;
      options.idxvv = 1:vacuum_objs.nv;
      if exist('vvgroup') & ~isempty(vvgroup)
         sz = size(vvgroup);
         if(min(sz)==1)
            options.idxvv = sort(setdiff(vvgroup(:),0));
         else
            options.idxvv = sort(setdiff(vvgroup(:,1),0));
         end
         options.vvgroup = vvgroup;
      else
         options.idxvv = 1:size(vvdata,2);
      end
      if(exist('cccirc')==1 & ~isempty(cccirc))
         if(isempty('ecnturn')~=1)
            idxec = 1:length(ecnturn);
         end
         idxfc = find(cccirc~=0);
         options.idxfc = setdiff(idxfc,idxec)-necoils;
      else
         options.idxfc = 1:size(fcdata,2);
      end
      plot_tok_geo(vacuum_objs,options);
   else
      iplteq=1; ipltflx=1;
      ipltbp=1; ipltfl=1;
      eval(tok_geo_plot_fn)
   end

   figure(2),clf,hold off
   if(strcmp(tok_geo_plot_fn,'generic'))
      options.iplteq=0; options.ipltflx=0;
      options.ipltbp=0; options.ipltfl=0;
      plot_tok_geo(vacuum_objs,options);
   else
      iplteq=0; ipltflx=0;
      ipltbp=0; ipltfl=0;
      eval(tok_geo_plot_fn)
   end
   if(~vacuum_model)
      title('Vacuum Field','FontSize',15)
      hold on
      psit = reshape(psivac,nz,nr);
      contour(rg,zg,psit,20,'w')
   end
end

% dzdis was set above, depending on scale parameter scldzdis
if(vacuum_model)
   drdis = zeros(1,ncc+nvv);
   drdip = 0;
   dzdip = 0;
   drdbetap = 0;
   dzdbetap = 0;
   drdli = 0;
   dzdli = 0;
   dfsdr = 0;
   dfsdz = 0;
elseif irzresp_output == 5 % gspert
   drdis = response.drdis;
   drdip = response.drdip;
   dzdip = response.dzdip;
   drdbetap = response.drdbetap;
   dzdbetap = response.dzdbetap;
   drdli = response.drdli;
   dzdli = response.dzdli;
else
   drdis = rzrig_data.drdis;
   drdip = rzrig_data.drdip;
   dzdip = rzrig_data.dzdip;
   drdbetap = rzrig_data.drdbetap;
   dzdbetap = rzrig_data.dzdbetap;
   drdli = rzrig_data.drdli;
   dzdli = 0;
   dfsdr = rzrig_data.dfsdr;
   dfsdz = rzrig_data.dfsdz;
end

if(~vacuum_model)

% Data for Ip circuit equation:

   [taup,Lp,Resp,eta] = calc_LR(rcur,a0,kap0,li_res,Te_res,Zeff_res);
   if(exist('Rp'))	% allow user override of plasma resistance
      Resp = Rp; 	%Resp = .5e-6;	% KLUGE - why is computed Resp so low???
   end
   if(~vacuum_objs.imks)
      Lp = Lp*1e6;     %H->uH
      Resp = Resp*1e6; %Ohm->uOhm
   end
%gradient in (vac)flux linked by plasma:
   mcIp = vacuum_objs.mpc'*equil_data.jphi(:)/sum(equil_data.jphi(:));
   mvIp = vacuum_objs.mpv'*equil_data.jphi(:)/sum(equil_data.jphi(:));

% Loop voltage Loop flux response data:

   if exist('lvdata')==1 & ~isempty(lvdata)  % dflvdz, dflvdr from rzrig
      if(irzresp_output==0)
         dflvdis = 0;  
         dflvdip = vacuum_objs.mph'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==1)
         dflvdis = rzrig_data.dflvdr*drdis;  
         dflvdip = rzrig_data.dflvdz*dzdip + vacuum_objs.mph'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==2)
         dflvdis = rzrig_data.dflvdz*dzdis;  
         dflvdip = rzrig_data.dflvdz*dzdip + vacuum_objs.mph'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==3)
         dflvdis = rzrig_data.dflvdz*dzdis + rzrig_data.dflvdr*drdis;  
         dflvdip = rzrig_data.dflvdz*dzdip + rzrig_data.dflvdr*drdip + vacuum_objs.mph'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==4)
         if upper(tokamak) == 'D3D'
	    dflvdis = rzrig_data.dflvdz*dzdis + rzrig_data.dflvdr*drdis; % For E coil response
	    dflvdis(3:end) = vst_data.dflvdis(3:end);
	 else
	    dflvdis = vst_data.dflvdis;
	 end
         dflvdip = rzrig_data.dflvdz*dzdip + rzrig_data.dflvdr*drdip + vacuum_objs.mph'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==5)
        dflvdis = mph'*response.dcphidis;
        dflvdip = mph'*response.dcphidip(:);
      end
   else
      dflvdis = [];
      dflvdip = [];
   end

% Flux Loop response data:

   if exist('fldata')==1 & ~isempty(fldata)
      if(irzresp_output==0)
         dfldis = 0;  
         dfldip = vacuum_objs.mpl'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==1)
         dfldis = rzrig_data.dfldr*drdis;  
         dfldip = rzrig_data.dfldr*drdip + vacuum_objs.mpl'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==2)
         dfldis = rzrig_data.dfldz*dzdis;
         dfldip = rzrig_data.dfldz*dzdip + vacuum_objs.mpl'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==3)
         dfldis = rzrig_data.dfldz*dzdis + rzrig_data.dfldr*drdis;  
         dfldip = rzrig_data.dfldz*dzdip + rzrig_data.dfldr*drdip + vacuum_objs.mpl'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==4)
         if upper(tokamak) == 'D3D'
	    dfldis = rzrig_data.dfldz*dzdis + rzrig_data.dfldr*drdis; % For E coil response
	    dfldis(3:end) = vst_data.dfldis(3:end);
	 else
	    dfldis = vst_data.dfldis;
	 end
         dfldip = rzrig_data.dfldz*dzdip + rzrig_data.dfldr*drdip + vacuum_objs.mpl'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==5)
        dfldis = mpl'*response.dcphidis;
        dfldip = mpl'*response.dcphidip(:);
      end
   else
      dfldis = [];
      dfldip = [];
   end
   
% Bprobe response data:

   if exist('bpdata')==1 & ~isempty(bpdata)
% dbpdz,dbpdr from rzrig
      if(irzresp_output==0)
         dbpdis = 0;  
         dbpdip = vacuum_objs.gpb'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==1)
         dbpdis = rzrig_data.dbpdr*drdis; 
         dbpdip = rzrig_data.dbpdr*drdip + vacuum_objs.gpb'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==2)
         dbpdis = rzrig_data.dbpdz*dzdis; 
         dbpdip = rzrig_data.dbpdz*dzdip + vacuum_objs.gpb'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==3)
         dbpdis = rzrig_data.dbpdz*dzdis + rzrig_data.dbpdr*drdis; 
         dbpdip = rzrig_data.dbpdz*dzdip + rzrig_data.dbpdr*drdip + vacuum_objs.gpb'*rzrig_data.cphi(:)/rzrig_data.ip0;  
       elseif(irzresp_output==4)
         if upper(tokamak) == 'D3D'
            dbpdis = rzrig_data.dbpdz*dzdis + rzrig_data.dbpdr*drdis;  % For E coil response
	    dbpdis(3:end) = vst_data.dbpdis(3:end);
	 else
	    dbpdis = vst_data.dbpdis;
	 end
         dbpdip = rzrig_data.dbpdz*dzdip + rzrig_data.dbpdr*drdip + vacuum_objs.gpb'*rzrig_data.cphi(:)/rzrig_data.ip0;  
      elseif(irzresp_output==5)
        dbpdis = gpb'*response.dcphidis;
	dbpdip = gpb'*response.dcphidip(:);
      end
   else 
      dbpdis = 0;
      dbpdip = 0;
   end

% MSE response data:

   if nmse>0
% dBrMSEdr, dBzMSEdr, dBrMSEdz, dBzMSEdz, from rzrig
      if(irzresp_output==0)
         dBrMSEdis = 0;
	 dBzMSEdis = 0;
	 dBrMSEdip = vacuum_objs.gmsebrp*rzrig_data.cphi(:)/rzrig_data.ip0; 
	 dBzMSEdip = vacuum_objs.gmsebzp*rzrig_data.cphi(:)/rzrig_data.ip0; 
      elseif(irzresp_output==1)
         dBrMSEdis = rzrig_data.dBrMSEdr*drdis;
	 dBzMSEdis = rzrig_data.dBzMSEdr*drdis; 
	 dBrMSEdip = rzrig_data.dBrMSEdr*drdip + vacuum_objs.gmsebrp*rzrig_data.cphi(:)/rzrig_data.ip0; 
	 dBzMSEdip = rzrig_data.dBzMSEdr*drdip + vacuum_objs.gmsebzp*rzrig_data.cphi(:)/rzrig_data.ip0; 
      elseif(irzresp_output==2)
         dBrMSEdis = rzrig_data.dBrMSEdz*dzdis;
	 dBzMSEdis = rzrig_data.dBzMSEdz*dzdis; 
	 dBrMSEdip = rzrig_data.dBrMSEdz*dzdip + vacuum_objs.gmsebrp*rzrig_data.cphi(:)/rzrig_data.ip0; 
	 dBzMSEdip = rzrig_data.dBzMSEdz*dzdip + vacuum_objs.gmsebzp*rzrig_data.cphi(:)/rzrig_data.ip0; 
      elseif(irzresp_output==3)
         dBrMSEdis = rzrig_data.dBrMSEdr*drdis + rzrig_data.dBrMSEdz*dzdis;
	 dBzMSEdis = rzrig_data.dBzMSEdr*drdis + rzrig_data.dBzMSEdz*dzdis; 
	 dBrMSEdip = rzrig_data.dBrMSEdr*drdip + rzrig_data.dBrMSEdz*dzdip + vacuum_objs.gmsebrp*rzrig_data.cphi(:)/rzrig_data.ip0; 
	 dBzMSEdip = rzrig_data.dBzMSEdr*drdip + rzrig_data.dBzMSEdz*dzdip + vacuum_objs.gmsebzp*rzrig_data.cphi(:)/rzrig_data.ip0; 
      elseif(irzresp_output==4) % use rzrig for now
         dBrMSEdis = rzrig_data.dBrMSEdr*drdis + rzrig_data.dBrMSEdz*dzdis;
	 dBzMSEdis = rzrig_data.dBzMSEdr*drdis + rzrig_data.dBzMSEdz*dzdis; 
	 dBrMSEdip = rzrig_data.dBrMSEdr*drdip + rzrig_data.dBrMSEdz*dzdip + vacuum_objs.gmsebrp*rzrig_data.cphi(:)/rzrig_data.ip0; 
	 dBzMSEdip = rzrig_data.dBzMSEdr*drdip + rzrig_data.dBzMSEdz*dzdip + vacuum_objs.gmsebzp*rzrig_data.cphi(:)/rzrig_data.ip0; 
      elseif(irzresp_output==5)
        dBrMSEdis = gmsebrp*response.dcphidis;
        dBzMSEdis = gmsebzp*response.dcphidis;
        dBrMSEdip = gmsebrp*response.dcphidip(:);
        dBzMSEdip = gmsebzp*response.dcphidip(:);
      end
   else 
      dBrMSEdis = 0;
      dBzMSEdis = 0;
      dBrMSEdip = 0;
      dBzMSEdip = 0;
   end
else
   dfldis = 0;
   dfldip = 0;
   dflvdis = 0;
   dflvdip = 0;
   dbpdis = 0;
   dbpdip = 0;
   dBrMSEdis = 0;
   dBzMSEdis = 0;
   dBrMSEdip = 0;
   dBzMSEdip = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble/modify data objects needed for state equation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(vacuum_model)
   ns = ncc+nvv;
   xmats = zeros(ns);
   xmatsp = zeros(ns,1);
   xmatps = zeros(1,ns);
   xmatpp = 0;
   xmatsb = zeros(ns,1);
   xmatsl = zeros(ns,1);
   xmatpb = 0;
   xmatpl = 0;
else
   if(irzresp_dynamic==0)				% no plasma motion
      ns = ncc+nvv;
      xmats = zeros(ns);
      xmatsp = zeros(ns,1);
      xmatps = zeros(1,ns);
      xmatpp = 0;
      xmatsb = zeros(ns,1);
      xmatsl = zeros(ns,1);
      xmatpb = 0;
      xmatpl = 0;

   elseif(irzresp_dynamic==1)				% rigid R motion only
      xmats = rzrig_data.xmatr;
      xmatsp = dfsdr*drdip;
      xmatps = rzrig_data.dfpdrC*drdis;
      xmatpp = rzrig_data.dfpdrC*drdip;
      xmatsb = dfsdr*drdbetap;
      xmatsl = dfsdr*drdli;
      xmatpb = dfpdr*drdbetap; 				% Wb/(unit betap)
      xmatpl = 1e-6*mu0*rcur*ip0/2; %dPsi_p/dli: Wb/unit-li

   elseif(irzresp_dynamic==2)				% rigid Z motion only
      ns = ncc+nvv;
% Must use dfsdz*dzdis here instead of xmatz, since dzdis may be rescaled.
      xmats = rzrig_data.dfsdz*dzdis;
      xmatsp = dfsdz*dzdip;
      xmatps = rzrig_data.dfpdzC*dzdis;
      xmatpp = rzrig_data.dfpdzC*dzdip;
      xmatsb = dfsdz*dzdbetap;
      xmatsl = zeros(ns,1);
      xmatpb = dfpdz*dzdbetap; 				% Wb/(unit betap)
      xmatpl = 1e-6*mu0*rcur*ip0/2; %dPsi_p/dli: Wb/unit-li

   elseif(irzresp_dynamic==3)				% rigid R and Z motion
% Must use dfsdz*dzdis here instead of xmatz, since dzdis may be rescaled.
      xmats = rzrig_data.dfsdz*dzdis + rzrig_data.xmatr;
      xmatsp = dfsdr*drdip+dfsdz*dzdip;
      xmatps = rzrig_data.dfpdrC*drdis+rzrig_data.dfpdzC*dzdis;
      xmatpp = rzrig_data.dfpdrC*drdip+rzrig_data.dfpdzC*dzdip;
      xmatsb = dfsdr*drdbetap + dfsdz*dzdbetap;
      xmatsl = dfsdr*drdli;
      xmatpb = dfpdr*drdbetap + dfpdz*dzdbetap;          % Wb/(unit betap)
      xmatpl = 1e-6*mu0*rcur*ip0/2; %dPsi_p/dli: Wb/unit-li

   elseif(irzresp_dynamic==4)
      xmats = rzrig_data.dfsdz*dzdis + rzrig_data.xmatr;
      xmatsp = dfsdr*drdip+dfsdz*dzdip;
      xmatps = rzrig_data.dfpdrC*drdis+dfpdzC*dzdis;
      xmatpp = rzrig_data.dfpdrC*drdip+dfpdzC*dzdip;
      xmatsb = dfsdr*drdbetap + dfsdz*dzdbetap;
      xmatsl = dfsdr*drdli;
      xmatpb = dfpdr*drdbetap + dfpdz*dzdbetap;          % Wb/(unit betap)
      xmatpl = 1e-6*mu0*rcur*ip0/2; %dPsi_p/dli: Wb/unit-li
      if upper(tokamak) == 'D3D'
	 xmats(3:end,3:end) = vst_data.xmats(3:end,3:end);
      else
	 xmats = vst_data.xmats;
      end
   elseif(irzresp_dynamic==5)
      xmats  = response.xmats;
      xmatsp = [vacuum_objs.mpc vacuum_objs.mpv]'*response.dcphidip(:)-[mcIp;mvIp]; % First term is both direct and indirect coupling
      xmatsb = [vacuum_objs.mpc vacuum_objs.mpv]'*response.dcphidbetap(:);
      xmatsl = [vacuum_objs.mpc vacuum_objs.mpv]'*response.dcphidli(:);
      xmatps = response.dpsipladis-[mcIp;mvIp]';
      xmatpb = response.dpsipladbetap;
      xmatpl = response.dpsipladli;
      xmatpp = response.dpsipladip-Lp;
   end
end

if iplcirc==0
   mxx = [[vacuum_objs.mcc vacuum_objs.mcv];[vacuum_objs.mcv' vacuum_objs.mvv]];
   rxx = diag([resc; vacuum_objs.resv]);
   if(vacuum_model)
      xmatx = 0; 
   else
      xmatx = xmats; 
   end
elseif iplcirc==1
   mxx = [[vacuum_objs.mcc vacuum_objs.mcv mcIp]; ...
	 [vacuum_objs.mcv' vacuum_objs.mvv mvIp]; ...
	 [mcIp' mvIp' Lp]]; 
   rxx = diag([resc; vacuum_objs.resv; Resp]);

   xmatx = [[xmats  xmatsp]; ...
	    [xmatps xmatpp] ];
else
   disp('Invalid iplcirc!!!')
   return
end
nxx = size(mxx,1);

if ~vacuum_model
temp = (xmatx-xmatx')./mean(mean(abs(xmatx)));
s1 = size(temp,1);
temp1 = temp(:);
idx = find(temp<tol_symmetry);
temp1(idx) = 0;
temp = reshape(temp1,s1,s1);
temp1 = max(max(abs(temp)));
if(temp1 > tol_symmetry) & verbose > 0
   figure(100001)
   plot(xmatx(end,:))
   hold on
   plot(xmatx(:,end),'r--')
   hold off
   title('large relative differences between xmatx and xmatx^T')
   fprintf(['ERROR build_tokamak_system: lstar not symmetric, ' ...
           'max relative error = %g\n'],temp1)
   fprintf('(Plot shows locations of non-symmetric entries.)\n');
   if iwait, wait, end
end
end

% Modifications to R and/or M matrix:
% Previously, we used Rckt and Lctk, which were "circuit size" but this 
% required user of code to understand internals of modeling. Now we use
% Rext and Lext, which only requires that they know the size of the objects 
% (M and R) they are actually passing in.

Rextra = [];
if(Rext==0)
   Rextra = zeros(ncc+nvv,1);
elseif(length(Rext)~=ncc & length(Rext)~=ncc+nvv)
   wait('ERROR build_tokamak_system: length(Rext) must = ncc or ncc+nvv')
   return;
else
   if(any(replace_Rext))
      if(length(replace_Rext)==1)
         replace_Rext = replace_Rext*ones(size(Rext));
      elseif(length(replace_Rext)~=length(Rext))
         wait('ERROR build_tokamak_system: length of replace_Rext must be 1 or length(Rext)');
         return;
      end
      Rextra = zeros(ncc+nvv,1);
      for k=1:length(Rext)
         Rextra(k) = Rext(k);
         if(replace_Rext(k))
            rxx(k,k) = 0;
         end
      end
   else
      if(length(Rext)==ncc)
         Rextra = [Rext; zeros(nvv,1)];
      else
         Rextra = Rext;
      end
   end
end
Lextra = [];
if(Lext==0)
   Lextra = zeros(ncc+nvv,1);
elseif(length(Lext)==ncc)
   Lextra = [Lext; zeros(nvv,1)];
elseif(length(Lext)==ncc+nvv)
   Lextra = Lext;
else
   wait('ERROR build_tokamak_system: length(Lext) must = ncc or ncc+nvv')
   return;
end
if(iplcirc)
   Rextra = [Rextra; 0]; Lextra = [Lextra; 0];
end
Rextra = diag(Rextra); Lextra = diag(Lextra);

if(use_cccirc)
   Rckt_extra = [];
   if(Rckt==0)
      Rckt_extra = zeros(ncx+nvx,1);
   elseif(length(Rckt)~=ncx & length(Rckt)~=ncx+nvx)
      wait('ERROR build_tokamak_system: length(Rckt) must = ncx or ncx+nvx')
      return;
   else
      if(replace_Rckt)
         Rckt_extra = zeros(ncx+nvx,1);
         if(length(Rckt)==ncx)
            for k=1:ncx
               Rckt_extra(k) = Rckt(k);
            end
            for k=1:ncc
               rxx(k,k) = 0;
            end
         else
            for k=1:ncx+nvx
               Rckt_extra(k) = Rckt(k);
            end
            for k=1:ncc+nvv
               rxx(k,k) = 0;
            end
         end
      else
         if(length(Rckt)==ncx)
            Rckt_extra = [Rckt; zeros(nvx,1)];
         else
            Rckt_extra = Rckt;
         end
      end
   end

   Lckt_extra = [];
   if(Lckt==0)
      Lckt_extra = zeros(ncx+nvx,1);
   elseif(length(Lckt)~=ncx)
      wait('ERROR build_tokamak_system: length(Lckt) must = reduced size')
      return;
   else
      Lckt_extra = [Lckt; zeros(nvx,1)];
   end

   if(iplcirc)
      Rckt_extra = [Rckt_extra; 0]; Lckt_extra = [Lckt_extra; 0];
   end
   Rckt_extra = diag(Rckt_extra); Lckt_extra = diag(Lckt_extra);
end

%------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now build the model (state equations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If plasma is included, incorporate effects in mutual inductance:

if(vacuum_model)
   L_star = mxx;
else
   L_star = mxx + xmatx;
end

if(use_cccirc)

   if(~isempty(cccirc)) 	% Circuit connections selected by cccirc
      if verbose>0
         disp('Warning: CC circuit connections being modified by cccirc!!!')
      end
      Pcc = zeros(ncc,ncx);  %Pcc maps from conn-circs to orig unconn-cir

      nturn = [ecnturn;vacuum_objs.fcnturn];
      for ii=1:ncx
        idx0=find(abs(cccirc)==ii);
        idx1=find(cccirc==ii);
        idx2=find(cccirc==-ii);
        if(vacuum_objs.iterminal)
           turns_ratio1 = 1;
           turns_ratio2 = 1;
        else
           turns_ratio1 = nturn(idx1)/sum(nturn(idx0));
           turns_ratio2 = nturn(idx2)/sum(nturn(idx0));
        end
        if ~isempty(idx1)
           Pcc(idx1,ii)=ones(length(idx1),1).*turns_ratio1;
        end
        if ~isempty(idx2)
           Pcc(idx2,ii)=-ones(length(idx2),1).*turns_ratio2;
        end
      end
      temp_file = [tokamak '_netlist.dat'];
      [netlist,inode1,inode2] = cccirc_to_netlist(cccirc,vacuum_objs,temp_file);

   else
      Pcc = eye(ncc);
      temp_file = [tokamak '_netlist.dat'];
      [netlist,inode1,inode2] = cccirc_to_netlist(cccirc,vacuum_objs,temp_file);
   end

%Pxx maps from full state vector *with* connected circs to unconn-circs:

   Pvv = eye(nvv);  
     if verbose > 1
        disp('VV circuit connections being modified by vvgroup...')
     end
     if(size(vvgroup,2)==1)       % distribute current according to resistance
        calc_vvfrac = 1;
        vvfrac = zeros(length(vvgroup),1);
     elseif(size(vvgroup,2)==2)   % user-specified sharing of current
        calc_vvfrac = 0;
        vvfrac = vvgroup(:,2);
     end
     vvid = vvgroup(:,1);
     nvg =max(abs(vvid));	% nvg = number of grouped vv elts
     Pvv = zeros(nvv,nvg);  
  
     for ii=1:nvg
        idx1=find(vvid==ii); 
        if ~isempty(idx1)
           if(calc_vvfrac)
              sum_rinv = sum(1./resv(idx1));
              vvfrac(idx1) = 1./resv(idx1)/sum_rinv;
           elseif(abs(sum(vvfrac(idx1))-1)>1e-3);
              fprintf(['WARNING build_tokamak_system: vessel group current ' ...
                        'sharing does not sum to 1\n  vessel indices = ']);
              fprintf('%d ',idx1);
              fprintf('\n');
              if iwait, wait, end
           end
           Pvv(idx1,ii)=vvfrac(idx1);
        end
        idx1=find(vvid==-ii);
        if ~isempty(idx1)
           if(calc_vvfrac)
              sum_rinv = sum(1./resv(idx1));
              vvfrac(idx1) = 1./resv(idx1)/sum_rinv;
           elseif(abs(sum(vvfrac(idx1))-1)>1e-3);
              fprintf(['WARNING build_tokamak_system: vessel group current ' ...
                        'sharing does not sum to 1\n  vessel indices = ']);
              fprintf('%d ',idx1);
              fprintf('\n');
              if iwait, wait, end
           end
           Pvv(idx1,ii)=-vvfrac(idx1);
        end
     end

   if(~isempty(vvcirc))
      temp = zeros(nvg,nvx);
      for ii=1:nvx
        idx1=find(vvcirc==ii);
        idx2=find(vvcirc==-ii);
        if ~isempty(idx1)
           temp(idx1,ii)=ones(length(idx1),1);
        end
        if ~isempty(idx2)
           temp(idx2,ii)=-ones(length(idx2),1);
        end
      end
      Pvv = Pvv * temp;
   end
   Pxx = [[Pcc zeros(ncc,nvx)]; [zeros(nvv,ncx) Pvv]];

   if(iplcirc)
      Pxx = [[Pxx zeros(size(Pxx,1),1)]; [zeros(1,size(Pxx,2)) 1]];
   end
% Amatrix:
   lstar = Pxx'*(L_star+Lextra)*Pxx + Lckt_extra;
   lstari = inv(lstar);
   temp = (Pxx'*(rxx+Rextra)*Pxx+Rckt_extra);
   temp1 = min(temp); [mm,ii] = min(temp1);
   if mm<0
      disp('ERROR build_tokamak_system: ')
      disp(['Consolidation of resistances has ' ...
        'resulted in NEGATIVE resistance value for circuit ' int2str(ii)])
      if iwait, pause, end
      return;
   end
   amat = -lstari*(Pxx'*(rxx+Rextra)*Pxx+Rckt_extra);
   nxx_cc = size(amat,1);    %nxx after cccirc connection...

% Bmatrix:
   vmat=zeros(nxx_cc,ncx); vmat(1:ncx,1:ncx)=eye(ncx);  %drive only circuits...
   bmat = lstari*vmat;

else	% connections specified by input netlist 
% L_star includes plasma, if iplcirc==1
   if(iplcirc)
      idx = strmatch('MIp',netlist.names);
      if(isempty(idx))
         wait('ERROR build_tokamak_system: netlist must contain name MIp if iplcirc=1');
         return;
      end
   end
   netlist_model = model_from_netlist(netlist,L_star+Lextra,...
		diag(rxx+Rextra),netlist_currents,netlist_voltages);
   rxx = netlist_model.rxx;   % now includes all resistance in series with coil
   Pxx = netlist_model.Pxx;
   ncx = sum(max(abs(Pxx(1:ncc,:))));	% number of coil states
   Pcc = Pxx(1:ncc,1:ncx);
   nvx = size(Pxx,2)-ncx-iplcirc;		% number of vessel states
   Pvv = Pxx(ncc+1:ncc+nvv,ncx+1:ncx+nvx);
   lstar = netlist_model.Mhat;
   amat = -inv(lstar) * netlist_model.Rhat;
   bmat = inv(lstar) * netlist_model.Vhat;
   nout_nl = size(netlist_model.Chat,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the output equations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We first construct the "reduced output" cmat and dmat, then
% use projections to generate full-size cmat_allPF and dmat_allPF.
% If there are voltage or current outputs from netlist, these must be 
% inserted after collapse to representations in terms of current states.

Pcci = pinv(Pcc);   %inverse to get new currents at output of cmat...
Pxxi = pinv(Pxx);

if iplcirc==0
   if verbose > 1
      disp(['cmat outputs=[Icx(ncx); psi_fl(nfl); Bp(nbp); I_rog(nrl); Vloop(nlv); Br(nmse); Bz(nmse); R; Z]'])
   end
   tmp1=zeros(ncx,nxx); tmp1(1:ncx,1:ncc)=Pcci*eye(ncc); %extract Ic f/ states  
   cmat = [tmp1; ([mlc mlv]+dfldis)];
   if(~isempty(vacuum_objs.gpb))
      cmat = [cmat; [gbc gbv]+dbpdis];
   end
   if nrl > 0  % add Rogowski loops
      tmp = rldata;
      tmp(1:ncc,:) = tmp(1:ncc,:) .* repmat(ccnturn,1,nrl);
      tmp(end,:) = []; % delete Ip reading
      cmat = [cmat; tmp'];
   end
   dmat = zeros(size(cmat,1),size(bmat,2));
   if nlv > 0 % add loop voltages
     cmat = [cmat; -([mhc mhv]+dflvdis)*Pxx*amat*Pxxi]; 
     dmat = [dmat; -([mhc mhv]+dflvdis)*Pxx*bmat]; 
   end
   if nmse > 0 % Add (2*nmse) MSE outputs
      cmat = [cmat; [gmsebrc gmsebrv]+dBrMSEdis; [gmsebzc gmsebzv]+dBzMSEdis];
      dmat = [dmat; zeros(2*nmse,size(dmat,2))];
   end
   if(use_netlist & nout_nl)
      cmat = [cmat; zeros(nout_nl,size(cmat,2))];
      dmat = [dmat; zeros(nout_nl,size(dmat,2))];
   end
   if(~vacuum_model)
      cmat = [cmat; drdis; dzdis];
      dmat = [dmat; zeros(2,size(dmat,2))];
   end
   cmat = cmat*Pxx;   %map currents to states
   if(use_netlist & nout_nl)
      cmat(end-nout_nl+1:end,:) = netlist_model.Chat;
      dmat(end-nout_nl+1:end,:) = netlist_model.Dhat;
   end
elseif iplcirc==1
   if verbose > 1
      disp('cmat outputs=[Icx(ncx); psi_fl(nfl); Bp(nbp); I_rog(nrl); Vloop(nlv); Br(nmse); Bz(nmse); R; Z; Ip]')
   end
   tmp1=zeros(ncx,nxx); tmp1(1:ncx,1:ncc)=Pcci*eye(ncc); %extract Ic f/ states
   tmp2=zeros(1,nxx); tmp2(end)=1;  %to extract Ip from state vector
   cmat = [tmp1; [[mlc mlv]+dfldis dfldip]];
   if(~isempty(vacuum_objs.gpb))
      cmat = [cmat; [[gbc gbv]+dbpdis dbpdip]];
   end 
   if(nrl > 0)  % add Rogowski loops
      tmp = rldata;
      tmp(1:ncc,:) = tmp(1:ncc,:) .* repmat(ccnturn,1,nrl);
      cmat = [cmat; tmp'];
   end
   dmat = zeros(size(cmat,1),size(bmat,2));
   if(nlv > 0) % add loop voltages
      cmat = [cmat;-[[mhc mhv]+dflvdis dflvdip]*Pxx*amat*Pxxi]; 
      dmat = [dmat;-[[mhc mhv]+dflvdis dflvdip]*Pxx*bmat]; 
   end
   if nmse > 0 % Add (2*nmse) MSE outputs
      cmat = [cmat; [[gmsebrc gmsebrv]+dBrMSEdis dBrMSEdip]; ...
                    [[gmsebzc gmsebzv]+dBzMSEdis dBzMSEdip] ];
      dmat = [dmat; zeros(2*nmse,size(dmat,2))];
   end
   if(use_netlist & nout_nl)
      cmat = [cmat; zeros(nout_nl,size(cmat,2))];
      dmat = [dmat; zeros(nout_nl,size(dmat,2))];
   end
   cmat = [cmat; ...
           [drdis drdip]; [dzdis dzdip]; tmp2];
   dmat = [dmat; zeros(3,size(dmat,2))];
   cmat = cmat*Pxx;   %map currents to states
   if(use_netlist & nout_nl)
      cmat(end-nout_nl-2:end-3,:) = netlist_model.Chat;
      dmat(end-nout_nl-2:end-3,:) = netlist_model.Dhat;
   end
end
if verbose > 1
   disp('(Icx(1:ncx) are cccirc-connected circuit currents)')
end

% Cmatrix with all unconnected PF currents output:

nyy = size(cmat,1);  %y vec with conn circ
%Pyy maps from output vector y with connected circs to y with unconn-circs:
Pyy = [[Pcc zeros(ncc,nyy-ncx)]; [zeros(nyy-ncx,ncx) eye(nyy-ncx)]];
cmat_allPF = Pyy*cmat;
dmat_allPF = Pyy*dmat;


if(vacuum_model)
   fmat = [];
   hmat_allPF = [];
   disturbances = '';
else

   disturbances = strvcat('betap', 'li', 'dbetap/dt', 'dli/dt');

% F matrix = multipliers of dbetap/dt and dli/dt to get psi_dot

   temp = [zeros(size(xmatsb,1),2) xmatsb xmatsl];
   if(iplcirc)
      temp = [temp; [zeros(1,2) xmatpb xmatpl]];
   end
   fmat = -inv(lstar) * Pxx' * temp;

% ****** Will this ALWAYS be correct for netlist modeling?  ******

% H matrix - multipliers of delta_betap and delta_li to output

% dcphidbetap units = A/m * m/beta = A/beta
   if irzresp_output == 5 % In this case use gspert response
     dcphidbetap = response.dcphidbetap(:); dcphidli = response.dcphidli(:);
   else % rzrig response
     dcphidbetap = (rzrig_data.dcdr(:)*drdbetap + rzrig_data.dcdz(:)*dzdbetap);
     dcphidli = (rzrig_data.dcdr(:)*drdli + rzrig_data.dcdz(:)*dzdli);
   end
%    tmp1=zeros(ncx,4);
%    hmat = [tmp1; [mpl'*dcphidbetap mpl'*dcphidli zeros(nfl,2)]];
%    if(~isempty(vacuum_objs.gpb))
%       hmat = [hmat; [gpb'*dcphidbetap gpb'*dcphidli zeros(nbp,2)]];
%    end
%    if(nrl > 0)  % add Rogowski loops
%       hmat = [hmat; zeros(nrl,size(hmat,2))];
%    end
%    if(nlv > 0) % add loop voltages
%      if iplcirc==1
%        hmat = [hmat;-[[mhc mhv]+dflvdis dflvdip]*Pxx*fmat];
%      else
%        hmat = [hmat;-[[mhc mhv]+dflvdis]*Pxx*fmat];
%      end
%    end
%    if nmse > 0
%       hmat = [hmat; [gmsebrp*dcphidbetap gmsebrp*dcphidli zeros(nmse,2)]];
%       hmat = [hmat; [gmsebzp*dcphidbetap gmsebzp*dcphidli zeros(nmse,2)]];
%    end
%    if(use_netlist & nout_nl)
%       hmat = [hmat; zeros(nout_nl,size(hmat,2))];
%    end
%    hmat = [hmat; [drdbetap drdli zeros(1,2)]; [dzdbetap dzdli zeros(1,2)]];
%    if iplcirc==1
%       hmat = [hmat; zeros(1,4)];
%    end
%    hmat_allPF = Pyy*hmat;
end

% Consistency check:

[vecs,vals] = eigsort(amat);
if verbose > 1
    disp(['After circuit connections:'])
    disp('Reconnected System is: amat,bmat,cmat,dmat...')
    disp('Reconnected System dimensions:')
    disp(['nx = ',int2str(size(amat,1))])
    disp(['nu = ',int2str(size(bmat,2))])
    disp(['ny = ',int2str(size(cmat,1))])
    disp(['#of driven circuits, ncx = ',int2str(ncx)])
    disp(['Max. OL eigenvalue for connected system= ',num2str(max(real(vals)))])
end
if ~all(isreal(vals)) & verbose > 0
  disp('ERROR build_tokamak_system: COMPLEX EIGENVALUES!!!!')
  if iwait, pause, end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equilibrium conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate equilibrium set of output values.

% Compute conductor contributions first. (Plasma contribution added below.)

% Determination of the "equilibrium loop voltage" is a bit confusing, since it
% seems that this value is not in fact even determined by the plasma 
% equilibrium, but rather by the external circuit.  What seems to be most
% consistent in terms of matching the coordinate transformations in the valid
% coords paper however, is to set the "equilibrium loop voltage" equal to 0.

if(vacuum_model)	% equilibrium of vacuum system = 0
   Ieq = zeros(ncx+nvx,1); 
   yeq = zeros(ncc+nfl+nbp+nrl+nlv+2*nmse,1); 
   if(use_netlist & nout_nl)
      yeq = [yeq; zeros(nout_nl,1)];
   end
else

% Set loop voltage yeq to NaN because loop voltage is not 
% defined by the equilibrium.
   lveq = nan*ones(size(mhc,1),length([cc0; vc0]));

% SHN: Setting loop voltage yeq to NaN blows the KSTAR simulator up
% Setting them to zero
%    lveq = 0.*ones(size(mhc,1),length([cc0; vc0]));


   if(iplcirc)
% This should provide correct Ieq regardless of whether cccirc or netlist is used
% to model, since full set of equil. currents = Pxx * Ieq in either case.
      Ieq = pinv(Pxx)*[cc0; vc0; ip0];
% NOTE that currents represented by Ieq may NOT be exactly equal to the corresponding 
% currents in cc0 if the appropriate circuit constraints are not imposed during equilibrium
% reconstruction (e.g. VFI constraint on D3D).
      Ic_eq = Pxx*Ieq; Ic_eq = Ic_eq(1:ncc);
      maxdiff = max(abs(Ic_eq - cc0));
      if(maxdiff > 20)
         fprintf('Ieq constrained by circuit connections differs from input (unconstrained) equilibrium.\n');
         fprintf('max diff = %d, average diff = %d\n',round(maxdiff),round(mean(abs(Ic_eq - cc0))));
      end

      yip = ip0;
      yeq = [[eye(ncc) zeros(ncc,nxx-ncc-1)]; ...	% (leave out plasma)
      		[vacuum_objs.mlc vacuum_objs.mlv]; ...  % full PF currents
		[vacuum_objs.gbc vacuum_objs.gbv]; ...
		rldata(1:end-1,:)'; ...
		lveq; ...
		[gmsebrc gmsebrv]; [gmsebzc gmsebzv]] * [cc0; vc0];

% Equilibrium branch currents are given by the netlist model, but branch voltage 
% is not defined by the equilbrium.
      [s1,s2]=size(yeq);
      I = [cc0; vc0; ip0];
      if(use_netlist & nout_nl)
         for i = 1:nout_nl
            if(strncmp(netlist_model.outputs(i,1),'V',1))
               yeq = [yeq; nan];
            else
               yeq = [yeq; netlist_model.Chat(i,:)*I(netlist_model.Istate_idx)];
            end
         end
      end
      yeq = [yeq; zeros(3,s2)];		% put r, z, and Ip outputs at end
   else
      Ieq = pinv(Pxx)*[cc0; vc0];
      yip = [];
      yeq = [[eye(ncc) zeros(ncc,nxx-ncc)]; ...
      		[vacuum_objs.mlc vacuum_objs.mlv]; ...
		[vacuum_objs.gbc vacuum_objs.gbv]; ...
		rldata(1:end-1,:)'; ...
   		lveq; ...
		[gmsebrc gmsebrv]; [gmsebzc gmsebzv]] * [cc0; vc0];

% Equilibrium branch currents are given by the netlist model, but branch voltage 
% is not defined by the equilbrium.
      [s1,s2]=size(yeq);
      I = [cc0; vc0; ip0];
      if(use_netlist & nout_nl)
         for i = 1:nout_nl
            if(strncmp(netlist_model.outputs(i,1),'V',1))
               yeq = [yeq; nan];
            else
               yeq = [yeq; netlist_model.Chat(i,:)*I(netlist_model.Istate_idx)];
            end
         end
      end
      yeq = [yeq; zeros(2,s2)];		% put r and z outputs at end
   end
end

% then add plasma contribution:

if(~vacuum_model)
   if(~isempty(mpl))
      yfl = mpl'*cphi(:); % flux loop contribution from plasma
   else
      yfl = [];
   end
   if(~isempty(gpb))
      ybp = gpb'*cphi(:); % B-probe contribution from plasma
   else
      ybp = [];
   end
   if nlv > 0
      ylv = zeros(size(mhc,1),1);	% loop voltage contribution from plasma
   else
      ylv = [];
   end
   if nmse > 0
      ymsebr = gmsebrp*cphi(:); % B-probe contribution from plasma
      ymsebz = gmsebzp*cphi(:); % B-probe contribution from plasma
   else
      ymsebr = []; ymsebz = [];
   end

   idx_R = ncx+vacuum_objs.nfl+vacuum_objs.nbp+1; 
   yeq(idx_R:idx_R+1)=0;

   yrog = zeros(nrl,1) + ip0;

   if(use_netlist & nout_nl)
      yckt = zeros(nout_nl,1);
   else
      yckt=[];
   end

% Add the plasma contribution to the vacuum contribution (in yeq).
   yeq = yeq + [zeros(ncc,1); yfl; ybp; yrog; ylv; ymsebr; ymsebz; yckt; rcur; zcur; yip];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Self-documentation for the model data structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output labels:

all_outputs = ''; 
coilnames_given = isfield(names,'coils');
if(coilnames_given)
   coilnames_given = ~isempty(names.coils);
end
if(coilnames_given)
   if(size(names.coils,1)~=vacuum_objs.nc)
      wait(['ERROR build_tokamak_system: must be ' int2str(vacuum_objs.nc) ...
		' names in names.coils'])
      return;
   end
end

if(necoils)
   if(coilnames_given)
      for i=1:necoils
         all_outputs = strvcat(all_outputs,names.coils(i,:));
      end
   elseif(isfield(vacuum_objs,'ecnames') & ~isempty(vacuum_objs.ecnames))
      all_outputs = strvcat(all_outputs,vacuum_objs.ecnames);
   else
      for i=1:necoils
         all_outputs = strvcat(all_outputs,['ecoil' int2str(i) ]);
      end
   end
end

if(coilnames_given)
   all_outputs = strvcat(all_outputs,names.coils(necoils+1:end,:));
elseif(~isempty(vacuum_objs.fcnames))
   all_outputs = strvcat(all_outputs,vacuum_objs.fcnames);
else
   if(~isempty(ecnturn))
      for i=1:length(ecnturn)
         all_outputs = strvcat(all_outputs,['ecoil' int2str(i) ]);
      end
   end
   all_outputs = strvcat(all_outputs,vacuum_objs.fcnames);
end
output_map.cc.block = 'control coil currents';
output_map.cc.ncc = ncc;
output_map.cc.start = 1;
output_map.cc.end = size(all_outputs,1);

if(isfield(names,'flux_loops') & ~isempty(names.flux_loops))
   all_outputs = strvcat(all_outputs,names.flux_loops);
elseif(~isempty(vacuum_objs.flnames))
   all_outputs = strvcat(all_outputs,vacuum_objs.flnames);
else
   all_outputs = strvcat(all_outputs,['psi_fl(' int2str(vacuum_objs.nfl) ')']);
end
output_map.fl.block = 'flux loops';
output_map.fl.nfl = nfl;
output_map.fl.start = output_map.cc.end+1;
output_map.fl.end = size(all_outputs,1);

if(isfield(names,'probes') & ~isempty(names.probes))
   all_outputs = strvcat(all_outputs,names.probes);
elseif(~isempty(vacuum_objs.bpnames))
   all_outputs = strvcat(all_outputs,vacuum_objs.bpnames);
else
   if(vacuum_objs.nbp>0)
     all_outputs = strvcat(all_outputs,['Bp(' int2str(vacuum_objs.nbp) ')']);
   end
end
output_map.bp.block = 'B-probes';
output_map.bp.nbp = nbp;
output_map.bp.start = output_map.fl.end+1;
output_map.bp.end = size(all_outputs,1);
last_end = output_map.bp.end;

output_map.rl.block = 'Rogowski loops';
output_map.rl.nrl = nrl;
output_map.rl.start = [];
output_map.rl.end = [];
if nrl > 0
   if(isfield(names,'rogowskis') & ~isempty(names.rogowskis))
      all_outputs = strvcat(all_outputs,names.rogowskis);
   elseif(isfield(vacuum_objs,'rlnames') && ~isempty(vacuum_objs.rlnames))
      all_outputs = strvcat(all_outputs,vacuum_objs.rlnames);
   else
      for k=1:nrl
         all_outputs = strvcat(all_outputs,['RL0' int2str(k)]);
      end
   end
   output_map.rl.start = last_end+1;
   output_map.rl.end = size(all_outputs,1);
   last_end = output_map.rl.end;
end

output_map.lv.block = 'loop voltages';
output_map.lv.nlv = nlv;
output_map.lv.start = [];
output_map.lv.end  = [];
if nlv > 0 
   if(isfield(names,'loop_voltages') & ~isempty(names.loop_voltages))
      all_outputs = strvcat(all_outputs,names.loop_voltages);
   elseif(isfield(vacuum_objs,'lvnames') && ~isempty(vacuum_objs.lvnames))
      all_outputs = strvcat(all_outputs,vacuum_objs.lvnames);
   else
      for k=1:nlv
         all_outputs = strvcat(all_outputs,['LV' int2str(k)]);
      end
   end
   output_map.lv.start = last_end+1;
   output_map.lv.end = size(all_outputs,1);
   last_end = output_map.lv.end;
end

output_map.mse.block = 'MSE';
output_map.mse.nmse = 2*nmse;
output_map.mse.start = [];
output_map.mse.end  = [];
if nmse > 0 
   if(isfield(names,'mse') & ~isempty(names.mse))
      all_outputs = strvcat(all_outputs,names.mse);
   elseif(isfield(vacuum_objs,'msenames') && ~isempty(vacuum_objs.msenames))
      for k=1:nmse
         all_outputs = strvcat(all_outputs,['Br' vacuum_objs.msenames(k,:)]);
      end
      for k=1:nmse
         all_outputs = strvcat(all_outputs,['Bz' vacuum_objs.msenames(k,:)]);
      end
   else
      for k=1:nmse
         all_outputs = strvcat(all_outputs,['BrMSE' int2str(k)]);
      end
      for k=1:nmse
         all_outputs = strvcat(all_outputs,['BzMSE' int2str(k)]);
      end
   end
   output_map.mse.start = last_end+1;
   output_map.mse.end = size(all_outputs,1);
   last_end = output_map.mse.end;
end

output_map.nl.block = 'netlist outputs';
output_map.nl.nnl = 0;
if(use_netlist & ~isempty(netlist_model.outputs))
   all_outputs = strvcat(all_outputs,netlist_model.outputs);
   output_map.nl.start = last_end+1;
   output_map.nl.end = size(all_outputs,1);
   last_end = output_map.nl.end;
else
   output_map.nl.start = [];
   output_map.nl.end = [];
end
if(~isempty(output_map.nl.end))
   output_map.nl.nnl = output_map.nl.end - output_map.nl.start + 1;
end

output_map.other.block = 'other outputs';
output_map.other.nother = 0;
if(~vacuum_model)
   all_outputs = strvcat(all_outputs,'R', 'Z');
   if(iplcirc)
      all_outputs = strvcat(all_outputs, 'Ip');
   end
   output_map.other.start = last_end+1;
   output_map.other.end = size(all_outputs,1);
else
   output_map.other.start = [];
   output_map.other.end = [];
end
if(~isempty(output_map.other.end))
   output_map.other.nother = output_map.other.end - output_map.other.start + 1;
end

% Input labels:
if(use_cccirc)
   reduced_outputs = ' ';
   idx = setdiff(unique(abs(cccirc)),0);	% exclude 0 (unconnected) coils
   for k=1:length(idx)
       idx1 = min(find(abs(cccirc)==idx(k)));
       reduced_outputs = strvcat(reduced_outputs,all_outputs(idx1,:));
   end
% idx = strmatch(reduced_outputs(end,:),all_outputs)
   idx = length(cccirc);
   for k=idx+1:size(all_outputs,1)
       reduced_outputs = strvcat(reduced_outputs, all_outputs(k,:));
   end
   reduced_outputs = reduced_outputs(2:end,:);	% remove beginning blank

   if(isfield(names,'inputs'))
      inputs = names.inputs;
   else
      inputs = ['V_' reduced_outputs(1,:)];
      for k=2:size(bmat,2)
         inputs = strvcat(inputs,['V_' reduced_outputs(k,:)]);
      end
   end
else
   reduced_outputs = ' ';

   inputs = netlist_model.inputs;
end

% State descriptions:

% ncx = nxx_cc - nvx - iplcirc;
states = '';
if(~isempty(ecdata))
   state_names = vacuum_objs.ecnames;
else
   state_names = '';
end
state_names = strvcat(state_names, vacuum_objs.fcnames);
for k=1:ncx
    idx = find(Pcc(:,k)==1);
    for j=1:length(idx)
       idx2 = find(Pcc(idx(j),:)~=0);
       if(length(idx2)==1)
          str = deblank(state_names(idx(1),:));
          break;
       end
    end
    str = [str ' current'];
    states = strvcat(states,str);
end
vessel_names= '';
 for k=1:vacuum_objs.nv
    vessel_names = strvcat(vessel_names,['vv' int2str(k)]);
 end
 for k=1:nvx
    idx = find(Pvv(:,k)~=0);
    mult = Pvv(idx(1),k);
    if(mult==1)
       multstr='';
    elseif(mult==-1)
       multstr='-';
    else
       multstr =[num2str(1/mult) '*'];
    end
    str = [multstr deblank(vessel_names(idx(1),:))];
    for j=2:length(idx)
       mult = Pvv(idx(j),k);
       if(mult==1)
          multstr='';
       elseif(mult==-1)
          multstr='-';
       else
          multstr =[num2str(1/mult) '*'];
       end
       str = [str '=' multstr deblank(vessel_names(idx(j),:))];
    end
    str = [str ' current'];
    states = strvcat(states,str);
 end
 if(iplcirc)
    states = strvcat(states, 'plasma current');
 end

% Construct output_signals object if any of basic "signal" objects exist:
if( isfield(vacuum_objs,'ecsignals') | isfield(vacuum_objs,'fcsignals') | ...
   	 isfield(vacuum_objs,'flsignals') | isfield(vacuum_objs,'bpsignals') | ...
    	 isfield(vacuum_objs,'rlsignals') | isfield(vacuum_objs,'lvsignals') )

   output_signals = cell(ncc+nfl+nbp+nrl+nlv,1);
   nout = 0;
   if( necoils>0 & isfield(vacuum_objs,'ecsignals') )
      output_signals(nout+[1:necoils],1:size(vacuum_objs.ecsignals,2)) = vacuum_objs.ecsignals;
   end
   nfc = ncc - necoils;

   nout = necoils;
   if( nfc>0 & isfield(vacuum_objs,'fcsignals') )
      output_signals(nout+[1:nfc],1:size(vacuum_objs.fcsignals,2)) = vacuum_objs.fcsignals;
   end

   nout = ncc;
   if( nfl>0 & isfield(vacuum_objs,'flsignals') )
      output_signals(nout+[1:nfl],1:size(vacuum_objs.flsignals,2)) = vacuum_objs.flsignals;
   end

   nout = ncc+nfl;
   if( nbp>0 & isfield(vacuum_objs,'bpsignals') )
      output_signals(nout+[1:nbp],1:size(vacuum_objs.bpsignals,2)) = vacuum_objs.bpsignals;
   end

   nout = ncc+nfl+nbp;
   if( nrl>0 & isfield(vacuum_objs,'rlsignals') )
      output_signals(nout+[1:nrl],1:size(vacuum_objs.rlsignals,2)) = vacuum_objs.rlsignals;
   end

   nout = ncc+nfl+nbp+nrl;
   if( nlv>0 & isfield(vacuum_objs,'lvsignals') )
      output_signals(nout+[1:nlv],1:size(vacuum_objs.lvsignals,2)) = vacuum_objs.lvsignals;
   end
   
   for j=1:size(all_outputs,1), for k=1:size(output_signals,2)
     if j>size(output_signals,1), output_signals{j,1} = 'none'; end
     if isempty(output_signals{j,k}), output_signals{j,k} = 'none'; end
   end, end
   
end

% Variable description structure to go into model structure:

description = struct( ...
'cc0', 'complete set of coil currents ', ...
'ncx', 'number of coils in reduced circuit', ...
'nvx', 'number of vessel elements in reduced circuit', ...
'amat', 'states [Ipf; Ivv] (Ip added if iplcirc=1)', ...
'bmat', 'inputs of control coils (V)', ...
'cmat', 'if iplcirc=0, outputs are [Icc(ncc);psi_fl(nfl);Bp(nbp);R;Z]. if iplcirc=1, Ip is added at end.', ...
'dmat', 'zeros(size(cmat,1),size(bmat,2));', ...
'fmat','matrix multiplying disturbance [delta_betap delta_li dbetap/dt dli/dt]', ...
'hmat','matrix multiplying disturbance [delta_betap delta_li dbetap/dt dli/dt]', ...
'Ieq', 'equilibrium current values for reduced system', ...
'yeq', 'equilibrium output values for reduced system', ...
'states','identity of state variables', ...
'outputs','identity of outputs of model', ...
'output_signals','names of signals that correspond to the output channels', ...
'output_map','start and end indices for blocks of diagnostics', ...
'mxx', 'conductor mutual inductance matrix (not including plasma motion effects)', ...
'rxx', 'conductor resistance matrix', ...
'xmatx', 'dpsi_conductors/dI_conductors caused by plasma motion', ...
'Pcc', 'maps from connected circuits to original unconnected circuits', ...
'Pxx', 'maps from full state vector with connected circs to unconn-circs', ...
'Resp', 'plasma resistance (Ohms)', ...
'Lp', 'plasma self inductance', ...
'Rextra','extra circuit resistance', ...
'Lextra','extra circuit inductance', ...
'lstar','Pxx^T*(mxx+xmatx)*Pxx', ...
'rzrig_data', 'data objects produced by rzrig', ...
'equil_data', 'data read from EFIT g-file', ...
'units', 'defines units used in data objects', ...
'flags','flags controlling how model was built');

if(isfield(tok_system,'netlist_model'))
   description.netlist_model = 'netlist model of system (see Engineering Physics Memo EPMmlw070120a)';
end

if(isfield(tok_system,'vmat'))
   description.vmat = 'mapping of PF voltages into PF coil circuit';
end

if(vacuum_objs.imks)
  if(vacuum_objs.iterminal)
     other_units = ...
	 'currents in Amps, flux in Wb, field in T, lengths in meters';
  else
     other_units = ...
     'currents in Amp-turns, flux in Wb/turn, field in T, lengths in meters';
  end
else
  if(vacuum_objs.iterminal)
     other_units = ...
	 'currents in MAmps, flux in Wb, field in T, lengths in meters';
  else
     other_units = ...
     'currents in MA-turns, flux in Wb/turn, field in T, lengths in meters';
  end
end
units = struct( ...
'rzrig_data', 'units defined in rzrig_data.units', ...
'equil_data', 'units defined in equil_data.gdef', ...
'others',other_units);

%save tok_system

if(vacuum_model)		% create temporary phony data for plasma
   cc0=0; Resp=0; Lp=0; rzrig_data=0; equil_data=0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put all model quantities into data structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flags = struct( ...
'imks', vacuum_objs.imks, ...
'iterminal', vacuum_objs.iterminal, ... 
'ichooseq', ichooseq, ...
'iplcirc', iplcirc, ...
'vacuum_model',vacuum_model, ...
'irzresp_dynamic', irzresp_dynamic, ...
'irzresp_output', irzresp_output);

tok_system = struct( ...	% build the structure
'tokamak', tokamak, ...
'config_name', 'default', ...
'cc0',cc0, ...
'vc0',vc0, ...
'ncx', ncx, ...
'nvx', nvx, ...
'amat', amat, ...
'bmat', bmat, ...
'cmat', cmat_allPF, ...
'dmat', dmat_allPF, ...
'fmat', fmat, ...
'Ieq', Ieq, ...
'yeq', yeq, ...
'inputs',inputs, ...
'states',states, ...
'outputs',all_outputs, ...
'output_signals',[], ...
'output_map',output_map, ...
'disturbance_inputs',disturbances, ...
'mxx', mxx, ...
'rxx', rxx, ...
'xmatx', xmatx, ...
'Pcc', Pcc, ...
'Pvv', Pvv, ...
'Pxx', Pxx, ...
'Resp',Resp, ...
'Lp',Lp, ...
'Rextra',Rextra, ...
'Lextra',Lextra, ...
'lstar',lstar, ...
'gspert_data', [], ...
'equil_data', equil_data, ...
'flags',flags, ...
'description',description, ...
'units',units);

if(exist('Rckt','var') & Rckt~=0)
   tok_system.Rckt=Rckt;
end
if(exist('Lckt','var') & Lckt~=0)
   tok_system.Lckt=Lckt;
end

if(exist('output_signals','var'))
   tok_system.output_signals = output_signals;
end
if(exist('config_name')==1)		% Allow re-definition of config_name...
   tok_system.config_name = config_name;
elseif(isfield(vacuum_objs,'config_name'))	% ...otherwise use vac objects.
   tok_system.config_name = vacuum_objs.config_name;
end
if(isempty(vc0))
   tok_system = rmfield(tok_system,'vc0');
end

if(vacuum_model)	% then remove the phony data from the structure	
   tok_system = rmfield(tok_system,'cc0');
   tok_system = rmfield(tok_system,'Resp');
   tok_system = rmfield(tok_system,'Lp');
   tok_system = rmfield(tok_system,'rzrig_data');
   tok_system = rmfield(tok_system,'gspert_data');
   tok_system = rmfield(tok_system,'equil_data');
   tok_system = rmfield(tok_system,'fmat');
   tok_system = rmfield(tok_system,'hmat');
   tok_system = rmfield(tok_system,'disturbance_inputs');
else
   if(irzresp_dynamic==5 | irzresp_output==5)
      tok_system.gspert_data = response;
   else
      tok_system = rmfield(tok_system,'gspert_data');
   end
end
if(exist('netlist_model'))	% for testing
   tok_system.netlist_model = netlist_model;
end
               
if(nargout>1)
   dbg_objs = struct( ...
      'lstar',lstar, ...
      'L_star',L_star, ...
      'Lextra',Lextra);
   if(exist('lstari'))
      dbg_objs.lstari=lstari;
   end
      if(~vacuum_model)
         dbg_objs.xmats = xmats;
         dbg_objs.xmatx = xmatx; 
         dbg_objs.dbg_rzrig = dbg_rzrig;
         dbg_objs.dfldis = dfldis;
         dbg_objs.dbpdis = dbpdis;
      end
end
% save temp.mat	% TESTING

if isfield(build_inputs,'scale_ip')
   tok_system = scale_equil_response(tok_system,scaleip);
end
