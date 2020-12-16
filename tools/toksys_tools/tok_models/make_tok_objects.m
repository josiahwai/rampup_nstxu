function make_tok_objects(make_tok_inputs,procedures,dummy)
 


%  USAGE:  make_tok_objects(make_tok_inputs,procedures)
%
%  PURPOSE: Script to calculate vacuum data objects (e.g. mutuals, Green fns)
%	for a tokamak system model. The "minimal" set of objects which must be
%	defined is VV, FC.
%
%  INPUTS (see TokSys_Users_Guide.pdf):
%   make_tok_inputs = data structure defining how to make data objects. 
%     This structure must contain:
%       tokamak = string defining tokamak being constructed, used in file names
%       config_name  = string label identifying machine configuration (Warning: 
%	 		non-alphanumeric chars in config_name may cause problems)
%       datadir = directory where files defining tokamak are located
%	vvdata_file  = Filename for vvdata data (must be of form *.data)
%		         (actually just has to have 4 characters after ".", and
%		         must have a name distinct from all the other *.data)
%	fcdata_file  = Fcoil data file (name must be of same form as vvdata)
%	fcnturn_file = file with PF turns 
%	nminvv       = # of conductors in *min* dimension of vacuum vessel elts
%		         (can be scalar or vector with length=#vessel elts)
%	nminfc       = # of conductors in *min* dimension of fcoils
%		         (can be scalar or vector with length=#fcoils)
%	Vessel resistance data, either:
%	  etav       = VV resistivity vector, uOhm-m  
%			OR
%	  vvres_file = file with vessel element (terminal) resistances (Ohms)
%	PF coil resistance data, either:
%	  etaf       = Fcoil resistivity vector, uOhm-m  
%				(NOTE: for SC coils, set etaf=1e-3*eta(Cu))
%			OR
%	  fcres_file = file with coil element (terminal) resistances (Ohms)
%	TD coil resistance data, either:
%	  etat       = TD coil resistivity vector, uOhm-m  
%			OR
%	  tdres_file = file with coil element (terminal) resistances (Ohms)
%     Optionally, it may contain:
%	ecdata_file  = Ecoil data file (name must be of same form as vvdata)
%	ecturn_file  = ecoil turn groupings
%	vvfrac_file  = file with vacuum vessel fractions of currents
%	fldata_file  = Flux loop data file 
%	bpdata_file  = Bprobe data file 
%	msedata_file  = MSE loci file 
%       rldata_file  = Rogowski loop data file
%       lvdata_file  = Loop voltage data file
%	limdata_file = file containing limiter definition data (npts x 2)
%	nr,nz        = # of grid elem. in r,z dir. (if~=0, must also enter
%			         values for rgmin,rgmax,zgmin,zgmax - see below)
%	rgmin,rgmax  = Min,max in major radial dimension on plasma grid [m]
%	zgmin,zgmax  = Min,max in vertical dimension on plasma grid [m]
%	E-coil resistance data (required if ecdata_file used), either:
%	  etae       = Ecoil resistivity vector, uOhm-m  
%			OR
%	  ecres_file = file with coil element (terminal) resistances (Ohms)
%	fcnames_file = file containing fcoil names
%	vvnames_file = file containing vessel element names
%	flnames_file = file containing flux loop names
%	bpnames_file = file containing Bprobe names
%	msenames_file = file containing names of MSE channels that view points in msedata
%       rlnames_file = file containing Rogowski loop names
%       lvnames_file = file containing loop voltage names
%	nminec       = # of conductors in *min* dimension of ecoils. Required if
%		ecdata_file is used. (scalar or vector of length=#ecoils)
%	nmingg       = scalar # of conductors in *min* dimension of plasma grid
%			Required if plasma grid inputs used.
%       plot_tok_geo_fn = string defining script to execute to plot the geometry
%	            of the device (e.g. 'plot_kstar_geo', 'plot_east_geo'...
%                   set to 'generic' (default) to use default to plot X-section
%       ecsignals_file = E coil signal names file
%       fcsignals_file = F coil signal names file
%       flsignals_file = flux loop signal names file
%       bpsignals_file = B probe signal names file
%       rlsignals_file = Rogowski loop signal names file
%       lvsignals_file = loop voltage signal names file
%       procedures   = array of strings defining procedures to operate on the
%			         data after generation (optional, default = [])
%
%  OUTPUTS: save files
%    <tokamak>_obj_<config_name> = contains structure tok_data_struct, units= terminal, MKS
%				(in file name, '/' replaced by '-', ' ' by '_', '\' by '')
 
%  RESTRICTIONS:
%	Note that the "minimal" set of objects which must be defined in 
%	order to save a save-set is VV, FC, and grid data.
%	fcnturn_file name must be such that when load it, will have fcnturn = 
%	  vector of turns (so e.g. fcnturn.mat which contains vector fcnturn,
%	  or fcnturn.dat ascii file containing vector of data).
%	Assumes all turns of E-coil are rectangular cross-section.
 
%  METHOD:  
%	Loads data files, defines 
%	geometry with define_<tokamak>_data, then uses mutind_fine to 
%	calculate mutuals. bgreens_fine to calculate fields at B-probes.
%
%  WRITTEN BY:  Mike Walker	ON	3/31/04
%	based on multiple versions of make*objects.m files by Dave Humphreys
%
%  VERSION @(#)make_tok_objects.m	1.35 02/23/15
%
%  MODIFICATIONS:
%     2007-04-19  JAL   (default e,fcnames), add: config_name,ccnturn,ccnames..
%                       now only 2 arguments with plot function in input, also
%                       input stored in out
%     2008-01-21  NWE   now loads rogowski & LV data and name files &
%                       calculates LV mutuals
%     2009-05-11  ASW   Adding MSE diagnostic
%     2009-10-27  ASW   More info about names such as ptd, pcs names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_testing = 0;

% TO DO:
% 1) Provide method for updating previously computed values (some calculations 
% take a long time to get accurate).
% 2) Remove Input: dummy, Remove default plot_tok_geo Remove output: fcnturn


   if nargin<=0
     disp('ERROR make_tok_objects: needs at least one input: make_tok_inputs')
     help make_tok_objects   
     return;
   elseif nargin==1
     procedures= [];
   elseif nargin==3
     disp('WARNING: make_tok_objects: now has only 2 inputs (struc,procedures)')
     wait('         Old call used 3 inputs (struc,plot,proceeds) Convert ?')
     procedures= dummy;
   end


   nprocs = size(procedures,1);

   if(~isfield(make_tok_inputs,'tokamak'))
      wait('ERROR make_tok_objects: input structure must contain "tokamak"')
      return;
   end
   if(~isfield(make_tok_inputs,'config_name'))
      wait('ERROR make_tok_objects: input structure must contain "config_name"')
      return;
   end
   if(~isfield(make_tok_inputs,'datadir'))
      wait('ERROR make_tok_objects: input structure must contain "datadir"')
      return;
   end
   if(~isfield(make_tok_inputs,'vvdata_file'))
      wait('ERROR make_tok_objects: input structure must contain "vvdata_file"')
      return;
   end
   if(~isfield(make_tok_inputs,'fcdata_file'))
      wait('ERROR make_tok_objects: input structure must contain "fcdata_file"')
      return;
   end
   if(~isfield(make_tok_inputs,'fcnturn_file'))
     wait('ERROR make_tok_objects: input structure must contain "fcnturn_file"')
      return;
   end
   if(~isfield(make_tok_inputs,'nminvv'))
      wait('ERROR make_tok_objects: input structure must contain "nminvv"')
      return;
   end
   if(~isfield(make_tok_inputs,'nminfc'))
      wait('ERROR make_tok_objects: input structure must contain "nminfc"')
      return;
   end
   if(isfield(make_tok_inputs,'ecdata_file'))
      if(~isempty(make_tok_inputs.ecdata_file) ...
	   		& ~isfield(make_tok_inputs,'nminec'))
         wait('ERROR make_tok_objects: input structure must contain "nminec"')
         return;
      end
      if(~isfield(make_tok_inputs,'etae') & ...
			 ~isfield(make_tok_inputs,'ecres_file'))
         wait('ERROR make_tok_objects: either etae or ecres_file must be input')
         return;
      end
   end
   if(isfield(make_tok_inputs,'etav') & isfield(make_tok_inputs,'vvres_file'))
      wait('ERROR make_tok_objects: only one of etav and vvres_file allowed')
      return;
   end
   if(~isfield(make_tok_inputs,'etav') & ~isfield(make_tok_inputs,'vvres_file'))
      wait('ERROR make_tok_objects: either etav or vvres_file must be input')
      return;
   end
   if(isfield(make_tok_inputs,'etaf') & isfield(make_tok_inputs,'fcres_file'))
      wait('ERROR make_tok_objects: only one of etaf and fcres_file allowed')
      return;
   end
   if(~isfield(make_tok_inputs,'etaf') & ~isfield(make_tok_inputs,'fcres_file'))
      wait('ERROR make_tok_objects: either etaf or fcres_file must be input')
      return;
   end
   if(isfield(make_tok_inputs,'etae') & isfield(make_tok_inputs,'ecres_file'))
      wait('ERROR make_tok_objects: only one of etae and ecres_file allowed')
      return;
   end

   struct_to_ws(make_tok_inputs);
   if(~exist('flsignals_file'))    
      flsignals_file = '';
   end
   if(~exist('bpsignals_file'))    
      bpsignals_file = '';
   end
   if(~exist('fcsignals_file'))    
      fcsignals_file = '';
   end
   if(~exist('ecsignals_file'))    
      ecsignals_file = '';
   end
   if(~exist('rlsignals_file'))    
      rlsignals_file = '';
   end
   if(~exist('lvsignals_file'))    
      lvsignals_file = '';
   end
   if(~exist('tdsignals_file'))    
      tdsignals_file = '';
   end
   if(~exist('slsignals_file'))    
      slsignals_file = '';
   end
   if(~exist('msesignals_file'))    
      msesignals_file = '';
   end
   if(~exist('ecdata_file'))  
      ecdata_file = '';
   end
   if(~exist('fldata_file'))    
      fldata_file = '';
   end
   if(~exist('bpdata_file'))  
      bpdata_file = '';
   end
   if(~exist('msedata_file'))  
      msedata_file = '';
   end
   if(~exist('limdata_file')) 
      limdata_file = '';
   end
   if(~exist('nr'))    
      nr = 0;
   end
   if(~exist('nz'))    
      nz = 0;
   end
   if(~exist('fcnames_file'))    
      fcnames_file = '';
   end
   if(~exist('ecnames_file'))    
      ecnames_file = '';
   end
   if(~exist('vvnames_file'))    
      vvnames_file = '';
   end
   if(~exist('flnames_file'))    
      flnames_file = '';
   end
   if(~exist('bpnames_file'))    
      bpnames_file = '';
   end
   if(~exist('msenames_file'))    
      msenames_file = '';
   end
   if(~exist('fcid_file'))
      fcid_file='';
   end
   if(~exist('vvid_file'))
      vvid_file='';
   end
   if(~exist('vvfrac_file'))
      vvfrac_file='';
   end
   if(~exist('ecturn_file'))
      ecturn_file='';
   end
   if(~exist('fcturn_file'))
      fcturn_file='';
   end
   if(~exist('turnfc_file'))
      turnfc_file='';
   end
   if(~exist('rlnames_file'))    
      rlnames_file = '';
   end
   if(~exist('rldata_file'))    
      rldata_file = '';
   end
   if(~exist('lvnames_file'))    
      lvnames_file = '';
   end
   if(~exist('lvdata_file'))    
      lvdata_file = '';
   end
   if(~exist('fcres_file'))    
      fcres_file = '';
   end
   if(~exist('ecres_file'))    
      ecres_file = '';
   end
   
   if exist('plot_tok_geo_fn')~=1 
      plot_tok_geo_fn='generic';
   end

% Prelims and Constants:
   mu0 = 0.4*pi;
   twopi = 2*pi;

% Derived Values:
 if nr&nz  %only do grid calcs if nr & nz ~=0...
   rg = linspace(rgmin,rgmax,nr)';
   zg = linspace(zgmin,zgmax,nz)';
   dr = rg(2)-rg(1);
   dz = zg(2)-zg(1);
   rgg = ones(nz,1)*rg';
   zgg = zg*ones(1,nr);
  %Define grid data object ggdata:
   ngg = length(zgg(:));
   ggdata = [zgg(:)'; rgg(:)'; dz*ones(1,ngg); dr*ones(1,ngg); zeros(2,ngg)];
 end
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Geometry Files, plot confirmation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nvv=0; nff=0; nex=0; nfl=0; nbp=0; nee=0; nlim=0; ncc=0; nmse=0;

if limdata_file
   load([datadir,limdata_file])
   nlim=size(limdata,2);
 else
   limdata=[];
end
if vvdata_file
   load([datadir,vvdata_file])
   [s1,nvv]=size(vvdata);
   if(s1~=6)
      wait('WARNING make_tok_objects: 1st dim. vvdata should be 6')
   end
   if(vvfrac_file), load([datadir,vvfrac_file]); end;
 else 
   vvdata = [];
end
if(exist('etav')==1 & max(size(etav)) > 1)
   if max(size(etav)) ~= size(vvdata,2)
      wait('ERROR make_tok_objects: length(etav) must = # of vessel elements')
      return;
   end
end
if fcdata_file
   load([datadir,fcdata_file]);
   [s1,nff] = size(fcdata);
   if(s1~=6)
      wait('WARNING make_tok_objects: 1st dim. fcdata should be 6')
   end
   load([datadir,fcnturn_file]);
   if(length(fcnturn)~=nff)
      wait('ERROR make_tok_objects: length(fcnturn) must match 2nd dim fcdata')
      return;
   end
 else
   fcdata = [];
end
if vvnames_file
  load([datadir,vvnames_file])
else
  vvnames = [];
end
if fcnames_file %jal04apr2007
   load([datadir,fcnames_file]);            % matfile with fcnames variable
elseif exist('nff') && ~isempty('nff')      % apply default PF1, PF2, ....
    if nff <= 0
      fcnames=[];
    elseif nff == 1
       fcnames= 'PF';
    else 
      fcnames= [];
      for ii=1:nff
	fcnames= strvcat(fcnames,['PF' int2str(ii)]);
      end
    end
else
    fcnames=[];
end
if fcturn_file
   load([datadir,fcturn_file]);
else
   fcturn = [];
end
if turnfc_file
   load([datadir,turnfc_file]);
else
   turnfc = [];
end

ecturn = [];	% only make ecturn non-empty if ecoil turns are "grouped"
if ecdata_file
   load([datadir,ecdata_file])
   [s1,s2] = size(ecdata);
   if(s1~=5)
      wait('ERROR make_tok_objects: ecdata must have 5 rows')
      return;
   end
%   if(s1>s2)
%     disp(['NOTE make_tok_objects: ecdata(' int2str(s1) ',' int2str(s2) ')'])
%     disp('WARNING make_tok_objects: Transposing to get 1 element per column')
%     ecdata = ecdata'; 
%     wait(['  ecdata transposed to:   ecdata(' int2str(s2) ',' int2str(s1) ')'])
%   end
   nex=size(ecdata,2); 
   nee = length(unique(ecdata(5,:)));
   if(ecturn_file)
      temp = load([datadir,ecturn_file],'ecturn');
      if(isfield(temp,'ecturn'))
         ecturn = temp.ecturn;
      else
         disp('ERROR make_tok_objects: ecturn file must must contain the variable ecturn')
         return
      end
      [s1,s2]=size(ecturn);
      if(s1<s2), ecturn = ecturn'; end;
   end
else
   ecdata = [];
end  
if ecnames_file  %jal04apr2007
   load([datadir,ecnames_file]);            % matfile with ecnames variable
   idx = unique(ecdata(5,:));
   if(length(idx)~=size(ecnames,1))
     wait(['ERROR make_tok_objects: size(ecnames) conflicts with ' ...
	 'number of ecoils in ecdata'])
     return;
   end
elseif exist('nee') && ~isempty('nee')      % apply default OH1,OH2,OH3...
    if nee <= 0
      ecnames=[];
    elseif nee == 1
       ecnames= 'OH';
    else 
       ecnames= [];
  for ii=1:nee
	 ecnames= strvcat(ecnames,['OH' int2str(ii)]);
       end
    end
else
    ecnames=[];
end
% flux loops
fldata = [];
if fldata_file
   load([datadir,fldata_file])
   nfl = size(fldata,2);
   if size(fldata,1) < 6
     fldata(6,1) = 0; % Need to be 6 rows
   end
end
if(exist('flzr_file','var'))
   wait('WARNING make_tok_objects: flzr file is no longer supported. Please convert to fldata file');
% If an flzr_file (obsolete) is specified, the dzfl and drfl must also
% be specified.
   clear flzr fldata
   load([datadir,flzr_file])
   if exist('fldata'), flzr=fldata; end
   nfl=size(flzr,2); 
   fldata = [flzr(1:2,:); zeros(4,size(flzr,2))];
%   fldata1 = [flzr(1:2,:); dzfl*ones(1,size(flzr,2)); ...
%			   drfl*ones(1,size(flzr,2)); ...
%			   zeros(4,size(flzr,2))];
end 
if flnames_file
   load([datadir,flnames_file]);   %matfile with flnames variable
   snames = size(flnames,1);
   sdata = size(fldata,2);
   if(snames~=sdata)
      disp('WARNING: size(flnames) does not match size(fldata)')
      if(snames<sdata)
          disp('Filling undefined names with defaults ')
          for k=snames+1:sdata
             flnames = strvcat(flnames,['fl' int2str(k)]);
          end
      else
          disp('Truncating names array')
          flnames = flnames(1:sdata,:);
      end
      wait;
   end
 else 
   flnames = [];
end
% Magnetic probes
if bpdata_file
   load([datadir,bpdata_file])
   nbp = size(bpdata,2); 
 else
   bpdata = [];
end
if bpnames_file
   load([datadir,bpnames_file]);   %matfile with bpnames variable
   snames = size(bpnames,1);
   sdata = size(bpdata,2);
   if(snames~=sdata)
      wait('WARNING: size(bpnames) does not match size(bpdata)')
      if(snames<sdata)
          disp('Filling undefined names with blanks')
          for k=snames+1:sdata
             bpnames = strvcat(bpnames,'   ');
          end
      else
          disp('Truncating names array')
          bpnames = bpnames(1:sdata,:);
      end
      wait;
   end
 else 
   bpnames = [];
end
% MSE diagnostics
if msedata_file
   load([datadir,msedata_file])
   nmse = size(msedata,2); 
 else
   msedata = [];
end
if msenames_file
   load([datadir,msenames_file]);
else 
   msenames = [];
end
% Rogowski loops
if rldata_file
   load([datadir,rldata_file]) 
   nrl = size(rldata,2); 
   if(size(rldata,1) ~= nee+nff+nvv+1) 
      wait(['ERROR make_tok_objects: size(rldata,1) != total'...
               ' # of conductors + Ip']);
     return;
   end
 else
   rldata = [];
   nrl = 0;
end
if rlnames_file
   load([datadir,rlnames_file]);  
   snames = size(rlnames,1);
   sdata = size(rldata,2);
   if(snames~=sdata)
      wait('WARNING: size(rlnames) does not match size(rldata)')
      if(snames<sdata)
          disp('Filling undefined names with blanks')
          for k=snames+1:sdata
             rlnames = strvcat(rlnames,'   ');
          end
      else
          disp('Truncating names array')
          rlnames = rlnames(1:sdata,:);
      end
      wait;
   end
 else 
   rlnames = [];
end
% loop voltage
if lvdata_file
   load([datadir,lvdata_file])    
   nlv=size(lvdata,2); 
   if size(lvdata,1) < 6
     lvdata(6,1) = 0; % Need to be 6 rows
   end
else
   lvdata = [];
end
if lvnames_file
   load([datadir,lvnames_file]);   
   snames = size(lvnames,1);
   sdata = size(lvdata,2);
   if(snames~=sdata)
      wait('WARNING: size(lvnames) does not match size(lvdata)')
      if(snames<sdata)
          disp('Filling undefined names with blanks')
          for k=snames+1:sdata
             lvnames = strvcat(lvnames,'   ');
          end
      else
          disp('Truncating names array')
          lvnames = lvnames(1:sdata,:);
      end
      wait;
   end
 else 
   lvnames = [];
   nlv = 0;
end
if(~isempty(fcres_file))
   temp = load([datadir,fcres_file],'fcres');
   if(isfield(temp,'fcres'))
      fcres = temp.fcres;
   else
      disp('ERROR make_tok_objects: fcres file must must contain the variable fcres')
      return
   end
end
if(~isempty(ecres_file))
   temp = load([datadir,ecres_file],'ecres');
   if(isfield(temp,'ecres'))
      ecres = temp.ecres;
   else
      disp('ERROR make_tok_objects: ecres file must must contain the variable ecres')
      return
   end
end

% If "signals" file is a mat file, then assume it is already cell array.
% If not, assume it is (text) string array.  String arrays will only
% work if there is only one column of names.

if(exist('flsignals_file') & ~isempty(flsignals_file))
   load([datadir,flsignals_file]);
   if(strcmp(flsignals_file(end-3:end),'.dat'))
      flsignals = strarr_to_cellarr(flsignals);
   end
   ssignals = size(flsignals,1);
   sdata = size(fldata,2);
   if(ssignals~=sdata)
      disp('WARNING: size(flsignals) does not match size(fldata)')
      if(ssignals<sdata)
          wait('Filling undefined signals with defaults ')
          for k=ssignals+1:sdata
             flsignals = strvcat(flsignals,'none');
          end
      else
          wait('Truncating signals array')
          flsignals = flsignals(1:sdata,:);
      end
   end
else
   flsignals=[];
end
if(exist('bpsignals_file') & ~isempty(bpsignals_file))
   load([datadir,bpsignals_file]);
   if(strcmp(bpsignals_file(end-3:end),'.dat'))
      bpsignals = strarr_to_cellarr(bpsignals);
   end
   ssignals = size(bpsignals,1);
   sdata = size(bpdata,2);
   if(ssignals~=sdata)
      disp('WARNING: size(bpsignals) does not match size(bpdata)')
      if(ssignals<sdata)
          wait('Filling undefined signals with defaults ')
          for k=ssignals+1:sdata
             bpsignals = strvcat(bpsignals,'none');
          end
      else
          wait('Truncating signals array')
          bpsignals = bpsignals(1:sdata,:);
      end
   end
else
   bpsignals=[];
end
if(exist('fcsignals_file') & ~isempty(fcsignals_file))
   load([datadir,fcsignals_file]);
   if(strcmp(fcsignals_file(end-3:end),'.dat'))
      fcsignals = strarr_to_cellarr(fcsignals);
   end
   ssignals = size(fcsignals,1);
   sdata = size(fcdata,2);
   if(ssignals~=sdata)
      disp('WARNING: size(fcsignals) does not match size(fcdata)')
      if(ssignals<sdata)
          wait('Filling undefined signals with defaults ')
          for k=ssignals+1:sdata
             fcsignals = strvcat(fcsignals,'none');
          end
      else
          wait('Truncating signals array')
          fcsignals = fcsignals(1:sdata,:);
      end
   end
else
   fcsignals=[];
end
if(exist('ecsignals_file') & ~isempty(ecsignals_file))
   load([datadir,ecsignals_file]);
   if(strcmp(ecsignals_file(end-3:end),'.dat'))
      ecsignals = strarr_to_cellarr(ecsignals);
   end
   ssignals = size(ecsignals,1);
   sdata = length(unique(ecdata(5,:)));
   if(ssignals~=sdata)
      disp('WARNING: size(ecsignals) does not match size(ecdata)')
      if(ssignals<sdata)
          wait('Filling undefined signals with defaults ')
          for k=ssignals+1:sdata
             ecsignals = strvcat(ecsignals,'none');
          end
      else
          wait('Truncating signals array')
          ecsignals = ecsignals(1:sdata,:);
      end
   end
else
   ecsignals=[];
end
if(exist('rlsignals_file') & ~isempty(rlsignals_file))
   load([datadir,rlsignals_file]);
   if(strcmp(rlsignals_file(end-3:end),'.dat'))
      rlsignals = strarr_to_cellarr(rlsignals);
   end
   ssignals = size(rlsignals,1);
   sdata = size(rldata,2);
   if(ssignals~=sdata)
      disp('WARNING: size(rlsignals) does not match size(rldata)')
      if(ssignals<sdata)
          wait('Filling undefined signals with defaults ')
          for k=ssignals+1:sdata
             rlsignals = strvcat(rlsignals,'none');
          end
      else
          wait('Truncating signals array')
          rlsignals = rlsignals(1:sdata,:);
      end
   end
else
   rlsignals=[];
end
if(exist('lvsignals_file') & ~isempty(lvsignals_file))
   load([datadir,lvsignals_file]);
   if(strcmp(lvsignals_file(end-3:end),'.dat'))
      lvsignals = strarr_to_cellarr(lvsignals);
   end
   ssignals = size(lvsignals,1);
   sdata = size(lvdata,2);
   if(ssignals~=sdata)
      disp('WARNING: size(lvsignals) does not match size(lvdata)')
      if(ssignals<sdata)
          wait('Filling undefined signals with defaults ')
          for k=ssignals+1:sdata
             lvsignals = strvcat(lvsignals,'none');
          end
      else
          wait('Truncating signals array')
          lvsignals = lvsignals(1:sdata,:);
      end
   end
else
   lvsignals=[];
end


if(fcid_file)
   load([datadir,fcid_file])
else
   fcid=[1:nff];
end
if(vvid_file)
   load([datadir,vvid_file])
else
   vvid=[];
end


%  3-D object initialization
   tddata = []; tdwireD = []; tdnturn = []; rest = []; tdnames = []; tdsignals = [];
   if(~exist('tddata_file'))  
      tddata_file = '';
   end
   if(~exist('tdwireD_file'))  
      tdwireD_file = '';
   end
   if(~exist('tdnturn_file'))  
      tdnturn_file = '';
   end
   if(~exist('tdres_file'))  
      tdres_file = '';
   end
   if(~exist('tdnames_file'))    
      tdnames_file = '';
   end
   if(~exist('tdsignals_file'))    
      tdsignals_file = '';
   end  
   if ~isempty(tdnturn_file)
     load([datadir,tdnturn_file]);
   end
   if ~isempty(tddata_file)
     load([datadir,tddata_file]);
     ntd = length(tddata);
     if ~exist('tdnturn','var'), tdnturn = ones(ntd,1); end
     Ptt = diag(tdnturn);
   else
     ntd = 0;
   end
   if ~isempty(tdwireD_file)
     load([datadir,tdwireD_file]);
   end
   if ~isempty(tdres_file)
     load([datadir,tdres_file]);
   end
   if ~isempty(tdnames_file)
     load([datadir,tdnames_file]);
   end
   if ~isempty(tdsignals_file)
     load([datadir,tdsignals_file]);
   end
   if ntd>0
     if ~isempty(bpdata) & size(bpdata,1)<5
       disp('ERROR make_tok_objects: A 5:th row with toroidal angle in degrees is required for bpdata');
       return
     end
     if ~isempty(msedata) & size(msedata,1)<5
       disp('ERROR make_tok_objects: A 5:th row with toroidal angle in degrees is required for msedata');
       return
     end
   end
   
   sldata = []; slnames = []; slsignals = [];
   if(~exist('sldata_file'))  
      sldata_file = '';
   end
   if(~exist('slnames_file'))    
      slnames_file = '';
   end
   if(~exist('slsignals_file'))    
      slsignals_file = '';
   end
   if ~isempty(sldata_file)
     load([datadir,sldata_file]);
     nsl = length(sldata);
   else
     nsl = 0;
   end
   if ~isempty(slnames_file)
     load([datadir,slnames_file]);
   end
   if ~isempty(slsignals_file)
     load([datadir,slsignals_file]);
   end

% ===================================================
% Confirmation plot:
  figure
  iplteq=0; ipltflx=0; ipltfl=0; ipltbp=0;
  ilabelfc=0; ilabelvv=0; ilabelfl=0; ilabelbp=0;
  if exist('psizr'), iplteq=1; ipltflx=1; end
  if(~isempty(fldata))
    ipltfl=1;   %turn on plotting of FL's
    ilabelfl=0; %turn off labeling of FL indices
  end
  if bpdata_file
    ipltbp=1;   %turn on plotting of BP's
    ilabelbp=0; %turn off labeling of BP indices  
  end
  if(~isempty(plot_tok_geo_fn))
     if(strcmp(plot_tok_geo_fn,'generic'))
        tok_data_struct.tokamak = tokamak;
        tok_data_struct.vvdata = vvdata;
        tok_data_struct.fcdata = fcdata;
        tok_data_struct.ecdata = ecdata;
        tok_data_struct.limdata = limdata;
        tok_data_struct.fldata = fldata;	
        tok_data_struct.bpdata = bpdata;	
        plot_tok_geo(tok_data_struct);
     else
        eval(plot_tok_geo_fn)
     end
  end
  
  if ~exist('ilabel3dplot','var'), ilabel3dplot = 0; end
% Confirmation plot for 3-D objects:
  if ~isempty(tddata) | ~isempty(sldata)
    figure
    set(gcf,'OuterPosition', [0 0 1200 1000])
    if(~isempty(tddata))
      plot_fils(tddata,'r');
      if ilabel3dplot & exist('tdnames','var')
	for j=1:size(tddata,2)
         h =  text(mean(tddata{j}(:,1)),mean(tddata{j}(:,2)),mean(tddata{j}(:,3)),tdnames(j,:));
	 set(h,'color','r');
	end
      end
    end
    if ~isempty(sldata)
      hold on
      plot_fils(sldata,'b');
      if ilabel3dplot & exist('slnames','var')
	for j=1:size(sldata,2)
          h = text(mean(sldata{j}(:,1)),mean(sldata{j}(:,2)),mean(sldata{j}(:,3)),slnames(j,:));
	  set(h,'color','b');
	end
      end
    end
    if exist('bpdata','var') & size(bpdata,1)>=5
      hold on
      for j=1:size(bpdata,2)
        Rbp = bpdata(2,j)+bpdata(4,j)*[-1 1]/2*cos(bpdata(3,j)*pi/180);
        xbp = Rbp*cos(bpdata(5,j)*pi/180);
	ybp = Rbp*sin(bpdata(5,j)*pi/180);
	zbp = bpdata(1,j)+bpdata(4,j)*[-1 1]/2*sin(bpdata(3,j)*pi/180);
        plot3(xbp,ybp,zbp,'k','linew',2)
	if ilabel3dplot & exist('bpnames','var')
          h = text(mean(xbp),mean(ybp),mean(zbp),bpnames(j,:));
	  set(h,'color','k');
	end
      end
    end
    if exist('msedata','var') & size(msedata,1)>=5
      hold on
      for j=1:size(msedata,2)
        Rms = msedata(2,j)+msedata(4,j)*[-1 1]/2*cos(msedata(3,j)*pi/180);
        xms = Rms*cos(msedata(5,j)*pi/180);
	yms = Rms*sin(msedata(5,j)*pi/180);
	zms = msedata(1,j)+msedata(4,j)*[-1 1]/2*sin(msedata(3,j)*pi/180);
        plot3(xms,yms,zms,'m','linew',2)
	if ilabel3dplot & exist('msenames','var')
          h = text(mean(xms),mean(yms),mean(zms),msenames(j,:));
	  set(h,'color','m');
	end
      end
    end
    title([upper(tokamak) ': 3D objects'])
  end
  drawnow
  
% ===================================================

%  wait('Paused: CR to continue...')

% Construct Ecoil sub-coil geometry from ecdata 
if nex
  ecdata1 = [];
  ecnturn = zeros(nee,1);	% force to be a column
  for i=1:nee
     idx1=find(ecdata(5,:)==i); ecnturn(i)=length(idx1);
     ecdata1 = [ecdata1 ecdata(1:4,idx1)];
  end
  ecdata1 = [ecdata1; zeros(2,sum(ecnturn))];
else
  ecnturn = [];
end

if(do_testing)
% Save test data for use in evaluating convergence of magnetics calculations.
% (Use function test_convergence.m.)

   test_data_file = [tokamak '_test_data.mat'];
   eval(['save ' test_data_file]);
   fprintf('test data saved in %s\n',test_data_file)
   wait
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Resistances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VV:
  if(exist('etav')==1)
     if(size(etav,2)>1), etav = etav'; end
     resv = twopi*vvdata(2,:)'.*etav./(vvdata(3,:).*vvdata(4,:))';
  elseif(exist('vvres_file')==1)  
     resv = vvres_file;
  else
     resv = [];
  end

% FC:
  if(exist('etaf')==1)
     if(size(etaf,2)>1), etaf = etaf'; end
     resf = twopi*fcdata(2,:)'.*etaf./(fcdata(3,:).*fcdata(4,:))';
  elseif(exist('fcres')==1)  
     resf = fcres;
  else
     resf = [];
  end

% EC:  
if nex
  if(exist('etae')==1)
     if(size(etae,2)>1), etae = etae'; end
     rese = [];
     for i=1:nee
       idx1=find(ecdata(5,:)==i); 
       rese1=sum(twopi*ecdata(2,idx1).*etae(i)./(ecdata(3,idx1).*ecdata(4,idx1)))';
       rese = [rese; rese1/ecnturn(i)^2];
     end
  elseif(exist('ecres')==1)  
     rese = ecres;
  else
     rese = [];
  end
end

% TD:
  if(exist('etat')==1)
     for j = 1:ntd
       rest(j,1) = 0; % Initialize this resistance to 0
       nfils = size(tddata{j},1); % Number of filaments
       Atd = pi*tdwireD(j)^2/4; % Cross sectional area
       for k = 1:nfils
         rest(j) = rest(j)+etat/Atd*norm(tddata{j}(k,1:3)-tddata{j}(k,4:6)); % micro-Ohms
       end
     end
  elseif(exist('tdres')==1)  
     rest = tdres; % Ohms
  else
     rest = [];
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate VV-VV mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nvv
  disp('Calculating MVV...')
  tic
% need finer discretization for non-rectangle shape elements in VV
  if length(nminvv)==1 %jal19apr2007
     nminvv= nminvv*ones(size(vvdata,2),1);
  elseif length(nminvv) ~= size(vvdata,2)
     wait('%ERROR: make_tok_objects: length(nminvv) must be 1 or length(vv)')
  end
  nmintgt = 2*nminvv;
  idx = union(find(vvdata(5,:)~=0),find(vvdata(6,:)~=0));
  nmintgt(idx) = 4*nminvv(idx);
  mvv = mutind_fine(vvdata,vvdata,nminvv,nmintgt,1);
  disp('MVV required: ')
  toc
  mvv = 0.5*(mvv+mvv');
else
  mvv = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate FF-VV mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nvv&nff
  disp('Calculating MFV...')
  tic
  mfv = mutind_fine(fcdata,vvdata,nminfc,nmintgt,1);
  disp('MFV required: ')
  toc
else
  mfv = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Ecoil-VV mutuals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nee&nvv
  disp('Calculating MEV...')
  tic
  mve1 = mutind_fine(vvdata,ecdata1,nminvv,nminec,1);

  mev = [];
  n2=0;
  for i=1:nee
     mev = [mev; mean(mve1(:,n2+1:n2+ecnturn(i))')];
     n2 = n2 + ecnturn(i);
  end
  disp('MEV required: ')
  toc
  clear mve1
else
  mev = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Plasma-VV mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ngg&nvv
  disp('Calculating MPV, GBR2V, BGZ2V ...')
  tic
  mpv = mutind_fine(vvdata,ggdata,nminvv,2*nmingg,1)';
  grddata = [zgg(:)'; rgg(:)'; zeros(1,ngg); zeros(1,ngg); zeros(1,ngg)];
  [gbr2v,gbz2v] = bgreens_fine(vvdata,grddata,nminvv,1);
  disp('MPV, GBR2V, GBZ2V required: ')
  toc
else
  mpv = []; gbr2v = []; gbz2v = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate FL-VV mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nfl&nvv
  disp('Calculating MLV...')
  tic
  mlv = mutind_fine(vvdata,fldata,nminvv,1,1)';
  disp('MLV required: ')
  toc
else
  mlv = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate BP-VV Green functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nbp&nvv
  disp('Calculating GBV...')
  tic
  [br,bz,gbv] = bgreens_fine(vvdata,bpdata,nminvv,1);
  disp('GBV required: ')
  toc
else
  gbv = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate MSE-VV Green functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nmse&nvv
  disp('Calculating GMSEV...')
  tic
  [gmsebrv,gmsebzv] = bgreens_fine(vvdata,[msedata;zeros(1,nmse);zeros(1,nmse)+0.02;zeros(1,nmse)],nminvv,1);
  disp('GMSEV required: ')
  toc
else
  gmsebrv = []; gmsebzv = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate FF-FF mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nff
  disp('Calculating MFF...')
  tic
% need finer discretization for non-rectangle shape elements
  if length(nminfc)==1 %jal19apr2007
     nminfc= nminfc*ones(size(fcdata,2),1);
  elseif length(nminfc) ~= size(fcdata,2)
     wait('%ERROR: make_tok_objects: length(nminfc) must be 1 or leng(fc)')
  end
  nmintgt = 2*nminfc;
  idx = union(find(fcdata(5,:)~=0),find(fcdata(6,:)~=0));
  nmintgt(idx) = 4*nminfc(idx);
  mff = mutind_fine(fcdata,fcdata,nminfc,nmintgt,1);
  disp('MFF required: ')
  toc
  mff = 0.5*(mff+mff');
else
  mff = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Ecoil-FF mutuals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nee&nff
  disp('Calculating MEF...')
  tic
  mfe1 = mutind_fine(fcdata,ecdata1,nminfc,nminec,1);

  mef = [];
  n2=0;
  for i=1:nee
     mef = [mef; mean(mfe1(:,n2+1:n2+ecnturn(i))')];
     n2 = n2 + ecnturn(i);
  end

  disp('MEF required: ')
  toc
  clear mfe1
else
  mef = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Plasma-FF mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ngg&nff
  disp('Calculating MPF, GBR2F, GBZ2F ...')
  tic
  mpf = mutind_fine(fcdata,ggdata,nminfc,2*nmingg,1)';
  grddata = [zgg(:)'; rgg(:)'; zeros(1,ngg); zeros(1,ngg); zeros(1,ngg)];
  [gbr2f,gbz2f] = bgreens_fine(fcdata,grddata,nminfc,1);
  disp('MPF, GBR2F, GBZ2F required: ')
  toc
else
  mpf = []; gbr2f = []; gbz2f = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Plasma-FL mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ngg&nfl
  disp('Calculating MPL...')
  tic
  mpl = mutind_fine(ggdata,fldata,nmingg,1,1);
  disp('MPL required: ')
  toc
else
  mpl = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate FL-FF mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nfl&nff   
  disp('Calculating MLF...')
  tic
  mlf = mutind_fine(fcdata,fldata,nminfc,1,1)';
  disp('MLF required: ')
  toc
else
  mlf = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate BP-FF Green functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nbp&nff
  disp('Calculating GBF...')
  tic
  [br,bz,gbf] = bgreens_fine(fcdata,bpdata,nminfc,1);
  disp('GBF required: ')
  toc
else
  gbf = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate MSE-FF Green functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nmse&nff
  disp('Calculating GMSEF...')
  tic
  [gmsebrf,gmsebzf] = bgreens_fine(fcdata,[msedata;zeros(1,nmse);zeros(1,nmse)+0.02;zeros(1,nmse)],nminfc,1);
  disp('GMSEF required: ')
  toc
else
  gmsebrf = []; gmsebzf = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Plasma-BP Green functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nbp&ngg
  disp('Calculating GPB...')
  tic
  [br,bz,tmp] = bgreens_fine(ggdata,bpdata,nmingg,1);
  gpb = tmp';
  disp('GPB required: ')
  toc
else
  gpb = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Plasma-MSE Green functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nmse&ngg
  disp('Calculating GPMSE...')
  tic
  [gmsebrp,gmsebzp] = bgreens_fine(ggdata,[msedata;zeros(1,nmse);zeros(1,nmse)+0.02;zeros(1,nmse)],nmingg,1);
  gpb = tmp';
  disp('GPMSE required: ')
  toc
else
  gmsebrp = []; gmsebzp = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate plasma-Ecoil mutuals 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nee&ngg
  disp('Calculating MPE, GBR2E, GBZ2E ...')
  tic
  mpe1 = mutind_fine(ecdata1,ggdata,1,2*nmingg,1)';

  mpe = [];
  n2=0;
  for i=1:nee
     mpe = [mpe mean(mpe1(:,n2+1:n2+ecnturn(i))')'];
     n2 = n2 + ecnturn(i);
  end

  grddata = [zgg(:)'; rgg(:)'; zeros(1,ngg); zeros(1,ngg); zeros(1,ngg)];
  [gbr2e1,gbz2e1] = bgreens_fine(ecdata1,grddata,1,1);
  gbr2e = []; gbz2e = [];
  n2=0;
  for i=1:nee
     gbr2e = [gbr2e mean(gbr2e1(:,n2+1:n2+ecnturn(i))')'];
     gbz2e = [gbz2e mean(gbz2e1(:,n2+1:n2+ecnturn(i))')'];
     n2 = n2 + ecnturn(i);
  end

  disp('MPE, GBR2E, GBZ2E required: ')
  toc
  clear mpe1
else
  mpe = []; gbr2e = []; gbz2e = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Ecoil self inductance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nee
  disp('Calculating MEE...')
  tic
  mee1 = mutind_fine(ecdata1,ecdata1,1,2*nminec,1);

  mee = []; n2=0;
  for i=1:nee
     mee = [mee mean(mee1(:,n2+1:n2+ecnturn(i))')'];
     n2 = n2 + ecnturn(i);
  end
  mee1 = mee; mee=[]; n2=0;
  for i=1:nee
     if size(mee1,1) > 1
        mee = [mee; mean(mee1(n2+1:n2+ecnturn(i),:))];
     else
        mee = [mee; mean(mee1(n2+1:n2+ecnturn(i)))];
     end
     n2 = n2 + ecnturn(i);
  end

  disp('MEE required: ')
  toc
  clear mee1
else
  mee = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate FL-Ecoil mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nfl&nee   
  disp('Calculating MLE...')
  tic
  mel = mutind_fine(ecdata1,fldata,1,1,1);

  mle = []; n2=0;
  for i=1:nee
     mle = [mle  mean(mel(n2+1:n2+ecnturn(i),:))'];
     n2 = n2 + ecnturn(i);
  end

  disp('MLE required: ')
  clear mel
  toc
else
  mle = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate BP-Ecoil Green functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nbp&nee
  disp('Calculating GBE...')
  tic
  [br,bz,gbe1] = bgreens_fine(ecdata1,bpdata,1,1);

  gbe = []; n2=0;
  for i=1:nee
     gbe = [gbe mean(gbe1(:,n2+1:n2+ecnturn(i))')'];
     n2 = n2 + ecnturn(i);
  end

  disp('GBE required: ')
  toc
else
  gbe = [];
end
clear gbe1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate MSE-Ecoil Green functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nmse&nee
  disp('Calculating GMSEE...')
  tic
  [gmsebre1,gmsebze1] = bgreens_fine(ecdata1,[msedata;zeros(1,nmse);zeros(1,nmse)+0.02;zeros(1,nmse)],1,1);

  gmsebre = []; gmsebze = []; n2=0;
  for i=1:nee
     gmsebre = [gmsebre mean(gmsebre1(:,n2+1:n2+ecnturn(i))')'];
     gmsebze = [gmsebze mean(gmsebze1(:,n2+1:n2+ecnturn(i))')'];
     n2 = n2 + ecnturn(i);
  end

  disp('GMSEE required: ')
  toc
else 
  gmsebre = []; gmsebze = [];
end
clear gmsebre1 gmsebze1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate plasma-plasma mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ngg
  disp('Calculating MPP...')
  tic
  [mpp,gbr2p,gbz2p] = calc_mpp(rgg,zgg,dr,dz,nz,nr);
  disp('MPP required: ')
  toc
else
  mpp = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate LV-Ecoil mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nlv&nee   
  disp('Calculating MHE...')
  tic
  meh = mutind_fine(ecdata1,lvdata,1,1,1);

  mhe = []; n2=0;
  for i=1:nee
     mhe = [mhe  mean(meh(n2+1:n2+ecnturn(i),:))'];
     n2 = n2 + ecnturn(i);
  end

  disp('MHE required: ')
  clear meh
  toc
else
  mhe = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate LV-FF mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nlv&nff   
  disp('Calculating MHF...')
  tic
  mhf = mutind_fine(fcdata,lvdata,nminfc,1,1)';
  disp('MHF required: ')
  toc
else
  mhf = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Plasma-LV mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ngg&nlv
  disp('Calculating MPH...')
  tic
  mph = mutind_fine(ggdata,lvdata,nmingg,1,1);
  disp('MPH required: ')
  toc
else
  mph = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate LV-VV mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nlv&nvv
  disp('Calculating MHV...')
  tic
  mhv = mutind_fine(vvdata,lvdata,nminvv,1,1)';
  disp('MHV required: ')
  toc
else
  mhv = [];
end

% Calculate mutuals for 3D objects

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Plasma-TD mutuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ngg&ntd
  disp('Calculating MPT ...') 
  tic
  mpt = mut_fine_fil(ggdata,fils2filcs(tddata),nmingg,1);
  disp('MPT required: ')
  toc
else
  mpt = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Ecoil-TD mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ntd&nee
  disp('Calculating MET...')
  tic
  met1 = mut_fine_fil(ecdata1,fils2filcs(tddata),nminec,1);
  met = [];
  n2=0;
  for i=1:nee
     met = [met; mean(met1(n2+1:n2+ecnturn(i),:))];
     n2 = n2 + ecnturn(i);
  end
  disp('MET required: ')
  toc
  clear met1
else
  met = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate FF-TD mutuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ntd&nff
  disp('Calculating MFT...')
  tic
  mft = mut_fine_fil(fcdata,fils2filcs(tddata),nminfc,1);
  disp('MFT required: ')
  toc
else
  mft = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate VV-TD mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ntd&nvv
  disp('Calculating MVT...')
  tic
  mvt = mut_fine_fil(vvdata,fils2filcs(tddata),nminvv,1);
  disp('MVT required: ')
  toc
else
  mvt = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate TD-TD mutuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ntd
  disp('Calculating MTT...')
  tic
  mtt =  mut_filcs2filcs(fils2filcs(tddata),fils2filcs(tddata),num2cell(tdwireD/2));
  disp('MTT required: ')
  toc
  mtt = 0.5*(mtt+mtt');
else
  mtt = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate FL-TD mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ntd&nfl
  disp('Calculating MLT...')
  tic
  mlt = mut_fine_fil(fldata,fils2filcs(tddata),1,1);
  disp('MLT required: ')
  toc
else
  mlt = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate LV-TD mutuals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ntd&nlv
  disp('Calculating MHT...')
  tic
  mht = mut_fine_fil(lvdata,fils2filcs(tddata),1,1);
  disp('MHT required: ')
  toc
else
  mht = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate SL-TD mutuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nsl&ntd
  disp('Calculating MST...')
  tic
  mst = mut_filcs2filcs(fils2filcs(sldata),fils2filcs(tddata));
  disp('MST required: ')
  toc
else
  mst = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate BP-TD Green functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nbp&ntd
  disp('Calculating GBT...')
  tic
  [gr, gz, gt] = green_paths2pts([bpdata(2,:)' bpdata(1,:)' pi/180*bpdata(5,:)'],fils2filcs(tddata)');
  gbt = -(gr.*cos(pi*bpdata(3,:)'*ones(1,ntd)/180) + gz.*sin(pi*bpdata(3,:)'*ones(1,ntd)/180));
  disp('GBT required: ')
  toc
else
  gbt = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate MSE-TD Green functions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nmse&ntd
  disp('Calculating GMSEBRT GMSEBZT...')
  tic
  [gmsebrt,gmsebzt] = green_paths2pts([msedata(2,:)' msedata(1,:)' pi/180*msedata(5,:)'],fils2filcs(tddata)');
  disp('GMSEBRT GMSEBZT required: ')
  toc
else
  gmsebrt = []; gmsebzt = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Plasma-SL mutuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ngg&nsl
  disp('Calculating MPS...')
  tic
  mps = mut_fine_fil(ggdata,fils2filcs(sldata),nmingg,1);
  disp('MPS required: ')
  toc
else
  mps = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate SL-Ecoil mutuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nsl&nee
  disp('Calculating MSE...(mutuals between Saddle loops and Ecoils)')
  tic
  mse1 = mut_fine_fil(ecdata1,fils2filcs(sldata),nminec,1)';
  mse = [];
  n2=0;
  for i=1:nee
     mse = [mse mean(mse1(:,n2+1:n2+ecnturn(i))')'];
     n2 = n2 + ecnturn(i);
  end
  disp('MSE required: ')
  toc
  clear mse1
else
  mse = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate SL-FF mutuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nsl&nff
  disp('Calculating MSF...')
  tic
  msf = mut_fine_fil(fcdata,fils2filcs(sldata),nminfc,1)';
  disp('MSF required: ')
  toc
else
  msf = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate SL-VV mutuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nsl&nvv
  disp('Calculating MSV...')
  tic
  msv = mut_fine_fil(vvdata,fils2filcs(sldata),nminvv,1)';
  disp('MSV required: ')
  toc
else
  msv = [];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge and Save Objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nff&nee
  ncc = nff + nee;
  resc = [rese; resf];
  mcc = [[mee mef]; [mef' mff]];
  mcv = [mev; mfv];
  mct = [met; mft];
  mlc = [mle mlf];
  mhc = [mhe mhf];
  msc = [mse msf];
  gbc = [gbe gbf];
  gmsebrc = [gmsebre gmsebrf];
  gmsebzc = [gmsebze gmsebzf];
  mpc = [mpe mpf];
  gbr2c = [gbr2e gbr2f];
  gbz2c = [gbz2e gbz2f];
  ccnames= strvcat(ecnames,fcnames); %jal18apr2007
% Up to this point, ecnturn can represent either the number of individual
% turns, or the number of "grouped" turns as in nstx.  If ecturn exists,
% use it to convert to actual physical turns - needed for correct
% conversion between lumped data objects and terminal data objects.
  if(~isempty(ecturn))
    for k=1:nee
      idx = find(ecdata(5,:)==k);
      ecnturn(k) = sum(ecturn(idx));
    end
  end
  ccnturn= [ecnturn; fcnturn];  %jal18apr2007
else 
 if nff
  ncc = nff;
  resc = [resf];
  mcc = [mff];
  mcv = [mfv];
  mct = [mft];
  mlc = [mlf];
  mhc = [mhf];
  msc = [msf];
  gbc = [gbf]; 
  gmsebrc = [gmsebrf]; 
  gmsebzc = [gmsebzf]; 
  mpc = [mpf]; 
  gbr2c = [gbr2f];
  gbz2c = [gbz2f];
  ccnames= fcnames;  %jal18apr2007
  ccnturn= fcnturn;  %jal18apr2007
 end
end

if(~exist('vvfrac'))
   vvfrac=[];
end

% Custom modifications of any of the calculated objects:
for k=1:nprocs
   disp(['NOTE: Evaluating procedure: ' procedures(k,:)])
   eval(procedures(k,:))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert everything to MKS and terminal units and save.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(exist('etav')==1)
   resv = resv * 1e-6;
end
mvv  = mvv * 1e-6;
mpv  = mpv * 1e-6;
mlv  = mlv * 1e-6;
mhv  = mhv * 1e-6;
gbv  = gbv * 1e-6;
mpl  = mpl * 1e-6;
mph  = mph * 1e-6;
gpb  = gpb * 1e-6;
gmsebrv  = gmsebrv * 1e-6;
gmsebzv  = gmsebzv * 1e-6;
gmsebrp  = gmsebrp * 1e-6;
gmsebzp  = gmsebzp * 1e-6;
mpp  = mpp * 1e-6;
gbr2p = gbr2p * 1e-6;
gbz2p = gbz2p * 1e-6;

turnmat = diag(ccnturn); %jal19apr2007

if(exist('etae')==1)
   resc(1:nee) = resc(1:nee) * 1e-6; 
   resc(1:nee) = resc(1:nee).*ccnturn(1:nee).^2;
end
if(exist('etaf')==1)
   resc(nee+1:end) = resc(nee+1:end) * 1e-6; 
   resc(nee+1:end) = resc(nee+1:end).*ccnturn(nee+1:end).^2;
end
mcv = mcv * 1e-6;  mcv = turnmat*mcv;
mcc = mcc * 1e-6;  mcc = turnmat*mcc*turnmat;
mpc = mpc * 1e-6;  mpc = mpc*turnmat;
if(~isempty(msc))
  msc = msc*turnmat;
end
if(~isempty(mlc))
   mlc  = mlc * 1e-6;  mlc = mlc*turnmat;
end
if(~isempty(mhc))
   mhc  = mhc * 1e-6;  mhc = mhc*turnmat;
end
if(~isempty(gbc))
   gbc  = gbc * 1e-6;  gbc = gbc*turnmat;
end;
if(~isempty(gmsebrc))
   gmsebrc  = gmsebrc * 1e-6;  gmsebrc = gmsebrc*turnmat;
   gmsebzc  = gmsebzc * 1e-6;  gmsebzc = gmsebzc*turnmat;
end
if(~isempty(gbr2c))
   gbr2c = gbr2c * 1e-6; gbr2c = gbr2c*turnmat;
   gbz2c = gbz2c * 1e-6; gbz2c = gbz2c*turnmat;
end
if(~isempty(gbr2v))
   gbr2v = gbr2v * 1e-6;
   gbz2v = gbz2v * 1e-6;
end
if ntd
  if(exist('etat')==1)
     rest = rest * 1e-6.*tdnturn; % tdwireD is effective cross-section diameter for *one* turn
  end
  mpt = mpt*Ptt;
  mct = turnmat*mct*Ptt;
  mvt = mvt*Ptt;
  mtt = Ptt*mtt*Ptt;
  if(~isempty(mlt))
    mlt = mlt*Ptt;
  end
  if(~isempty(mht))
    mht = mht*Ptt;
  end
  if(~isempty(mst))
    mst = mst*Ptt;
  end
  if(~isempty(gbt))
    gbt = gbt*Ptt;
  end
  if(~isempty(gmsebrt))
    gmsebrt = gmsebrt*Ptt;
  end
  if(~isempty(gmsebzt))
    gmsebzt = gmsebzt*Ptt;
  end
end

   datafile = [tokamak '_obj_tmp_mks'];
   savestr = ['save ' datafile  ...
        ' make_tok_inputs' ...
   	' ecdata fcdata vvdata tddata fldata sldata bpdata msedata rldata lvdata limdata ' ...
	' tokamak fcnames ecnames vvnames tdnames flnames slnames bpnames msenames rlnames lvnames ccnames ' ...
	' rg zg rgg zgg nr nz fcid vvid ecturn fcturn turnfc ' ...
	' nvv ncc ntd nfl nsl nbp nrl nlv nmse resv resc rest fcnturn ecnturn ccnturn tdnturn vvfrac' ...
        ' mvv mcv mpv mlv mhv gbv mcc mpc mlc mhc gbc mpl mph gpb mtt mct mvt mpt mlt mht mst gbt mps msv msc' ...
	' mpp gbr2p gbz2p gbr2c gbz2c gbr2v gbz2v nmingg nminfc nminvv' ...
	' gmsebrv gmsebzv gmsebrc gmsebzc gmsebrp gmsebzp gmsebrt gmsebzt' ...
	' flsignals bpsignals fcsignals ecsignals tdsignals rlsignals lvsignals slsignals'];
   if(exist('nminec'))
      savestr = [savestr ' nminec'];
   end
   if(exist('config_name'))
      savestr = [savestr ' config_name'];
   end
   eval(savestr)

   imks=1; iterminal=1;
   tok_data_struct=make_tok_data_struct(datafile,imks,iterminal); 
   savestr = ['save ' tokamak '_obj'];
   temp = config_name;
   idx = findstr(' ',temp); temp(idx) = '_';
   idx = findstr('/',temp); temp(idx) = '-';
   idx = findstr('\',temp); temp(idx) = '';
   if(exist('config_name'))
      savestr = [savestr '_' temp];
   end
   savestr = [savestr ' tok_data_struct'];
   eval(savestr)
   unix(['rm ' datafile '.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   disp('All done.')
