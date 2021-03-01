function tok_data_struct = load_tok_objects(tok,config_name,efit_grid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  SYNTAX:  tok_data_struct = load_tok_objects(tok,config_name, grid)
%
%  PURPOSE:  Load vacuum data objects for selected device configuration.
%		Provides general access to load_<device>_object.m scripts
%
%  INPUT:
%	tok = string defining tokamak device or requests for help, eg help or help d3d
%       config_name = string defining configuration name
%	efit_grid = string identifying the EFIT grid resolution
%
%  OUTPUT:
%       tok_data_struct = vacuum objects data structure
%
%  CALLING SEQUENCE:
% 	tok_data_struct = load_tok_objects('help')
% 	tok_data_struct = load_tok_objects('help d3d')
%	
%	Most Current TokSys Models as of 20121005 (need to confirm these...)
%	   tok_data_struct = load_tok_objects('d3d','current')
%
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Matthew J. Lanctot 
%
%  MODIFICATIONS: 
%     	MJL 2012/10/05	Created.Provide access to all load_<device>_object.m scripts
%	DAH 2014/7/11 	Updated EAST current to 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  @(#)load_tok_objects.m	1.2 10/05/12

if nargin==0
   config_name = 'help';
end
tok_data_struct = [];	% in case of error in execution

gatools_root = getenv('GATOOLS_ROOT');

switch tok
   case 'help'
      disp('Usage:')
      disp('   tok_data_struct = load_tok_objects(<tok>,<config>,<efit_grid>);')
      disp('   Put <> entries in single quotes')
      disp('Available tok devices:')
      disp('   D3D, NSTX, EAST, KSTAR, ITER, FDF')
      disp('Help:')
      disp('   Call with string tok=help (in single quotes)')
      disp('   For device specific help (say d3d), call with string tok=help d3d  (in single quotes)')

   case 'help d3d'
      disp('Usage:')
      disp('   tok_data_struct = load_tok_objects(<tok>,<config>,<efit_grid>);')
      disp('   (Put <config> in single quotes...)')
      disp('Available configurations:')
      disp('   after_ADP:')
      disp('   with_ADP:')
      disp('Available efit grids:')
      disp('	with_ADP: 3365, 6565')
      disp('	after_ADP: 3333, 3365,')
      disp('                   6565, 6565_biggrd, biggrd_6565_fine,')
      disp('		   RTEFIT, RTEFIT6565,')
      disp('		   129')	   

   case 'd3d'
     switch config_name
	case 'after_ADP'
   	   grids={'3333','3365','6565',...
   		  '6565_biggrd', 'biggrd_6565_fine',...
   		  'RTEFIT', 'RTEFIT6565', '129'};
   	   if sum(strcmpi(efit_grid,grids)) < 1      
   	     disp(['Requested efit_grid is not available for config_name ' config_name])
   	     disp('Returning 129 grid')
   	     efit_grid = '129';
   	   end   
   	   eval(['load ' gatools_root '/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_' efit_grid '.mat;'])
   
	case 'with_ADP'
   	   grids={'3365','6565'};
   	   if sum(strcmpi(efit_grid,grids)) < 1      
   	     disp(['Requested efit_grid is not available for config_name ' config_name])
   	     disp('Returning 6565 grid')
   	     efit_grid = '6565';
   	   end   
   	   eval(['load ' gatools_root '/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_' efit_grid '.mat;'])

	case 'current'
	   eval(['load ' gatools_root '/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_129.mat;'])
	   
	otherwise
	   wait('ERROR load_tok_objects: invalid configuration name for d3d')
	   tds=load_tok_objects('help d3d');	 
     end % d3d

   case 'help nstx'
      disp('Usage:')
      disp('   tok_data_struct = load_nstx_objects(<config>);')
      disp('   (Put <config> in single quotes...)')
      disp('Available configurations:')
      disp('   02feb: 20020207: 28 Bp probes, 66 Flux loops, 19 VF coils, 62 VV cells')
      disp('   05apr: 20050402: 70 Bp probes, 61 Flux loops, 25 VF coils, 68 VV cells')
      disp('Available efit grids:')
      disp('	02feb: 6565')
      disp('	05apr: 3333, 6565')

   case 'nstx'
      efit_grid='';
      switch config_name
      	case '02feb'
	   grids={'6565'};
	   if sum(strcmpi(efit_grid,grids)) < 1      
	     disp(['Requested efit_grid is not available for config_name ' config_name])
	     disp('Returning 6565 grid')
	     efit_grid = '6565';
	   end   
	   eval(['load ' gatools_root '/tokamaks/nstx/make/02feb/nstx_obj_02feb_' efit_grid '.mat;'])
	
	case '05apr'
	   grids={'3333','6565'};
	   if sum(strcmpi(efit_grid,grids)) < 1      
	     disp(['Requested efit_grid is not available for config_name ' config_name])
	     disp('Returning 6565 grid')
	     efit_grid = '6565';
	   end    
	   eval(['load ' gatools_root '/tokamaks/nstx/make/05apr/nstx_obj_05apr_' efit_grid '.mat;'])
	   
	case 'current'
	   eval(['load ' gatools_root '/tokamaks/nstx/make/05apr/nstx_obj_05apr_6565.mat;'])
	
	otherwise
	   wait('ERROR load_tok_objects: invalid configuration name for nstx')	
	   tds=load_tok_objects('help nstx');
      end %nstx   

   case 'help kstar'
      disp('Usage:')
      disp('   tok_data_struct = load_kstar_objects(<config>);')
      disp('   (Put <config> in single quotes...)')
      disp('Available configurations:')
      disp('   2008: 2008: 74 Bp probes, 5 flux loops, 18 VF coils, 150 VV cells ')
      disp('   2009: 2009: 115 Bp probes, 45 flux loops, 18 VF coils, 166 VV cells')
      disp('   2010: 2010: 128 Bp probes, 45 flux loops, 18 VF coils, 184 VV cells ')
      disp('   2010_Sabbagh: 2010 geometry but with Sabbaghs EFIT model (differences?)')
      disp('   3333: 2010: 2010 geometry but on 3333 EFIT grid and 136 Bp probes, 174 VV cells')
      disp('   6565: 2010: 3333 but on 6565 EFIT grid')

   case 'kstar'
      switch config_name
         case '2008'
            eval(['load ' gatools_root '/tokamaks/kstar/make/kstar_obj_mks_struct_' config_name '.mat;'])

         case '2009'
            eval(['load ' gatools_root '/tokamaks/kstar/make/kstar_obj_mks_struct_' config_name '.mat;'])  
            
         case '2010'
            eval(['load ' gatools_root '/tokamaks/kstar/make/kstar_obj_mks_struct_' config_name '.mat;'])

         case '2010_Sabbagh'
            eval(['load ' gatools_root '/tokamaks/kstar/make/kstar_obj_mks_struct_' config_name '.mat;'])
            
         case '3333'
            eval(['load ' gatools_root '/tokamaks/kstar/make/kstar_obj_mks_struct_' config_name '.mat;']) 

         case '6565'
            eval(['load ' gatools_root '/tokamaks/kstar/make/kstar_obj_mks_struct_' config_name '.mat;'])
	    
	 case 'current'
	    eval(['load ' gatools_root '/tokamaks/kstar/make/kstar_obj_mks_struct_6565.mat;'])

         otherwise
            wait('ERROR load_tok_objects: invalid configuration name for kstar')
	    tds=load_tok_objects('help kstar');
      end %kstar

    case 'help east'
            disp('Usage:')
            disp('   tok_data_struct = load_east_objects(<config>);')
            disp('   (Put <config> in single quotes...)')
            disp('Available configurations:')
            disp('   2008: 2008: 38 Bp probes, 35 flux loops, 16 VF coils, 90 VV cells ')
            disp('   2010: 2010: Geometry unchanged but .mat structure has additional null entries')
            disp('   2012: 2012: Geometry unchanged but .mat structure has revised setup (no null entries)')
            disp('Available efit grids:')
            disp('   2008: 3333, 6565') 
            disp('   2010: 3333, 6565')   
            disp('   2012: 3333, 6565')  
	    disp('   2014: 3333, 6565')
      
   case 'east'
      switch config_name
         case '2008'
            grids={'3333','6565'};
            if sum(strcmpi(efit_grid,grids)) < 1      
              disp(['Requested efit_grid is not available for config_name ' config_name])
              disp('Returning 6565 grid')
              efit_grid = '6565';
            end   
            eval(['load ' gatools_root '/tokamaks/east/make/east_obj_2008_' efit_grid '.mat;'])

         case '2010'
            grids={'3333','6565'};
            if sum(strcmpi(efit_grid,grids)) < 1      
              disp(['Requested efit_grid is not available for config_name ' config_name])
              disp('Returning 6565 grid')
              efit_grid = '6565';
            end   
            eval(['load ' gatools_root '/tokamaks/east/make/east_obj_2010_' efit_grid '.mat;'])
           
         case '2012'
            grids={'3333','6565'};
            if sum(strcmpi(efit_grid,grids)) < 1      
              disp(['Requested efit_grid is not available for config_name ' config_name])
              disp('Returning 6565 grid')
              efit_grid = '6565';
            end   
            eval(['load ' gatools_root '/tokamaks/east/make/east_obj_2012_' efit_grid '.mat;']) 
	    
        case '2014'
            grids={'3333','6565'};
            if sum(strcmpi(efit_grid,grids)) < 1      
              disp(['Requested efit_grid is not available for config_name ' config_name])
              disp('Returning 6565 grid')
              efit_grid = '6565';
            end   
            eval(['load ' gatools_root '/tokamaks/east/make/east_obj_2014_' efit_grid '.mat;']) 
	    	    
	 case 'current'
	    eval(['load ' gatools_root '/tokamaks/east/make/east_obj_2014_6565.mat;'])

         otherwise
            wait('ERROR load_tok_objects: invalid configuration name for east')
	    tds=load_tok_objects('help east');
      end   %east

   case 'help iter'
      disp('Usage:')
      disp('   tok_data_struct = load_iter_objects(<config>);')
      disp('   (Put <config> in single quotes...)')
      disp('Available configurations:')   
      disp('	standard: 20080125: New geometry. 48 Bp probes, 24 flux loops, 12 VF coils, 113 VV cells')
      disp('	standard: 20080125: New geometry. Has double wall and plate inside first wall below midplane on LFS')	   
      disp('	hairpin:  20080130: Standard with 6 additional VF coils on LFS just outside first wall')      
      disp('	vvopt1:   20080131: Standard with 4 plates sections insdie first wall: 2 above & below midplane on HFS and LFS')
      disp('	cc_ins_vv:20080215: Standard with 7 additional VF coils: 3 inside first wall and 4 outside second wall')
      disp('	elmvv02:  20080331: Standard with 3 ELM coils (2 VF coils each) located inside first wall at, above, below the midplane')
      disp('	elmvv03:  20080419: Standard with 2 additional VF coils at Z=+3.45 and Z=-2.63')
      disp('	elmvac03: 20080419: Similar to elvv03 but coil positions are shifted')
      disp('	2010:	  20100625: New geometry. Subtle changes in vertical coil locations, limiter.  2 VF coils inside 1st wall.')
      disp('	2010v3p3: 20100819: 2010 with slight change in vertical position of PS6')
      disp('	2011:	  20110504: 2010v3p3 with 14 additional Bp probes and modified sensor geometry')	    
      disp('	vs3a:	  20110505: Partial 3rd wall for vertical stability')
   
   case 'iter'
      switch config_name
         case '2011'
            eval(['load ' gatools_root '/tokamaks/iter/make/2011/iter_obj_2011;'])

         case '2010v3p3'
            load([gatools_root '/tokamaks/iter/make/2010v3p3/iter_obj_2010v3p3.mat'])

         case '2010'
            eval(['load ' gatools_root '/tokamaks/iter/make/2010/iter_obj_2010;'])

         case 'cc_ins_vv'
              eval(['load ' gatools_root '/tokamaks/iter/make/cc_ins_vv/iter_obj_mks_struct;'])

         case 'elmvac03'
              eval(['load ' gatools_root '/tokamaks/iter/make/elmvac03/iter_obj_mks_struct;'])

         case 'elmvv02'
              eval(['load ' gatools_root '/tokamaks/iter/make/elmvv02/iter_obj_mks_struct;'])

         case 'elmvv03'
              eval(['load ' gatools_root '/tokamaks/iter/make/elmvv03/iter_obj_mks_struct;'])

         case 'hairpin'
              eval(['load ' gatools_root '/tokamaks/iter/make/hairpin/iter_obj_mks_struct;'])

         case 'standard'
              eval(['load ' gatools_root '/tokamaks/iter/make/standard/iter_obj_mks_struct;'])

         case 'vs3a'
              eval(['load ' gatools_root '/tokamaks/iter/make/vs3a/iter_obj_mks_struct;'])

         case 'vvopt1'
              eval(['load ' gatools_root '/tokamaks/iter/make/vvopt1/iter_obj_mks_struct;'])
	      
	 case 'current'
	      eval(['load ' gatools_root '/tokamaks/iter/make/2011/iter_obj_2011;'])

         otherwise
            wait('ERROR load_tok_objects: invalid configuration name for iter') 
	    tds=load_tok_objects('help iter');  
       end %iter  

   case 'help fdf'
	disp('Usage:')
	disp('   tok_data_struct = load_fdf_objects(<config>);')
	disp('   (Put <config> in single quotes...)')
	disp('Available configurations:')
	disp('   fdf11: 2011: 16 Bp probes, 29 flux loops, 22 VF coils, 52 VV cells')
 
   case 'fdf'
     switch config_name
	case 'fdf11'
   	   eval(['load ' gatools_root '/tokamaks/fdf/make/fdf11/fdf_obj_fdf11.mat;'])
	   
	case 'current'
   	   eval(['load ' gatools_root '/tokamaks/fdf/make/fdf11/fdf_obj_fdf11.mat;'])
      
      	otherwise
      	   wait('ERROR load_tok_objects: invalid configuration name for fdf')
	   tds=load_tok_objects('help fdf');
     end %fdf
      
   otherwise
	   wait('ERROR load_d3d_objects: invalid configuration name')
   end %tok
