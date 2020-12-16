function plot_tok_geo_config(tok)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  SYNTAX:  plot_tok_geo_config(tok)
%
%  PURPOSE:  call plot_tok_geo for multiple vacuum data objects
%
%  INPUT: tok: name of the device
%       
%
%  OUTPUT: Figure plotted to screen and .ps file
%
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  MJL 20121002
%
%  MODIFICATIONS: 
%	So far only ITER and D3D cases are implemented
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% Inits
  %path(['/u/lanctot/matlab/toksysv0.0/analysis'],path)
  
  %% Plot options
  %	iblackbg = flag: 0= white figure background, [1]=black figure background
%       ipltpsi  = flag: 1= contour equilibrium psizr, >1= #contours to plot
%				(default = 1 if psizr exists, else 0)
%       ipltB    = flag: 1= contour equilibrium |B|, >1= #contours to plot [0]
%       ipltJ    = flag: 1= contour equilibrium jphi, >1= #contours to plot [0]
%     (Only one of ipltpsi, ipltB, ipltJ may be specified; all need equil_data.)
%	ipltlim  = flag: 1= plot limiter switch [1]
%	ipltfl   = flag: 1= plot Flux Loops [0]
%	ipltbp   = flag: 1= plot B-Probe indices [0]
%	ipltgap  = flag: 1= plot gap locations [0]
%	ilabeleq = flag: 1= contour label flux values [0]
%	ilabelfc = flag: 1= label PF's with indices, 
%                       >1 => use name labels, with fontsize=ilabelfc [0]
%	ilabelcc = flag: 1= label PF coil STATES with indices (cannot also set ilabelfc) 
%       Pcc      = matrix that computes PF currents from states (required if ilabelcc=1)
%	ilabelvv = flag: 1= label VV elements with indices [0]
%	ilabelfl = flag: 1= label FL's with indices [0]
%                       >1 => use name labels, with fontsize=ilabelfl [0]
%	ilabelbp = flag: 1= label BP's with indices,
%                       >1 => use name labels, with fontsize=ilabelbp [0]
%	ilabelgap= flag: 1= label gaps with indices [0]
%	idxvv    = indices of vacuum vessel to plot (default = all)
%	idxfl    = indices of flux loops to plot (default all, used only if ipltfl)
%	idxbp    = indices of Bprobes to plot (default all, used only if ipltbp)
%	vvgroup  = grouping vector for vacuum vessel elements, used only if 
%			ilabelvv=1 to define vessel element labels
%       ileftright: bit 1 = plot left side, bit 0 = plot right side (default)

  options = struct('iblackbg',0,...
  		   'ipltpsi',0,... 	   	   
  		   'ipltB',0,...   
  		   'ipltJ',0,...   
  		   'ipltlim',1,...   
  		   'ipltfl',0,...       
  		   'ipltbp',1,...       
  		   'ipltgap',0,...      
  		   'ilabeleq',0,...     
  		   'ilabelfc',1,...     	     
  		   'ilabelcc',0,...     
  		   'Pcc',0,...          
  		   'ilabelvv',0,...     
  		   'ilabelfl',0,...          
  		   'ilabelbp',0,...          
  		   'ilabelgap',0); %,...    
  		  % 'idxvv',,...        
  		  % 'idxfl',,...        
  		  % 'idxbp',,...        
  		  % 'vvgroup',);

  %% Get tok and config
  switch tok
  
   case 'd3d'
        d3d_startup
        configs = {'with_ADP',...
        	  'after_ADP'};
	configs=configs(end:-1:1);
	pconfigs=[1,2];	
	efit_grid='6565';
		  
  	%% Plots
  	sconfigs = size(configs);
	spconfigs = size(pconfigs);
  	for k = 1:spconfigs(2)
   	  tok_data_struct = load_d3d_objects(configs{pconfigs(k)},efit_grid);  %curly braces returns contents of cell array
   	  plot_tok_geo(tok_data_struct,options)
	  
	  stitle=title([tok ' Geometry -' configs{pconfigs(k)}]);
	  %set(stitle,'Position',[1,1])
   	  print('-dpsc',configs{k})
	  
	%k=waitforbuttonpress
	%hold on
  	end

   case 'nstx'
        nstx_startup
        configs = {'02feb',...
        	   '05apr'};
	pconfigs=[1,2];	
	efit_grid='6565';
		  
  	%% Plots
  	sconfigs = size(configs);
	spconfigs = size(pconfigs);
  	for k = 1:spconfigs(2)
   	  tok_data_struct = load_nstx_objects(configs{pconfigs(k)},efit_grid);  %curly braces returns contents of cell array
   	  plot_tok_geo(tok_data_struct,options)
	  
	  stitle=title([tok ' Geometry -' configs{pconfigs(k)}]);
	  %set(stitle,'Position',[1,1])
   	  print('-dpsc',configs{k})
	  
	%k=waitforbuttonpress
	%hold on
  	end   
   
   case 'east'
        east_startup
        configs = {'2008',...		
        	   '2010',...
		   '2012'};		
	pconfigs=[1,2,3];			
	efit_grid='6565';
		  
  	%% Plots
  	sconfigs = size(configs);
	spconfigs = size(pconfigs);
  	for k = 1:spconfigs(2)
   	  tok_data_struct = load_east_objects(configs{pconfigs(k)},efit_grid);  %curly braces returns contents of cell array
   	  plot_tok_geo(tok_data_struct,options)
	  
	  stitle=title([tok ' Geometry -' configs{pconfigs(k)}]);
	  %set(stitle,'Position',[1,1])
   	  print('-dpsc',configs{k})
	  
	%k=waitforbuttonpress
	%hold on
  	end      
	
   case 'kstar'
       kstar_startup
        configs = {'2008',...		
        	   '2009',...		
		   '2010',...
		   '2010_Sabbagh',...
		   '3333',...
		   '6565'};		
	pconfigs=[1];		      
	efit_grid='6565';		
		  			
  	%% Plots
  	sconfigs = size(configs);
	spconfigs = size(pconfigs);
  	for k = 1:spconfigs(2)
   	  tok_data_struct = load_kstar_objects(configs{pconfigs(k)},efit_grid);  %curly braces returns contents of cell array
   	  plot_tok_geo(tok_data_struct,options)
	  axis([0 5 -4.5 4.5]);
	  stitle=title([tok ' Geometry -' configs{pconfigs(k)}]);
	  %set(stitle,'Position',[1,1])
   	  print('-dpsc',configs{k})
	  
	%k=waitforbuttonpress
	%hold on
  	end      
   
   case 'iter'
       iter_startup
       configs = {'vs3a',...
       		  '2011',...
  		  '2010v3p3',...
  		  '2010',...
  		  'elmvac03',...
  		  'elmvv03',...
  		  'elmvv02',...
		  'cc_ins_vv',...
		  'vvopt1',...
  		  'hairpin',...
  		  'standard'};
	configs=configs(end:-1:1);
	pconfigs=[10,11];
		  
  	%% Plots
  	sconfigs = size(configs);
	spconfigs = size(pconfigs);
  	for k = 1:spconfigs(2)
   	  tok_data_struct = load_iter_objects(configs{pconfigs(k)});  %curly braces returns contents of cell array
   	  plot_tok_geo(tok_data_struct,options)
	  
	  stitle=title([tok ' Geometry -' configs{pconfigs(k)}]);
	  %set(stitle,'Position',[1,1])
   	  print('-dpsc',configs{k})
	  
	%  k=waitforbuttonpress
	%hold on
  	end

   case 'fdf'
       fdf_startup
        configs = {'fdf11'};		
	pconfigs=[1];		      
	efit_grid='6565';		
		  			
  	%% Plots
  	sconfigs = size(configs);
	spconfigs = size(pconfigs);
  	for k = 1:spconfigs(2)
   	  tok_data_struct = load_fdf_objects(configs{pconfigs(k)},efit_grid);  %curly braces returns contents of cell array
   	  plot_tok_geo(tok_data_struct,options)
	  axis([0 5 -4.5 4.5]);
	  stitle=title([tok ' Geometry -' configs{pconfigs(k)}]);
	  %set(stitle,'Position',[1,1])
   	  print('-dpsc',configs{k})
	  
	%k=waitforbuttonpress
	%hold on
  	end         
  end




