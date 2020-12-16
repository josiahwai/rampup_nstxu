 %
%  SYNTAX:  make_model_script
%
%  PURPOSE:  Script to build output objects needed to build isoflux models.
%
%  INPUT:
%	shotnum	= shot number 
%	time    = equilibrium time (ms)
%       tokamak = name of device to build model for
%       tok_data_struct = vacuum data objects for this device
%	out_irzresp = output irzresp value (1, 2, or 3)
%       eqdir   = directory where efit equilibrium files are located
%       efit_source = tree in mdsplus to get equilibrium, e.g. EFIT01, EFITRT
%               (only one of either eqdir or efit_source should be specified)
%
%  OUTPUT: 
%	out_objs_<shotnum>_<time>.mat  - mat file containing isoflux objects
%	out2_objs_<shotnum>_<time>.mat - mat file containing magnetics objects
% 
%  RESTRICTIONS: Requires existence of script mdl_parms_<shotnum>_<time>.m
%	in your matlab path.
 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	5/11/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)make_model_script.m	1.3 09/01/09

%	dyn_irzresp = dynamic irzresp value (1 or 2)
%	model_dir   = directory to store output *.mat files containing
%		output objects (default = './')

if(~exist('eqdir') & ~exist('efit_source'))
  wait('ERROR make_model_script: either eqdir or efit_source must be specified')
end

outflag=3;		% generate both output and output2 data objects
output_irzresp = out_irzresp;	% rigid r and/or z in output equations
eval(['mdl_parms_' int2str(shotnum) '_' int2str(time)])
build_model
