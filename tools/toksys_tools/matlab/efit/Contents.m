 %
% Matlab EFIT Area (Plus MDS-EFIT Read)
%
% read_afile      Reads a0-eqdsk file and puts variables in Matlab env
% read_gfile      Reads g0-eqdsk file and puts variables in Matlab env
% read_gfile_func Same as read_gfile but function call & structure output
% read_gfile_tok  read gfile structure for particular machines d3d,kstar...
% read_mfile_func read mfile in netcdf format into struct (ONLY for Matlab 7 V14)
% plot_efit       Plots efit flux contours
% get_efit_data1  Gets data for lots of EFIT slices and plots it
%                 Note: Conflict with daves get_efit_data in /plasma_models 7/03
% efit_movie      Make movie from EFIT flux data
% run_efit        Runs EFIT from inside Matlab much as EFIT runs from unix
% get_efit_profile Get EFIT current profile parametrization from esave.dat
% read_esave      Read the esave.dat file (fortran mex function)
% efit_lp         Calculate plasma inductance from gfile information
% plasma_field    Calculate plasma flux and field at arbitrary rpt,zpt points
% shot_from_gfile returns shot,time,dir (or tree,server) based on gfile name
%
% New Routines to extract EFIT "eq." structure from MDS: (work in progress)
%
% cc_efit_to_tok  converts data fetched from mdsplus or gfile to toksys units
% eq_mds          Fetch an EFIT tree out of mdsplus
% eq_mk_time_last Modifies eq. to put time index as last index in all vectors
% gdata_to_eq     Reads gdata from read_gfile_tok and puts into "eq." structure
% read_mds_eqdsk  Return data like output of read_gfile_tok.m from mdsplus
% read_mds_eq_func Return data like output of read_gfile_func.m from mdsplus
% equil_to_gdata  Identical output as read_gfile_func but reads from mds+
%
% Below eq. routines work in progress:
% eq_time_lim    Reduces time data in sructure eq. within limits; also makes var
% eq_plasma_cur  Generates the plasma current from flux inside eq.
% plasma_current Generates plasma current and current density from efit psi
% inside_plasma  Function Determines all points inside plasma
% eq_ga_env      Makes standard G (& later A file) variables from eq. into env.
%
% Misc EFIT scripts
%
% calc_fpcurr    Emulates fpcurr in efit (ffprime calculation from polynom.)
% calc_ppcurr    Emulates ppcurr in efit (pprime calculation from polynom.)
% calc_bsffel    Emulates bsffel in efit (iparm'th term of poly representation)
% calc_bsppel    Emulates bsppel in efit (iparm'th term of poly rep. of pprime)
% gfile_def      Generates structure of EFIT gfile variable definitions 
% test_get_efit_profile - Example for get_efit_profile.m,current reconstruction
%
% Emmanuel Joffrin's EFIT tools:
%
% fluxfun - function to extract profile data from EFIT and make flux mappings
% fluxmoy - called by fluxfun...
% isoflux - called by fluxfun...
