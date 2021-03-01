% link_fine1.m: Information and Matlab compile of fine1.c wrapper S-function to fine1.f
%
% Note: fine points to this file for running
%
% =============================================
% file structure:
% fine.m  -      help and wrapper for fine1.mex
% fine1.f -      vector Interface code for fine.for
% fine.f  -      fine code
% fineg.f -      Mex gateway routine
% fine1.script - how to link files in matlab
%
% =============================================
% To link use:
%  mex fine1.f fine.f fine1g.f
% or: 
% fmex fine1.f fine.f fine1g.f
% =============================================

 
%If the libf2c.a library is not on the library path, 
%you need to add the path to the mex process explicitly 
%with the -L command. 
%For example:
%mex -L/usr/local/lib/ -lf2c sfun_atmos.c sfun_atmos_sub.o

mex -L/m/linux/matlabR2014b/bin/glnxa64/ -lmex -lmx -lmat fine1.F fine.f
