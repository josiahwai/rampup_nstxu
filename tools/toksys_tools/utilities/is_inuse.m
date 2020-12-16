function [inuse] = is_inuse(name)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  
%
%  PURPOSE:  Determine if a name is already in use as a function name
%
%
%  INPUT:
%	name = string defining name to be tested
%
%  OUTPUT:
%	inuse = returns:
%		1 if already in use as a function
%		0 if not in use
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	1/29/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = which(name);

if strcmp(s,'')
   inuse = 0;
else
   inuse = 1;
end
