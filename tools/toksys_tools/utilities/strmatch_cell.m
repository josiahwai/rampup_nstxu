function idx = strmatch_cell(str,cellarray,exact)
 %
%  SYNTAX:  idx = strmatch_cell(str,cellarray)
%
%  PURPOSE:  Replacement for  matlab strmatch function, which doesn't work 
%	correctly for cell arrays.
%
%  INPUT:
%	str = string
%	cellarray = cell array (n x 1) of strings
%
%  OUTPUT:
%	idx = indices of matching entries in cellarray

%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	6/18/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)strmatch_cell.m	1.1 06/21/10

if(nargin<2)
   wait('ERROR strmatch: strmatch_cell requires at least 2 arguments')
   return;
end

if(nargin==3)
   if(~strcmp(exact, 'exact'))
      wait('ERROR strmatch: 3rd argument can only be "exact"')
      return;
   end 
   ii = strcmp(str,cellarray);
   idx = find(ii);
else
   idx = strmatch(str,cellarray);
end


