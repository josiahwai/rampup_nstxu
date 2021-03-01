function make_psc_file(figures,filename)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  make_psc_file(figures,filename)
%
%  PURPOSE:  Make a color postscript file for printing or viewing with "gs".
%
%  INPUT:
%	figures = list (vector) of figure numbers
%	filename = name of postscript file to create
%
%  OUTPUT: postscript file
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	3/22/00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(['!rm '  filename]);
for k=figures
   figure(k)
   eval(['print -dpsc2 -append ' filename]);
end
