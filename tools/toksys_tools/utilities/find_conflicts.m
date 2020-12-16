 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  find_conflicts
%
%  PURPOSE:  Find potential conflicts between variable names in your
%	environment and function names in your search path.
%
%  INPUT: none
%
%  OUTPUT: printed warnings of potential conflicts
%
%  RESTRICTIONS:  
%
%  METHOD:  Entire environment is saved at beginning of execution and restore
%	when complete; this may take some time if large workspace is used.
 
%  WRITTEN BY:  Mike Walker 	ON 	2/5/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save xxtempxx
clear tempcell tempstr
tempcell = who;
for k=1:size(tempcell)
   tempstr = tempcell{k};
   clear(tempcell{k})
   if ~strcmp(which(tempstr),'')
      fprintf('potential conflict between variable %s and function\n %s\n\n',...
					tempstr, which(tempstr))
   end
end
clear tempcell tempstr
load xxtempxx

