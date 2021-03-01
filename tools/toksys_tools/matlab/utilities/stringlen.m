function len = stringlen(string)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  len = stringlen(string)
%
%  PURPOSE: Find length of string as defined by locating last nonblank character
%
%  INPUT:
%	string = text string to find length of
%
%  OUTPUT:
%	len = length of character string
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	??
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = length(string);   	% if no blank spaces found
for k=length(string):-1:1
   if ~strcmp(' ',string(k))
      len = k;
      return;
   end;
end;

      
