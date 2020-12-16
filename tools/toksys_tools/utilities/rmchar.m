function result = rmchar(string,char)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:   result = rmchar(string,char)
%
%  PURPOSE:  Remove all instances of a specified character from a string.
%
%  INPUT:
%	string = string to remove character from 
%	char   = character to remove
%
%  OUTPUT:
%	result = string after removing char
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	2/10/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = findstr(char,string);
if(~isempty(idx))
   result = string(1:idx(1)-1);
   for k=2:length(idx)
      result = [result string(idx(k-1)+1:idx(k)-1)];
   end
end
result = [result string(idx(end)+1:end)];
