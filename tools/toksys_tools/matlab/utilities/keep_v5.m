 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  
%
%  PURPOSE:  
%
%
%  INPUT:
%	keepstr = string containing names of variables to keep
%
%  OUTPUT:
%
%  RESTRICTIONS:
 
%  METHOD:  
%
%
%  WRITTEN BY:  Mike Walker 	ON 	5/4/95
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ttxxxxx = who;	% creates a cell array

len = length(ttxxxxx);
kkkxxx = blanks(30); % 30 is max variable name size
kkkxxx(1:7) = 'keepstr';

for k=1:len
   temp = ttxxxxx{k}
   if isempty(findstr(keepstr,temp)) & any(temp~=kkkxxx(1:length(temp)))
      eval(['clear ' temp])
   end
end

clear kkkxxx keepstr ttxxxxx temp len k
