function string = intvec2str(intvec)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX: string = intvec2str(intvec)
%
%  PURPOSE:  Convert vector of integers to string.  Puts commas between 
%	entries.
%
%  INPUT:
%	intvec = vector of integers
%
%  OUTPUT:
%	string = vector converted into string
%
%  WRITTEN BY:  Mike Walker 	ON 	??/94
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


string = [];
len = length(intvec);

for k=1:len-1
   string = [string int2str(intvec(k)) ', '];
end;     
string = [string int2str(intvec(len))];
