function strout = insert_string(str2insert,string,indices)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  strout = insert_string(str2insert,string,indices)
%
%  PURPOSE:  Insert a specified substring within a given string at one
%	or more locations.
%
%  INPUT:
%	str2insert = string to insert
%	string     = string to insert it into
%	indices    = indices in string at which to insert (length>=1) 
%			(inserts BEFORE the character located at each index)
%
%  OUTPUT:
%	strout     = string resulting from the insertions
%
%  RESTRICTIONS:
 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	11/7/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin < 3)
   disp('insert_string ERROR: must be 3 input arguments')
   return;
end

ind = indices;

nstrings_to_insert = length(ind);

% handle cases where indices are not in order

idx = find(ind==1);
if(~isempty(idx))
   strout = str2insert;
   nstrings_to_insert = nstrings_to_insert-1;
   ind = setdiff(ind,1);
else
   strout = '';
end

i1=1;
while(nstrings_to_insert)
   i2 = min(ind);
   strout = [strout string(i1:i2-1) str2insert];
   nstrings_to_insert = nstrings_to_insert-1;
   ind = setdiff(ind,i2);
   i1 = i2;
end

strout = [strout string(i1:end)];


