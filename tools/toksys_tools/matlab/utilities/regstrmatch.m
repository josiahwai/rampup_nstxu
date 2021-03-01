function indices = regstrmatch(string_pattern,stringarr)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  
%
%  PURPOSE:  Regular string match a'la UNIX string matching.
%
%  INPUT:
%
%
%  OUTPUT:
%
%  RESTRICTIONS: Only one of * or % is permitted in string_pattern.
 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	9/24/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

istar = findstr(string_pattern,'*');
if isempty('istar')
   nstar = 0;
else
   nstar = length(istar);
end
iper = findstr(string_pattern,'%');
if isempty('iper')
   nper=0;
else
   nper = length(iper);
end
if nper & nstar
   disp('only one of * or % is permitted')
end

strarr = stringarr;

if 0
strarr = [stringarr char(32*ones(size(stringarr,1),1))];	% 32 is blank
for k=1:size(strarr,1)
   ii = strlen(strarr(k,:));
   strarr(k,ii+1:ii+1) = '%';
end
end

lpattern = length(string_pattern);

strwcmp

if nstar					% an asterisk was used

   if istar(1)~=1 
      index = strmatch(string_pattern(1:istar(1)-1),strarr);
      strarr = strarr(index,:);
   end
   icnt=0;
   for k=1:size(strarr,1)
      if strwcmp(string_pattern,strarr(k,:))
	 icnt = icnt+1;
	 indices(icnt) = index(k);
      end
   end

elseif nper					% a percent sign was used

else						% else just regular strmatch

   indices = strmatch(string_pattern,stringarr);

end
