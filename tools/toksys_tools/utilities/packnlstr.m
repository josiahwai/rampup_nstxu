function str_out = packnlstr(str_in)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  
%
%  PURPOSE:  Pack namelist string (remove all the extraneous separators).
%
%  INPUT:
%	str_in = string to pack
%
%  OUTPUT:
%	str_out = string after packing
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	10/17/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove any tabs in string
   idx = findstr(str_in,'	');
   str_in(idx) = ' ';

str_out = str_in;

ichange = 1;
while(ichange)
   ichange = 0;
   [str_out,ichange] = xpackx(str_out,' (','(',ichange);
   [str_out,ichange] = xpackx(str_out,'( ','(',ichange);
   [str_out,ichange] = xpackx(str_out,' =','=',ichange);
   [str_out,ichange] = xpackx(str_out,') ',')',ichange);
end

function [packed_str,ichange] = ...
	 xpackx(str_to_pack,str_to_replace,replacement,ichange)
% utility function

   tempstr = str_to_pack;
   packed_str = '';
   idx = findstr(tempstr,str_to_replace);
   if(~isempty(idx)) ichange = 1; end;
   while(~isempty(idx))
      packed_str = [packed_str tempstr(1:idx(1)-1) replacement];
      tempstr = tempstr(idx(1)+2:end);
      idx = findstr(tempstr,str_to_replace); 
   end
   packed_str = [packed_str tempstr];
