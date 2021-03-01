function write_ascii_table(table,filename,format)
 %
%  SYNTAX:  write_ascii_table(table,filename,format)
%
%  PURPOSE: Write out table to ascii file for easy import into MSWord
%
%  INPUT:
%	table = table to write
%	filename = name of file to write table into
%	format = can be either string or integer.
%		string = defines format for write to file
%		integer = defines number of digits to write out for each entry
%
%  OUTPUT:
%	file with formatted data (name=filename)
 
%  RESTRICTIONS:
%
%  METHOD:  The matlab function dlmwrite almost does what you need, but doesn't
%	allow you to specify data format.
%
%  WRITTEN BY:  Mike Walker 	ON 	6/6/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%shift = 10^format;
%temp = round(shift * table)/shift;

if(isstr(format))	% if a string, must be a format specification
   fid = fopen(filename,'w');
   for k=1:size(table,1)
      fprintf(fid,format,table(k,:))
   end
   fclose(fid);
else
   temp = table;
   dlmwrite(filename,temp,'\t',0,0,format);
end
