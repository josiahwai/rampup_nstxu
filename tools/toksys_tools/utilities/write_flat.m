 function dum = write_flat(dataobj,output_file,flag)
%
% function dum = write_flat(dataobj,output_file,flag)
%
%   Matlab m-file to write dataobj to flatfile output_file.flat in
%   flat format = single ascii number per line followed by CR.
%
% Usage:
%   >>write_flat(data,'filename');
%
% Inputs:
%   dataobj = data object desired to write to file output_file in flat format
%   output_file = (quoted) string denoting name of output file
%   flag (optional) = 1 to write matrix dimensions in first 2 lines, 0 to omit
%		(in this case, 1st line = # of rows, 2nd line=# of cols)
% Restrictions:
%   None

   fid = fopen([output_file,'.flat'],'w');
   dum = fid;
   if nargin>2, fprintf(fid,'%4i\n',size(dataobj)); end
   fprintf(fid,'%17.10e\n',dataobj);  %used to be 13.6
   fclose(fid);
   
