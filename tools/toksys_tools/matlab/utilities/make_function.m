function make_function(infile,outfile)
 %
%  USAGE: make_function(infile,outfile)
%
%  PURPOSE: make a matlab script a function that returns a structure
%
%  INPUTS:
%	The script is read from infile.
%
%  OUTPUTS: 
%	The function is written into outfile.
 
%  RESTRICTIONS: 
%
%  METHOD:  
%
%  VERSION 
%
%  WRITTEN BY:  Matthew J. Lanctot on March 19 2013
%
%  MODIFICATION HISTORY:
% 	2013-03-19	Created
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen([infile '.m'],'rt');
fid2=fopen([outfile '.m'],'wt');
fprintf(fid2,['function ws = ' outfile '\n']);
id=0;
a=fgetl(fid);
while(ischar(a))
  id=id+1;
  fprintf(fid2,'%s ; \n',a); 
  a=fgetl(fid);
end
fprintf(fid2,['ws=ws2struct();\n']);
fclose(fid2);
fclose(fid);
