function truncate_header(infile,outfile,nlines)
%   @(#)truncate_header.m       1.1       14/02/14
%  USAGE: truncate_header
%
%  PURPOSE: remove header from a file
%
%  INPUTS:
%
%  OUTPUTS: 
%
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
fid=fopen(infile,'rt');
fid2=fopen(outfile,'wt');
id=0;
a=fgets(fid);
while(ischar(a))

  id=id+1;
  
  if id <= nlines
   a=fgets(fid);
   continue
  else
   fprintf(fid2,a);
  end
  
  a=fgets(fid);
  
end
