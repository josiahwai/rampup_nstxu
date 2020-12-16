function whatfunction(name_string)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  whatfunction(name_string)
%
%  PURPOSE:  Find a function name given only a partial name string.
%		(See also lookfor, search_path)
%
%  INPUT:
%	name_string = function name string to match
%
%  OUTPUT:  names of functions containing name_string
%
%  RESTRICTIONS:  Bug in MATLAB restricts when this command will work.  
%	Uses the "unix" command to find m-files in a directory and sometimes
%	this command doesn't work.
 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	9/28/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = path;
k2 = findstr(P,':');
k1=[1 k2+1];
k2=[k2-1 length(P)];

ndirectories = length(k2);
for k=1:ndirectories
   directory =  P(k1(k):k2(k))
   [ier,files] = unix(['ls -1 ' directory '/*.m']);

   if ier~=0
      disp([' Unix command failed for directory ' directory])
      return;
   else
      file_ends = findstr('.m',files)+2;
      file_starts = [1 file_ends(1:length(file_ends)-1)+1];
      for j=1:length(file_ends)
         mfile = files(file_starts(j):file_ends(j));
         kk = findstr(mfile,'/');
	 mfile = mfile(kk(length(kk))+1:length(mfile));
         kk = findstr(mfile,'.m');
         mfile = mfile(1:kk-1);
         if length(name_string) < length(mfile) & findstr(mfile,name_string)
	    disp(mfile)
         end
      end
   end
end
