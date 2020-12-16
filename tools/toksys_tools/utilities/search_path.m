function results = search_path(string,path_depth,search_depth)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  results = search_path(string,path_depth,search_depth)
%
%  PURPOSE:  String search of the matlabpath for documentation about a function,
%	or model, or whatever, whose documentation contains the search string.
%	(See also lookfor, whatfunction)
%
%  INPUT:
%	string       = string to search for
%	path_depth   = one of the following (optional argument):
%			- integer specifying how many directories in
%			path to search
%			- string specifying last directory to search 
%				(e.g. 'GAcontrol')
%			- 'all', i.e. search all directories (default)
%	search_depth = 'shallow' or 'deep' (optional, default='shallow')
%		   'shallow' = look at Contents.m files only
%		   'deep'    = search both Contents.m files and all m-files
%				(deep search takes awhile)
%
%  OUTPUT:
%	results = string containing documentation lines where strings matched
%		and files where found
% 
%  RESTRICTIONS:  Bug in MATLAB restricts when the 'deep' option will work.  
%	Uses the "unix" command to find m-files in a directory and sometimes
%	this command doesn't work.
 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	9/27/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
   search_depth='shallow';
end
if nargin < 2
   path_depth = 'all';
elseif isempty(path_depth)
   path_depth = 'all';
end

P = path;
k2 = findstr(P,':');
k1=[1 k2+1];
k2=[k2-1 length(P)];

if ischar(path_depth)
   if strcmp(path_depth,'all')
      ndirectories = length(k2);
   else
      k=findstr(path_depth,P);
      [dummy,ndirectories]=max(find(k>k1));
   end
else
   ndirectories = path_depth;
end

results = char(32*ones(1,80));
line_cnt = 0;
for k=1:ndirectories
   directory =  P(k1(k):k2(k));
   temp = [directory ':'];
   line_cnt = line_cnt+1;
   results(line_cnt,1:length(temp)) = temp; 
   fid = fopen([directory '/Contents.m']);

   if fid>0
	 line = 0;
	 while line~=-1
	    line = fgets(fid);
	    if(any(findstr(line,string)))
	       line_cnt = line_cnt+1;
	       results(line_cnt,1:length(line)) = line;
	    end
	 end
	 fclose(fid);
   else
      line_cnt = line_cnt+1;
      temp = '   No Contents.m file.';
      results(line_cnt,1:length(temp)) = temp;
   end

   disp(['starting deep search of ' directory ' ...'])
   if strcmp(search_depth,'deep')
      [ier,files] = unix(['ls -1 ' directory '/*.m']);
      if ier~=0
	 disp(' Unable to do deep search because unix command failed')
	 return;
      else
         file_ends = findstr('.m',files)+2;
         file_starts = [1 file_ends(1:length(file_ends)-1)+1];
	 for j=1:length(file_ends)
	    mfile = files(file_starts(j):file_ends(j));
	    if isempty(findstr(mfile,'Contents.m'))
	       fid = fopen(mfile,'r');
               if fid>0
	          line_cnt = line_cnt+1;
	          results(line_cnt,1:length(mfile)) = mfile;
	          line = 0;
      	          while line~=-1
	             line = fgets(fid);
	             if(any(findstr(line,string)))
	                line_cnt = line_cnt+1;
	                results(line_cnt,1:length(line)) = line;
	             end
	          end
	          fclose(fid);
	       else
                  line_cnt = line_cnt+1;
                  temp = ['   Unable to open ' mfile];
                  results(line_cnt,1:length(temp)) = temp;
	       end
	    end
	 end
      end
   end
end
