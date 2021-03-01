 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  read_namelist
%
%  PURPOSE:  Read a FORTRAN namelist file and create variables in the
%	matlab workspace defined by variables in the namelist.
%
%  INPUT: [default]
%    filename = name of namelist file to read
%    to_upper = 1; % makes all variables UPPER case [0]
%    to_upper =-1; % makes all variables lower case [0]
%
%  OUTPUT:
%    varlist = list of variables read from namelist file
%
%    All variables found in namelist(s) within the file are made into
%    matlab variables. Duplicate names are overwritten with subsequent data.
%
%  RESTRICTIONS:
%	(1) won't work if '=' signs in quotes in any of data
%	(2) won't work for multiple (*) string variables 
 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	10/16/02
%  JAL3/30/2004 to include to_upper functionality of original routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO:
% (1) put in code to ignore the special characters $, &, /, and = if they
%	follow an odd number of quote (') symbols (i.e. they are embedded
%	in text strings).
% (2) make multiple entries (*) work with string variables.

debug=0;
if exist('to_upper')~=1 to_upper=0; end %default is accept case in namelist

varlist = [];
fid = fopen(filename,'r');

icnt = 0;
linecnt = 0;
separators = [' ';','];	% space, comma; tabs will be replaced with spaces

line1=0;
nlcnt=0;
nlnames = '';
while line1~=-1
   line1 = fgets(fid);

% find start of namelist

   totstr = '';
   if(~isempty(findstr(line1,'$')) | ~isempty(findstr(line1,'&')))
      fprintf('Begin processing %s\n',line1)
      nlcnt = nlcnt+1;
      temp = deblank(line1);
      if(size(nlnames,2)~=0)
         if(length(temp)>size(nlnames,2))
            bb = blanks(length(temp)-size(nlnames,2));
            nlnames = [nlnames bb*ones(size(nlnames,1))];
         end
      end
      nlnames(nlcnt,1:length(temp)) = temp;
      while line1~=-1
         line1 = fgets(fid);

% find end of namelist

         idx = union(findstr(line1,'$'),findstr(line1,'&'));
         if(~isempty(idx))
           if(idx(1)>2)
               totstr = [totstr line1(2:idx(1)-1) ' '];	% save string in between
           end
           break
         end;
         if(length(line1)>1 & strcmp(line1(2),'/'))
           break
         end;

         totstr = [totstr line1(2:end) ' '];		% save string in between
      end
   end

% If a namelist is found and is nonempty, then parse the string in between.

   if(~isempty(totstr))	
      totstr = packnlstr(totstr);	%get rid of extraneous spaces, tabs, ...

      eq_idx = findstr(totstr,'=');	% find all equal signs
      ll = length(eq_idx);

      for k=1:ll			% loop over all equal signs found

         var = deblank(totstr(1:eq_idx(1)-1))	% variable is string before =

         idx1 = findstr(var,'(');	% if variable contains parens, then...
         if(~isempty(idx1))
            idx2 = findstr(var,')');
            idx3 = findstr(var,',');
            if(isempty(idx3))
               istart = sscanf(var(idx1+1:idx2-1),'%d',1);	% index to insert values
            else
               var(idx3) = ' ';
               istart = sscanf(var(idx1+1:idx2-1),'%d',Inf);	% n-dim	
            end
            var = var(1:idx1-1);
         else
            istart = 1;			% else insert values starting at index 1
         end

         if(length(eq_idx) > 1)		% if more than one = sign present, ...

            mlen = eq_idx(2);
            partstr = totstr(1:mlen);
            sep_idx = union(findstr(partstr,' '),findstr(partstr,','));
% remove indices for separators "protected" by parens
            idx1 = findstr(partstr,'(');
            idx2 = findstr(partstr,')');
            if(length(idx2)>length(idx1))
               fprintf('mismatch in parens %s\n',partstr);
            elseif(~isempty(idx1))
               idx3 = [];
               for i=1:length(idx1)
                  idx3=union(idx3,find(idx1(i)<sep_idx & sep_idx<idx2(i)));
               end
               sep_idx = setdiff(sep_idx,sep_idx(idx3));
            end

            dd = diff(sep_idx);
            for i = length(sep_idx)-1:-1:1
               if(dd(i)>1), break;, end;
            end
   
            if(isempty(dd))
               max_idx = sep_idx(1);
            else
               max_idx = sep_idx(i+1);	% assumes no space between var and "="
            end
   
            datastr = totstr(eq_idx(1)+1:max_idx);
            totstr = totstr(max_idx+1:end);
            eq_idx = findstr(totstr,'=');

         else				% if last equal sign, remainder is data

            mlen = length(totstr);
            datastr = totstr(eq_idx(1)+1:end);
         end

         if debug			% debugging only
           var, datastr
           wait
         end

% var = variable name, datastr = string containing data, are available here

         idx = findstr(datastr,'''');		% check for strings in data

         if(~isempty(idx))			% if strings in data
            idata = 0;
            data = '';
            len = length(idx);
            if(modulo(len,2))	% if unfinished quote
               fprintf('ERROR - non-matching quote in %s\n',var)
            end
            for k=1:floor(len/2);
               data0=sscanf(datastr(idx(2*k-1)+1:idx(2*k)-1),'%s');
               idata= idata + 1; % count data added to array
               data(idata,1:length(data0)) = data0;
            end
   
         else					% else if no strings in data
   
            data = [];
            [token,remainder] = strtok(datastr,separators);
            while(~isempty(remainder))
               idx = findstr(token,'*');	% see if multiple values
               if(~isempty(idx))
                  nvals = str2num(token(1:idx(1)-1));
                  value = str2num(token(idx(1)+1:end));
                  data = [data; value*ones(nvals,1)];
               else
                  data = [data; str2num(token)];
               end
               [token,remainder] = strtok(remainder,separators);
            end
         end
         li = length(istart);
         if(isempty(strmatch(var,varlist)))
            varlist = strvcat(varlist,var);
         end
         if(li>1)
            str= [var '(' int2str(istart(1))];
            for i=2:li
               str= [str ',' int2str(istart(i))];
            end
         else
            iend = istart + size(data,1)-1;
            str= [var '(' int2str(istart) ':' int2str(iend) ',:'];
         end
         if to_upper==1
            str=upper(str);
	 elseif to_upper==-1
	    str=lower(str);
   	 end
         str = [str ') = data' ';'];
	 eval(str);
      end
      if(~isempty(findstr(line1,'$')))
         fprintf('End processing %s\n',line1)
      end
   end
end

fclose(fid);
clear separators sep_idx token remainder totstr value var varc nvals
clear mlen ll linecnt line1 k idx icnt i temp

