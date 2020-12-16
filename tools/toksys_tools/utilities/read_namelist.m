 function [nlnames,nlstr,nlall]= read_namelist(filename,to_upper);
 %
%  SYNTAX:
%       FUNCTION CALL:
%          [nlnames,nlstr,nlall]= read_namelist(filename,to_upper);
%
%       SCRIPT CALL (compatable with old script)
%          read_namelist; % SCRIPT FORM reads & writes to caller area
%
%  PURPOSE:  Read a FORTRAN namelist file and create variables in the
%	     matlab workspace defined by variables in the namelist.
%
%  INPUT: [default]
%    filename = name of namelist file to read
%    to_upper = 1; % makes variables UPPER case; -1= LOWER case;s [0]
%
%    NOTE: if any input variable missing it looks 1st in "caller" then in "base" 
%
%  OUTPUT:
%
%    nlnames= Name of all namelists found in file
%    nlstr=   Structure of namelist and variables
%    nlall=   total list of all namelists and variable (for write_namelist)
%
%    NOTE: if output variable "nlstr" is not present, routine puts all
%          namelist variables in "caller" routine. Same for "nlall"
%          (duplicate variable names in "caller" routine are overwritten)
%
%    NOTE: To get all structure namelist variables into environment:
%        [names, isstruc]= struct_to_ws(nlstr); % Do 1st level
%        if any(isstruc)                         % Loop to do 2nd level
%          id= find(isstruc);
%          for ii= 1:length(id)
%            str= ['[names1,isstruc1]= struct_to_ws(' char(names(id(ii))) ');'];
%            eval(str)
%          end  % for
%        end    % if any
%%        clear nlstr
%
%  RESTRICTIONS:
%	won't work if '=' signs in quotes in any of string data
%	won't work for multiple (*) string variables 
 
%  METHOD: Good Template for converting a script to a dual function/script type 
%
%  WRITTEN BY:  Jim Leuer 	ON 	10/16/02
%  Major part taken from Walkers old version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)read_namelist.m	1.3 09/30/13

% TODO:
% (1) put in code to ignore the special characters $, &, /, and = if they
%	follow an odd number of quote (') symbols (i.e. they are embedded
%	in text strings).
% (2) make multiple entries (*) work with string variables.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if nargin<=0 % READ INPUT FROM CALLING or BASE Areas
    ex= 1;
    evalin('caller','filename;','ex=0;');
    if ex
      filename= evalin('caller','filename;');
    else                           % look in base operating space
      ex= 1;
      evalin('base','filename;','ex=0;');
      if ex
        filename= evalin('base','filename;');
      else % Not found in calling or base
        disp('%ERROR read_namelist: No variable ''filename'' in caller or base');
        return
      end
    end
    ex= 1;
    evalin('caller','to_upper;','ex=0;');
    if ex
      to_upper= evalin('caller','to_upper;');
    else                                % look in base operating space
      ex= 1;
      evalin('base','to_upper;','ex=0;');
      if ex
        to_upper= evalin('base','to_upper;');
      else % Not found in calling or base
        to_upper= 0;
      end
    end
   elseif nargin<=1
    ex= 1;
    evalin('caller','to_upper;','ex=0;');
    if ex
      to_upper= evalin('caller','to_upper;');
    else                                % look in base operating space
      ex= 1;
      evalin('base','to_upper;','ex=0;');
      if ex
        to_upper= evalin('base','to_upper;');
      else % Not found in calling or base
        to_upper= 0;
      end
    end   
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open File filename
  [fid, mesg] = fopen(filename,'r');
  if fid<=0
     disp(['%ERROR read_namelist: couldnt open file: ', filename]);
     disp(mesg)
     nlnames= -1;
     nlstrc= [];
     nlall= -1;
     return
   end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  icnt = 0;
  linecnt = 0;
  separators = [' ';','];	% space, comma; tabs will be replaced with spaces

  line1=0;
  nlcnt=0;
  nlnames = '';
  nlall= ''; % all namelist concatanated into one character array
  while line1~=-1
   line1 = fgets(fid);
   if line1==-1 break, end
   
% find start of namelist (note $ or & must be in column 1 or 2 to be namelist
   line0= remove_space(line1,1); % remove leading blank
   if ~isempty(line0)
      line0= line0(1);              % 1st should be $ or & if namelist
   end
   totstr = '';
   if (line0=='$' | line0=='&')
%      if(~isempty(findstr(line0,'$')) | ~isempty(findstr(line1(1),'&')))
      fprintf('Begin processing: %s\n',line1)
      nlcnt = nlcnt+1;
      nlname= remove_space(line1,0);
      nlsym=  nlname(1);     % $ or & symbol used to start name list name
      nlname= nlname(2:end); % name
      nlnames= strvcat(nlnames,nlname);
      nlall=   strvcat(nlall,deblank(line1)); % load in namelist start line
% find end of namelist
      while line1~=-1
         line1 = fgets(fid);
         idx = union(findstr(line1,'$'),findstr(line1,'&')); % problem for str with &
         if(~isempty(idx))
           if(idx(1)>2)
               totstr = [totstr line1(2:idx(1)-1) ' '];	% save string in between
           end
           break
         end;
         if(length(line1)>1 & strcmp(line1(2),'/')) % break on back slash
           break
         end;

         totstr = [totstr line1(2:end) ' '];		% save string in between
      end
   end

   nlterm= deblank(line1); % termination line of namelist

% If a namelist is found and is nonempty, then parse the string in between.

   if(~isempty(totstr))	
%      str= [nlname '_list= '''';']; % initialize namlist_list storage variable nul
%      eval(str); % store list of all variables in namelist_list

      totstr = packnlstr(totstr);	%get rid of extraneous spaces, tabs, ...

      eq_idx = findstr(totstr,'=');	% find all equal signs
      ll = length(eq_idx);

      for k=1:ll			% loop over all equal signs found

         var = remove_space(totstr(1:eq_idx(1)-1),0);	% variable is string before =

         idx1 = findstr(var,'(');	% if variable contains parens, then...
         if(~isempty(idx1))
            idx2 = findstr(var,')');
            istart = sscanf(var(idx1+1:idx2-1),'%d',1);	% index to insert values
            var = var(1:idx1-1);
         else
            istart = 1;			% else insert values starting at index 1
         end

         if(length(eq_idx) > 1)		% if more than one = sign present, ...

            mlen = eq_idx(2);
            sep_idx = union(findstr(totstr(1:mlen),' '), ...
   		findstr(totstr(1:mlen),','));
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

% var = variable name, datastr = string containing data, are available here
% If you search datastr based on '' then entire namelist is not parsed correctly
% Due to presence of newlines, ect 
% Try searching datastr based on comma separator instead. Result: Yes this works
% However in mhdin.dat file, some string variables still not read correctly like MPNAM2
% In this case, last two names are missing when read

         idx = findstr(datastr,'''');		% check for strings in data
         idatatype= 0;                          % 1= data, 2= string, 3=T/F
         if(~isempty(idx))			% if strings in data
	    idatatype= 2;                       % 1= data, 2= string, 3=T/F
            idata = 0;
            data = '';
            len = length(idx);
            if(modulo(len,2))	% if unfinished quote
               fprintf('ERROR - non-matching quote in %s\n',var)
            end
	    idx = findstr(datastr,','); %data is before each comma
	    iis = 1;	    
            for k=1:length(idx)+1;  %+1 is for last entry
	       if k <= length(idx)
	         tmp = datastr(iis:idx(k));
	       else
	         tmp = datastr(iis:end);
	       end
	       idxx = findstr(tmp,'''');
	       if length(idxx)==2    	 %there are two quotes      
                 data0= tmp(idxx(1)+1:idxx(2)-1);
	       elseif length(idxx)==1                     %there is only quote  
	         if idxx < length(tmp)/2  %the quote is before the string
		   data0 = strtrim(tmp(idxx:end));
		   data0 = strrep(data0,' ','');
		   data0 = strrep(data0,'\n','');
		   data0 = strrep(data0,'''','');		   
		 else                        %the quote is after the string
		   data0 = strtrim(tmp(1:idxx));
		   data0 = strrep(data0,' ','');
		   data0 = strrep(data0,'\n','');
		   data0 = strrep(data0,'''','');
		 end
	      
	       end
               idata= idata + 1; % count data added to array
	       if k <= length(idx), iis = idx(k)+1;end
               data= strvcat(data,data0);
	       idxx = [];
            end
         else					% else if no strings in data
   
%          check for logical
           idx_log = union(findstr(upper(datastr),'T'), ...
                     findstr(upper(datastr),'F'));
           if ~isempty(idx_log)                  % Logical True False
	    idatatype= 3;                       % 1= data, 2= string, 3=T/F
            data = logical([]);
            [token,remainder] = strtok(datastr,separators);
            while(~isempty(remainder))
              idtrue = findstr(upper(token),'T');
              idfals = findstr(upper(token),'F');
              if ~isempty(idtrue)
                data0 = 1;
              elseif ~isempty(idfals)
                data0 = 0;
              else
		disp([' %ERROR: data is logical but doesnt contain T or F'])
                data0 = [];
              end
              idx = findstr(token,'*');	% see if multiple values
              if(~isempty(idx))
                 nvals  = str2num(token(1:idx(1)-1));
              else
                 nvals= 1;
              end
              data = [data; logical(data0*ones(nvals,1))];
              [token,remainder] = strtok(remainder,separators);
            end % while
	   else                                 % numeric data
	    idatatype= 1;                       % 1= data, 2= string, 3=T/F
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
            end % while
	   end  % if idx_log, else
         end    % if idx, else
         iend = istart + size(data,1)-1;
         varc= var;
         if to_upper==1
           varc=upper(varc);
	 elseif to_upper==-1
	   varc=lower(varc);
	 end
% We have varc, data, and idatatype here can now output to structure or 'caller'
         if nargout<=1
          if idatatype==1 % Data
	    dataa= num2str(data',16);
            str= [varc '(' int2str(istart) ':' int2str(iend) ',:) =' ...
	        ' [' dataa ']'';'];
 	    evalin('caller',str);
          else          %  Character
            for ii= istart:iend
               dataa= double(data(ii-istart+1,:));
	       dataa= num2str(dataa);
               str= [varc '(' int2str(ii) ',:) = char([' dataa ']);'];
	       evalin('caller',str);
            end
          end
	 else
           str= [varc '(' int2str(istart) ':' int2str(iend) ',:) = data' ';'];
           str= ['nlstr.' nlname '.' str];
           eval(str);
	 end
%         str= [nlname '_list= strvcat(' nlname '_list,''' varc ''');'];
%	 eval(str); % store list of all variables in namelist_list
         nlall= strvcat(nlall, varc); % store all names 
      end
      nlall= strvcat(nlall, nlterm); % terminator for namelist
      if(~isempty(findstr(line1,'$')))
         fprintf('End processing %s\n',line1)
      end
   end
  end

  if ~isempty(nlall)
    if nargout<=0
      assignin('caller','nlall',nlall);
      assignin('caller','nlnames',nlnames);
    elseif nargout>=2
      nlstr.filename= filename;
      nlstr.nlall= nlall;
      nlstr.nlnames= nlnames;
    end
  end

  fclose(fid);

  return

% =============================================
% testing
  clear
  filename= 'namelist_in'
   read_namelist;
  [nlnames1,nlstr1,nlall1]= read_namelist(filename);
  struct_to_ws(nlstr1);
  struct_to_ws(IN_1);
  struct_to_ws(IN_2);
  struct_to_ws(IN_3);


  read_namelist;

  filename ='/u/leuer/efit/east/run_ref1/m010003.02013'
  read_namelist(filename);

  filename ='/u/leuer/efit/east/run_ref1/m010003.02013'  
  [nlnames1,nlstr1,nlall1]= read_namelist(filename);

  [names, isstruc]= struct_to_ws(nlstr1);
   if any(isstruc)                   % Do 2nd level
      id= find(isstruc);
      for ii= 1:length(id)
         str= ['[names1,isstruc1]= struct_to_ws(' char(names(id(ii))) ');'];
         eval(str)
      end  % for
   end    % if any
  

  for ii=1:size(nlnames1,1)
      struct_to_ws(['nlstr1.' remove_space(nlnames1(ii,:),0)]);
  end
