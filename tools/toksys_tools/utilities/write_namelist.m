 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX: 
%        write_namelist % runs like script with variableread from calling routn.
%         
%
%  PURPOSE:  write FORTRAN namelist file . Typically used after a read_namelist
%            has been read in and data changed and then written out.
% 
%
%  INPUT: [default]
%    nlall=        namelist names and variables (typical from read_namelist)
%    filename_ot = name of namelist file to written
%    to_upper    = 1; % Variables: 1-UPPER CASE, -1=LOWER CASE, [0=no change]
%    n_per_line  =    number of numeric data to put on single line [5]
%
%  OUTPUT:
%    All variables found in namelist(s) within the file are made into
%    matlab variables. Duplicate names are overwritten with subsequent data.
%    nlnames=       Character Array contaning names of each namelist read
%    NAMELIST_list= list of each variable in NAMELIST
%                   NAMELIST is name of each item in nlnames
%                   ex. if file contains &IN then IN_list contains variable nms.
%    nlall=         total list of all namelists and variable (write_namelist)
%
%  RESTRICTIONS:
%
 
%  METHOD:  
%
%  WRITTEN BY:  Jim Leuer 	ON 	12/04/06
%  taken in part from old write_namelist and new Mike read_namelist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if exist('idebug')     ~= 1      idebug= 1;     end
  if exist('n_per_line') ~= 1      n_per_line= 5; end
  if exist('to_upper')   ~= 1      to_upper=   0; end


  fid = fopen(filename_ot,'w');

  icnt = 0;
  linecnt = 0;
  separator = char(' ');       % seperator for output 'space
  

  line1=0;
  nlcnt=0;
  nlnames = '';
% =======================================
  icnt= 0;
  nlines= size(nlall,1);
    
  while icnt < nlines
    icnt= icnt + 1;
    line1 = nlall(icnt,:);

%    if line1==-1 break, end
   
% find start of namelist

   totstr = '';
   if(~isempty(findstr(line1,'$')) | ~isempty(findstr(line1,'&')))
      fprintf('Begin Namelist processing %s\n',line1)

      fprintf(fid,'%s\n',deblank(line1)); % **** OUTPUT NAMELIST NAME TO FILE

      nlcnt =   nlcnt+1;
      nlstart=  deblank(line1); % start full name list name with & or $
      nlname=   remove_space(line1,0);
      nlsym=    nlname(1);     % $ or & symbol used to start name list name
      nlname=   nlname(2:end); % name
% find end of namelist
      while icnt <= nlines
         icnt= icnt + 1;
         line1 = nlall(icnt,:);
         idx = union(findstr(line1,'$'),findstr(line1,'&')); % problem for str with &
         if(~isempty(idx))
           break
         end;
         if(length(line1)>1 & strcmp(line1(2),'/')) % break on back slash
           break
         end;
         totstr = strvcat(totstr,line1);
      end
   end

   nlend= deblank(line1); % termination line of namelist

% If a namelist is found and is nonempty, then parse the string in between.

   if ~isempty(totstr)

      ll = size(totstr,1);

      for kk=1:ll			     % loop over all variables

        var = remove_space(totstr(kk,:),0); % variable 
        istart = 1;			     % starting index

%        find if type =1 number, 2= logical, 3= char array, All others are ERROR
        type= -1;
	str= [    'if isnumeric(' var ') type= 1; ' ...
              'elseif islogical(' var ') type= 2; ' ...
	      'elseif ischar('    var ') type= 3; ' ...
	      'else   type= 0; '...
	      'disp([''%ERROR Dont know var type: '' ' '''' var '''' ']); ' ...
	      'end'];
	eval(str);
        str= [' val= ' var ';'];  % set "val" to values of variable
	eval(str)
        [ni,nj]= size(val);
	nval= ni*nj;              % total number of variables (or characters)

        if nval > 0 % data includes info (if ni*nj==0 null set dont print)	 
 	    
         varc= var;
         if to_upper==1
           varc=upper(varc);
	 elseif to_upper==-1
	   varc=lower(varc);
	 end
         
%        ===============================
         if type==3			% CHARACTER ARRAY
 
           fprintf(fid,' %s= \n',varc); % ***** Print out variable 'name = <cr>' 
 
           for ii=1:ni;
               fprintf(fid,' ''%s''\n',val(ii,:)); % ***** Print one string/line
           end
   
         elseif type==2		        % LOGICAL

           str= char('F'*ones(1,nval)); % array of false
	   id= find(val);
	   str(id)= char('T'*ones(1, length(id))); % overwrite True states
	   strd= char('.'*ones(1,nval));
	   strs= char(separator*ones(1,nval));
	   strr= [strd; str; strd; strs]; % see strr(:)'

           if nval <= 3*n_per_line % print on one line:	   
             fprintf(fid,' %s= %s\n',varc,strr(:)'); % *****  name = .T. .F. ,,,
           else
             fprintf(fid,' %s=',varc);
             for ii= 1:n_per_line:nval
                i_end= min(nval,ii+n_per_line-1);
                fprintf(fid,' %s\n',val(ii:i_end));    %     T. .F. ,,,
             end
	   end

         elseif type==1		        % NUMERIC

           if nval <= n_per_line % print on one line:	   
             fprintf(fid,' %s= ',varc); % *****  name = 
             fprintf(fid,' %16.10g',val(:)); % *****  data
             fprintf(fid,'\n'); % return
           else
             fprintf(fid,' %s=',var);
             for ii= 1:n_per_line:nval
                i_end= min(nval,ii+n_per_line-1);
                fprintf(fid,' %16.10g',val(ii:i_end));    %     T. .F. ,,,
                fprintf(fid,'\n'); % return
             end
	   end

         else
%	   disp('ERROR NO TYPE')
	 end    % if type   
 	end     % if nval > 0     
       end      % for kk

       fprintf(fid,'%s\n',nlend); % ***** Print Ending of NL & <cr>' 
       if(~isempty(findstr(line1,'$')))
         fprintf('End processing %s\n',line1)
       end
    end         % ~isempty(totstr)
end             % while icnt <= nlines


fclose(fid);
%clear separators sep_idx token remainder totstr value var varc nvals
%clear mlen ll linecnt line1 k idx icnt i temp str

return

% testing
  clear
  filename= 'namelist_in'
  read_namelist
  
  str1='this is a new string'
  num1= [1 2 3];
  
  filename_ot= 'nl_ot.dat'
  write_namelist
  
 fclose(fid) 
  write_namelist
