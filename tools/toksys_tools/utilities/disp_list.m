  function strot= disp_list(strin)
%
% PURPOSE: DISP_LIST displays list of character arrays with prepended row number
%                    primarily for display of row number/string identification
%
% SYNTAX:
%                disp_list(strin); % displays strot
%        strot = disp_list(strin); % outputs array for future display

% INPUT: [default]
%        strin=   array of strings to display
%
% Output: (if no output given it displays str_array)
%       strot  =   strin with row id in first column(s) followed by space
%
% RESTRICTIONS: strin must be a character array or cell array

% Jim Leuer 18May07

% ---------------------------------------------------------------------------

  if (nargin <= 0)
     disp('%ERROR disp_list needs 1 argument')
     help disp_list
     str_array= char([]);
     return
   end

  if iscellstr(strin)  % convert "cellstr" to "char_array"
     strin= char(strin);
  end

  if ~isstr(strin) 
     disp('%ERROR disp_list Argument must be char_array or cellstr')
     str_array= char([]);
     return
  end

  len= size(strin,1);
  sp= char(' '*ones(1,len))';
  strot= [int2str((1:len)') sp strin];
  if nargout==0  
    disp(strot);
  end

  return

% test stuff

  strin= ['str1  ';...
          'str2  ';...
          'str34 ';...
          'str   '];
%  id= strsmatch(strin,strbase)
%  id= strsmatch(strin,strbase,[],1)                                 
%  id= strsmatch(strin,strbase,'exact',0)                                 
