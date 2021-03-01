  function id = strmatch_end(strin,strbase,exact)
%
% PURPOSE: STRSATCH finds STRIN at END of STRBASE and returns pointer 
%
%          This is like strmatch except it looks at end of each row for strin
%          blanks at end of array are ignored. Good for finding things like:
%          ZBBBS in charter array containing row with: eq.RESULTS.GEQDSK.ZBBBS
%
% SYNTAX:
%        id = strmatch_end(strin,strbase)         % match 1st characters of strin
%        id = strmatch_end(strin,strbase,'exact') % exact match , no warnings
%
% INPUT: [default]
%        strin=   array of strings to find index in strbase
%        strbase= array of strings to search for strin
%        exact=   'exact' then strings must be exact, {optional} []
%
% Output:
%        id=    pointer in strbase array for strings in strin
%               id = -1 problem in arguments
%               [] number indicates element of strin not found in strbase
%
% NOTE: Previous versions had 'exact' by default. New version does not use
%       'exact' as default YOU MUST INCLUDE IT AS AN INPUT VARIABLE
% NOTE: Now strin & strbase can be cell strings or string arrays

% Jim Leuer 4-14-99
% 1-20-03 exact added (note old routine had exact as default)

% ---------------------------------------------------------------------------

  if (nargin <= 1)
     disp('%ERROR STRMATCH_END needs 2 arguments')
     help strsmatch_end
     id= -1;
     return
  elseif nargin <= 2
     exact=[];
  end

  if iscellstr(strin)  % convert "cellstr" to "char_array"
     strin= char(strin);
  end
  if iscellstr(strbase) % convert "cellstr" to "char_array"
     strbase= char(strbase);
  end

  if (~isstr(strin) | ~isstr(strbase))
     disp('%ERROR STRSMATCH: Both arguments must be char_array or cellstr')
     help strsmatch
     id= -1;
     return
  end

  strin_leng= length(strin(:,1));
  id= zeros(strin_leng,1);

% use reverse order to do check
 if strcmp(exact,'exact') % exact string compare  
   id= strmatch(removes_space(strin(end:-1:1),3), ...
       removes_space(strbase(:,end:-1:1),3),'exact'); 
 else
   id= strmatch(removes_space(strin(end:-1:1),3), ...
       removes_space(strbase(:,end:-1:1),3)); 
 end

  
  return

% test stuff

%  strin= ['str22  '];

%  strbase= ['123 str  ';...
%            '123 str0 ';...
%            '123 str1 ';...
%            '123 str2 ';...
%            '123 str3 ';...
%            '123 str2 ';...
%            '123 str22';...
%            '321 str4 '];
%
%  id= strmatch_end(strin,strbase)
%  id= strsmatch(strin,strbase,'exact')                                 

%  strin= cellstr(strin); % test for use of cellstr
%  strbase= cellstr(strbase); % test for use of cellstr
