  function id = strsmatchs(strin,strbase,exact,idebug)
%
% PURPOSE: STRSMATCHS finds pointers in character array STRBASE for array STRIN.
%
%          This is like STRSMATCH except it will return all matches in strbase
%          for a particular row character in strin. If more than one entry is
%          found in strbase then the output has entries in strbase.
%
% SYNTAX:
%        id = strsmatchs(strin,strbase)         % match 1st characters of strin
%        id = strsmatchs(strin,strbase,'exact') % exact match of characters
%
% INPUT: [default]
%        strin=   array of strings to find index in strbase
%        strbase= array of strings to search for strin
%        exact=   'exact' then strings must be exact, {optional} []
%        idebug=  1; turn on warning messages [0]
%
% Output:
%        id=    pointer in strbase array for strings in strin
%               id = -1 problem in arguments
%               0 number indicates element of strin not found in strbase
%
% NOTE: Previous versions had 'exact' by default. New version does not use
%       'exact' as default YOU MUST INCLUDE IT AS AN INPUT VARIABLE
% NOTE: Now strin & strbase can be cell strings or string arrays
%
%  SEE: strmatch, /users/leuer/matlab/util/strsmatch

% Jim Leuer 4-14-99
% 1-20-03 exact added (note old routine had exact as default)

% ---------------------------------------------------------------------------
  if (nargin <= 1)
     disp('%ERROR STRSMATCHS needs 2 arguments')
     help strsmatchs
     id= -1;
     return
  end

  if (nargin <= 2)
     exact=[];
     idebug= 0;
  end

  if (nargin <= 3)
     idebug= 0;
  end

  if iscellstr(strin)  % convert "cellstr" to "char_array"
     strin= char(strin);
  end
  if iscellstr(strbase) % convert "cellstr" to "char_array"
     strbase= char(strbase);
  end
  
  if (~isstr(strin) | ~isstr(strbase))
     disp('%ERROR STRSMATCH: Both arguments must be strings')
     help strsmatch
     id= -1;
     return
  end

  strin_leng= length(strin(:,1));

 if strcmp(exact,'exact') % exact string compare  
  for ii=1:strin_leng
    str= remove_space(strin(ii,:));
    idd= strmatch(str,strbase,'exact');
    if isempty(idd)
      if idebug 
         disp(['%WARNING strsmatch: couldnt strin in strbase: ',str]);
      end
      id(ii,1)= 0;
    else
      lidd= length(idd);
      id(ii,1:lidd)= idd';
    end
  end 
 else     % match only 1st characters of string
  for ii=1:strin_leng
    str= remove_space(strin(ii,:));
    idd= strmatch(str,strbase);
    if isempty(idd)
      if idebug 
         disp(['%WARNING strsmatch: couldnt strin in strbase: ',str]);
      end
      id(ii,1)= 0;
    else
      lidd= length(idd);
      id(ii,1:lidd)= idd';
    end
  end 
 end

  if length(id(1,:))>=2
    if idebug 
     disp(['%NOTE strsmatchs: Output has more than one column: ',...
              int2str(length(id(1,:)))]);
    end
  end
  
  return

% test stuff

%  strin= ['str1  ';...
%          'str2  ';...
%          'str34 ';...
%          'str   '];

%  strbase= ['str  ';...
%            'str0 ';...
%            'str1 ';...
%            'str2 ';...
%            'str3 ';...
%            'str2 ';...
%            'str22';...
%            'str4 '];

% old:
%  id0= strsmatch(strin,strbase,'exact')
%  id1= strsmatch(strin,strbase)

%  id= strsmatchs(strin,strbase,'exact')
%  id2= strsmatchs(strin,strbase)
%  id3= strsmatchs(strin,strbase,[],1)
                                    
