
  function strot = strcombine(str1,str2)
%
% PURPOSE: strcombine finds single string from combination of 1st identical
%          characters of str1 & str2 (or all rows of strs array).
%          ex: str1= 'id1 ', str2= 'id2'  => strot='id'
%
% SYNTAX:
%        strot = strcombine(str1,str2); % 2 single row character array
%        strot = strcombine(str);       % multi row character array
%
% INPUT: [default]
%        str1, str2=  single row character array
%        strs=        multi row character array
%
% Output:
%        strot=  single row string containing identical starting characters
%
% Note: Trailing blanks are not considered characters
%       function is case sensitive so use upper/lower to get identical case

% Jim Leuer 04Apr'07

% ---------------------------------------------------------------------------

  if nargin == 0
     help strcombine
     strot=char([]);
     return
  elseif nargin == 1
     strs= str1;
  else
     strs= strvcat(str1,str2);
  end

  if ~isstr(strs)
     disp('%WARNING STRCOMBINE: input arguments should be strings')
  end

  [nrow,ncol]= size(strs);
  
  if nrow <= 1  
    strot= deblank(strs);
  else
    [ir,ic]= find(diff(deblank(strs))~=0);
    id= min(ic)-1;
    if id>=1
     strot= strs(1,1:id);
    else
     strot=char([]);
    end    
  end
	 

  return
    
% testing
%  format compact
%  strs= ['pf1  ';'pf2  ';'pf3 1';'pf1  '];  
%  strot= strcombine(strs)
%  str1= strs(1,:); str2= strs(1:2,:);
%  strot= strcombine(str1,str2)
    
