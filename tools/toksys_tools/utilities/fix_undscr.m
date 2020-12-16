  function [string_ot]= fix_undscr(string_in)

% PURPOSE: put \ infront of any _ so it plots correctly when using TEX
%
% SYNTAX:  [string_ot]= fix_undscr(string_in);
%
% INPUT:
%       string_in=  string containing _ which are converted to \_
%
% OUTPUT:
%       string_ot=  string containing \_ which prints ok in TEX
% 
% Note: Now works on string arrays

% Jim Leuer 1-19-98
% jal 12/15/03 made so works on arrays
% --------------------
% start checking

  if nargin <= 0 
    disp('%ERROR: fix_undscr: Needs atleast 1 input argument')
    help fix_undscr
    string_ot= [];
  end

  if isempty(string_in)
     string_ot= string_in;
     return
  end
  
% --------------------
% Add \ infront of every _

 string_ot= char([]);
 for ii=1:size(string_in,1)
  str= string_in(ii,:);
  id= find(str == '_');
  if isempty(id)
    string_ot= strvcat(string_ot,str);
  else
    is= 1;
    str_ot= char([]);
    for iii= 1:length(id)
      str_ot= [str_ot, str(is:id(iii)-1), '\_'];
      is= id(iii) + 1;
    end
    if is <= length(str)
     str_ot= [str_ot, str(is:end)];
    end
   string_ot= strvcat(string_ot,str_ot);
  end
 end
  return
