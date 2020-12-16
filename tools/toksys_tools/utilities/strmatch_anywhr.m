 function [irow, jcol]= strmatch_anywhr(str,strs,flag);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:
%       FUNCTION CALL:
%           irow=        strmatch_anywhr(str,strs); % rows containing str
%          [irow, jcol]= strmatch_anywhr(str,strs,'all'); % return all matches
%
%       
%  PURPOSE:  Returns rows of strs which have str anywhere in line
%
%  INPUT: [default]
%         str=	    single row of text (ex: 'F2a')
%         strs=     string array to look for str (ex: strvcat('PF1a','PF2a')
%         flag=     'all'; % returns all matches [remove duplicate row entrys]
%
%  OUTPUT:
%         irow=     row index of strs which have str (ex: 2)
%         jcol=     colum index of 1st occurance of strs (ex: 2)
%
%  RESTRICTIONS:
%         str and strs must be type character; str is a single row 
%	
 
%  METHOD: 
%            executes strfind(text,patterns(ii,:)) ii= 1,... 
%
%  WRITTEN BY:  Jim Leuer 	ON 	2May2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if nargin<=1 %
    disp('%ERROR strmatch_anywhr: needs 2 input arguments (str,strs)');
    return
  end

  if size(str,1)>=2 | ~ischar(str)
    disp('%ERROR strmatch_anywhr: STR input must be single charcter row');
    return
  end

  if ~ischar(strs)
    disp('%ERROR strmatch_anywhr: STRS input must be charcter array');
    return
  end

  if nargin<=2 %
    flag= char([]);
  end

  dum= strs';
  dum= dum(:)';
  id= strfind(dum,str)';
  jlen= size(strs,2);
  jcol= mod(id-1,jlen)+1;      % find column J of match
  irow= (id - jcol)/jlen + 1;
  
% default is to remove duplicate rows unless flag=='all' 
  if ~strcmp(lower(flag),'all')
    idd= diff(irow);
    if any(~idd)
      idd= find([1;idd]);
      irow= irow(idd);
      jcol= jcol(idd);
    end
  end
  
  return

% =============================================
% testing
 
   str=  'stx';
   str=  't';
   str=  'x';
   strs= strvcat('d3d','eastt','nstx','iter')
   [irow,icol]= strmatch_anywhr(str,strs)
   [irow,icol]= strmatch_anywhr(str,strs,'ALL')
  
  
