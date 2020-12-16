 function [jtext,ipatt]= strsfind(text,patterns);
 %
%  SYNTAX:
%          [jtext,ipatt] = strsfind(text,patterns);
%       
%  PURPOSE:  Multiple "patterns" version of "strfind"
%            Finds patterns (1 or more rows) in single text (single row)
%            NOTE: trailing blanks in patterns are removed
%
%  INPUT: 
%         text=	    single row of text (ex: '/home/leuer/tokamaks/nstx/efit/')
%         patterns= patterns to search "text" (ex: ['east'; 'nstx'; 'iter'...])
%                   NOTE: trailing blanks in patterns(ii,:) are removed
%
%  OUTPUT:
%         jtext=      index in text(j) of start of pattern
%         ipatt=      row of patterns assoicated with a each particular jtext id
%              note: jtext and ipatt have same length
%
%  RESTRICTIONS:
%         text and patterns must be type character; text is a single row 
 
%  METHOD: 
%            executes strfind(text,patterns(ii,:)) ii= 1,... 
%
%  WRITTEN BY:  Jim Leuer 	ON 	19Mar2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if nargin<=1 %
    disp('%ERROR strsfind: needs 2 input arguments (text,patterns)');
    return
  end

  if size(text,1)>=2 | ~ischar(text)
    disp('%ERROR strsfind: TEXT input must be single charcter row');
    return
  end
  if ~ischar(patterns)
    disp('%ERROR strsfind: PATTERNS input must be charcter array');
    return
  end

% Start loop
  jtext= [];
  ipatt= [];
  num= size(patterns,1);  
  for ii=1:num
    dum= strfind(text,deblank(patterns(ii,:)));
    jtext= [jtext dum];
    ipatt= [ipatt ii*ones(1,length(dum))];
  end

  return

% =============================================
% testing
 
%  text= pwd;
%  patterns= strvcat('d3d','east','nstx','iter')
%  [itext]= strsfind(text,patterns)
%  [itext, jpatt]= strsfind(text,patterns)
   itext= strsfind(pwd,patterns)
  
  
