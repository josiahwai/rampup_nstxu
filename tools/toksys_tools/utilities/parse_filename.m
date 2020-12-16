  function [dir,nam,ext]= parse_filename(string_in)

% PURPOSE: determine the directory, name and extention of a file
%
% SYNTAX:  [dir,nam,ext]= parse_filename(string_in)
%
% INPUT:
%       string_in=  string containing full file name (ex: /dir/subdir/name.ext)
%
% OUTPUT: (example: /dir/subdir/name.ext)
%       dir=  directory name of file with last /     (ex: /dir/subdir/)
%       nam=  name of file without extention         (ex: name)
%       ext=  extention name with pre .              (ex: .ext)
% 
% CAUTION: Only works for single row vector of ASCII characterics

% Jim Leuer 1-98

% --------------------
% initialization

  if isunix dirsym='/'; else dirsym= '\'; end; % take care of Unix/PC dir names 

% --------------------
% Input checking

  if nargin <= 0 
    disp('%ERROR: parse_filenam: Needs at least 1 input argument')
    help parse_filenam
    dir=[]; nam=[]; ext=[];
    return
  end

  if isempty(string_in)
     dir=[]; nam=[]; ext=[];
     return
  end

  if ischar('string_in')~=1
     disp('%CAUTION: string_in: Input string_in must be a character string')
     dir=[]; nam=[]; ext=[];
     return
  end
  
% ------------------------------------------------------------
% start parsing

% last slash represents the end of the directory:
  id=  find(string_in==dirsym);
  if isempty(id)
     dir=[];
  else
     dir= string_in(1:id(end));
     if id(end) >= length(string_in)
        nam=[]; ext=[];
        return
     end
     string_in= string_in(id(end)+1:end);
  end    

  if isempty(string_in)
     nam= []; ext= [];
     return
  end

% last period represents the extention
  id= find(string_in=='.');
  if isempty(id)
     ext=[];
  else
     ext= string_in(id(end):end);
     if id(end) >= 2
        string_in= string_in(1:id(end)-1);
     else
        string_in= [];
     end
  end    

% what remains is the name
  nam= string_in;

  return
