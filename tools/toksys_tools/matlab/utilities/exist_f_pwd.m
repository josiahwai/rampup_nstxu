  function [iopt]= exist_f_pwd(name)

% performs an "exist(name,'file')" only on the present working directory pwd
%
% SYNTAX:
%  [iopt]= exist_f_pwd(name);
%
% INPUT:
%         name	= string of file name to search for in the pwd
%
% OUTPUT:
%         iopt  = same output as exist for file 0 no exist 2 file

% Jim Leuer 1-99

  if nargin <= 0
     disp('%Error: exist_f_pwd takes a single string argument, output =  -99')
     help exist_f_pwd
     iopt=-99;
     return
  end

  if ~isstr(name)
     disp('%Caution: exist_f_pwd needs STRING argument, output = -99')
     help exist_f_pwd
     iopt=-99;
     return
  end

% construct full name

  name_new= [pwd,'/',name];

  iopt= exist(name_new,'file');

  return

