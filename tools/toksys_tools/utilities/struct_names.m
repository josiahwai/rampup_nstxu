  function [var_names, str_names]=struct_names(struc, var_names, str_names, pre_name)
%
%  PURPOSE:  Returns names of all variables and substructures within structure
%
%  SYNTAX:
%         [var_names, str_names]= struct_names(struc); % Normal Execution
%
%  INPUT: <default>
%     struc	structure object to find internal variable and structure names
%
%  OUTPUT:
%     var_names	full name of all variables in structure (ex. eq.GEQDSK.PSIRZ) 
%     str_names	full name of all substructures (ex. eq.GEQDSK)
%
%  NOTE: TO GET ALL VARIABLES IN STRUCTURE TO WORKSPACE see struct_to_ws_all
%
%  NOTE: Uses Recursive execution of struct_names to get all variables and structures
%        in recursive operation is includes additional inputs:
%        [var_names, str_names]= struct_names(struc, var_names, str_names, pre_name);
%        where: 
%        var_names	input var_names to append output
%        str_names	input str_names to append output 
%        pre_name	is name of previous structure to prepend on str_names
%

%
%  WRITTEN BY:  Jim Leuer    ON      22Apr2008
% ==========================================================================
% 
% Defaults
  if nargin==0
     disp('% STRUCT_NAMES needs at least one argument of type "structure"')
     help struct_names
     var_names= [];
     str_names= [];
     return
  elseif nargin==1
     var_names= [];
     str_names= [];
     pre_name=  [];
  elseif nargin==2
     str_names= [];
     pre_name=  [];
  elseif nargin==3
     pre_name=  [];
  end

  if ~isstruct(struc)
      disp('% STRUCT_NAMES argument must be of type "structure" ')
      return
  end

  struc_name= inputname(1);
% -----------------------------------------------
% Read Names
% -----------------------------------------------

  namest= fieldnames(struc); % structure field names    (cell_str)
  for ii=1:length(namest) % min(3,length(namest))
     str= ['struc' '.' char(namest(ii)) ];
     str2= [pre_name struc_name '.' char(namest(ii)) ];
     eval(['iss= isstruct(' str ');']);
     if iss         % structure
       str_names= strvcat(str_names,str2);
       strr= ['[var_names, str_names]= struct_names(' str ...
               ' ,var_names, str_names, str2 );'];
       eval(strr)   % variable    
     else
       var_names= strvcat(var_names,str2);
     end
  end
       
   return
 
% =============================

% Testing 

   addpath /home/leuer/matlab/efit/
   shot= 131498;
   tree='EFIT01';
   server= 'DIII-D';
   eq= eq_mds(shot);
   [var_names, str_names]=struct_names(eq)
   [var_names, str_names]=struct_names(st)
