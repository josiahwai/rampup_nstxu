  function [names, str_names]= struct_to_ws_all(s,toupper,putwhere)
%
%  PURPOSE:  Places all variables found in structure in WS
%            works through all substructures to find ALL variables
%            (no structures copied to WS only variables)
%
%  SYNTAX:
%         [names, str_names]= struct_to_ws_all(struc); %
%
%  INPUT: <default>
%    s        = Any Matlab Structure
%    toupper  = 1;    %make all variables UPPER CASE, -1=lower case, <0>=no change
%    putwhere = 'base'; % puts variables in base WS; <'caller'>=in calling routine
%
%  OUTPUT:
%     names	name of all variables made in environment (ex. PSIRZ) 
%     str_names	full name of all substructures found in struc (ex. eq.GEQDSK)
%
%  NOTE: No structures copied to WS (only variables are copied)
%  NOTE: Uses Recursive execution of struc_names to get all variables and structures
%
%  SEE: struct_to_ws for only 1st level (variables and structures copy to environment)

%
%  WRITTEN BY:  Jim Leuer    ON      22Apr2008
% ==========================================================================
% 
% Defaults
  if nargin==0
     disp('% STRUCT_TO_WS_ALL needs at least one argument of type "structure"')
     help struct_to_ws_all
     names=[];
     return
  elseif nargin==1
     toupper= 0;
     putwhere= 'caller';
  elseif nargin==2
     putwhere= 'caller';
  end

  if ~isstruct(s)
      disp('% STRUCT_TO_WS_ALL argument must be of type "structure" ')
      return
  end

  names= [];
% -----------------------------------------------
% Get list of all Variables and Structures
% -----------------------------------------------
  [var_names, str_names]= struct_names(s);
  if ~isempty(var_names)
     id= strfind(var_names(1,:),'.');
     vn= [ char(('s'*ones(size(var_names,1),1)))  var_names(:,id(1):end)];
     for ii=1:size(vn,1)
         id= strfind(vn(ii,:),'.');
	 nam= remove_space(vn(ii,id(end)+1:end));
	 if toupper==1
	    nam= upper(nam);
	 elseif toupper==-1
	    nam= lower(nam);
	 end
         eval(['data =' vn(ii,:) ';']);
         assignin(putwhere,nam,data);
	 names= strvcat(names,nam);
     end
  else
     disp(['% STRUCT_TO_WS_ALL found no variables  in:' inputname(1)])
     names= [];
     str_names= [];
  end
       
  return
 
% =============================

% Testing 
   st1= struct('v1',1,'v2','var2'); 
   st2= struct('v3',3,'v4','var4','st1',st1);
   st3= struct('v5',5,'v6','var6','st2',st2);
   [var_names, str_names]=struc_names(st3)
   [names, str_names]= struct_to_ws_all(st3)
   [names, str_names]= struct_to_ws_all(st3,1,'base');

   addpath /home/leuer/matlab/efit/
   shot= 131498;
   tree='EFIT01';
   server= 'DIII-D';
   eq= eq_mds(shot);
   [var_names, str_names]=struc_names(eq)
   [names, str_names]= struct_to_ws_all(eq,[],'base');
   [names, str_names]= struct_to_ws_all(eq,-1,'base');
