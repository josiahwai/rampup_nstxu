  function [names, isstruc]=struct_to_ws(struc,toupper,putwhere,verbose,list,sym)
%
%  PURPOSE:  Pulls 1st level field names of structure into calling routine work
%            space or optionally the main workspace.
%            Arbitrarly loads only those symbols contained in list
%
%  SYNTAX:
%          [names,isstruc]= struct_to_ws(struc);
%          [names,isstruc]= struct_to_ws(struc, toupper, putwhere, verbose);
%          [names,isstruc]= struct_to_ws(struc,toupper,putwhere,verbose,list,sym);
%
%  INPUT: [default] 
%    struc    = Any Matlab Structure
%    toupper  = 1; make all variables UPPER CASE, -1=lower case, [0]=no change
%    putwhere = 'base'; puts variables in base WS; ['caller']=in calling routine
%    verbose  = 1; prints out all variables created [0]=no print
%
%    list     = list of variables to put in environment []= all 1st lev. variables
%    sym      = symbol(s) to append to end of variable name []
%               Note: if symapp(1) =='+' then the remaining symapp(2:end) is
%                     "pre-appended" to symbol name instead of "post-appended"
%
%  OUTPUT:
%           All variables in struc.* are placed into workspace "caller" or "base"
%
%           names=   list of names of all variables placed in workspace
%           isstruc= column vector of length names: 1= struct 0= variable
%
%                 Example Structure:	eq.AFILE.BCENTR(1,5), eq.NW
%                 Function Call:	names= struct_to_ws(eq)
%                 Into Workspace:	AFILE.BCENTR(1,5), NW
%
%  NOTE: Can be used 2nd, 3rd, ... time to get all variables into workspace
%
%  Example of putting in variables 1&2 levels down:           
%        [names, isstruc]= struct_to_ws(eq); % Do 1st level
%        if any(isstruc)                   % Do 2nd level
%          id= find(isstruc);
%          for ii= 1:length(id)
%            str= ['[names1,isstruc1]= struct_to_ws(' char(names(id(ii))) ');'];
%            eval(str)
%          end  % for
%        end    % if any

%
%  WRITTEN BY:  Jim Leuer    ON      12/14/05
% ==========================================================================
% 
% Defaults

  if nargin==0
     disp('% STRUCT_TO_WS needs at least one argument of type "structure"')
     help struct_to_ws
     names=[];
     return
  elseif nargin==1
     toupper= 0;
     putwhere= 'caller';
     verbose= 0;
     list= [];
     sym=  [];
  elseif nargin==2
     putwhere= 'caller';
     verbose= 0;
     list= [];
     sym=  [];
  elseif nargin==3
     verbose= 0;
     list= [];
     sym=  [];
  elseif nargin==4
     list= [];
     sym=  [];
  elseif nargin==5
     sym=  [];
  end

  if ~isstruct(struc)
      disp('% STRUCT_TO_WS argument must be of type "structure" ')
      return
  end

  if isempty(toupper)
     toupper= 0;
  end
  if isempty(putwhere)
     putwhere= 'caller';
  end
  if isempty(verbose)
     verbose= 0;
  end

  if ~strcmp(putwhere,'base')  & ~strcmp(putwhere,'caller')
     disp('% STRUCT_TO_WS argument: putwhere must be ''base'' or ''caller'' ')
     return
  end

  sym_add= 0; % sym_add=0 none, sym_add=1= post, sym_add=-1 pre-add
  if ~isempty(sym)
     sym_add= 1;
     if sym(1)=='+'
        sym_add=-1;
	sym= sym(2:end);
	if isempty(sym) sym_add=0; end
     end
  end
    
% -----------------------------------------------
% Read Names
% -----------------------------------------------
 isstruc= [];
 namest= fieldnames(struc); % structure field names    (cell_str)
 nameso= namest; % Output names with toupper correction (cell_str)

% If 'list' exists deselect names 
  if ~isempty(list); % deselect names if list is present
    id= strsmatch(list,nameso);
    idd= find(id);
    nameso= nameso(id(idd));
    namest= nameso;
    iddd= find(id==0);
    if ~isempty(iddd)
      disp('% CAUTION: struct_to_ws found strings in "list" but not in "struc"')
    end  
  end
    
  names= char([]); % final output names
  if ~isempty(nameso)
   isstruc= zeros(length(nameso),1);
   if toupper==0
   elseif toupper==1
       nameso= upper(nameso);
   elseif toupper==-1
       nameso= lower(nameso);
   end  
   for ii=1:length(nameso)
       namest1= char(namest(ii));
       nameso1=  char(nameso(ii));
       data= getfield(struc,namest1);
       if isstruct(data)
          isstruc(ii)= 1;
       end
       if sym_add
          if sym_add==1
	     nameso1= [nameso1 sym];
	  elseif sym_add==-1
	     nameso1= [deblank(sym) nameso1];
	  end
       end
       assignin(putwhere,nameso1,data);
       names= strvcat(names,nameso1);
       if verbose
          if isstruc(ii)
	    disp(['% STRUCT_TO_WS put STRUCTURE in workspace: ',nameso1])
	  else
            disp(['% STRUCT_TO_WS put variable in workspace: ',nameso1])
	  end
       end
    end % for
   end   % ~isempty
 
   return
 
% =============================

% Testing 

   addpath /thor/leuer/matlab/efitgold_new/
   shot= 113363;
   tree='EFIT01';
   server= 'NSTX';
   eq= eq_mds(shot, tree,  server);
   names= struct_to_ws
   names= struct_to_ws(eq);
   [names, isstruc]= struct_to_ws(eq,0);
   [names, isstruc]= struct_to_ws(eq,-1,'caller');
   names= struct_to_ws(RESULTS,0,'base',1);
   names= struct_to_ws(GEQDSK,1,1);
   [names, isstruc]= struct_to_ws(GEQDSK,-1);
% test of list and sym:
   list= fieldnames(AEQDSK)
   names= struct_to_ws(AEQDSK,0,'base',1,list(1:4));
   names= struct_to_ws(AEQDSK,0,'base',1,list(5:8),'_1')
   
