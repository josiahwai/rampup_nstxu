
% save_mat.m: Script to save list of variables defined in var_list
%             into a Matlab binary file (.mat format) defined by SAVE_FILE.
%             Routine will identify variables defined in var_list which are not
%             in environment and will continue saving without these variables.
%             var_list can be list of variables or file .ptn with list of
%             variables.
%             Use LOAD to restore variables from saved file.
%
% SYNTAX:
%        save_mat
%
% INPUT:  [default]
%        save_file=   name of save file to create [save.mat]
%        var_list=    list of variables to save or name of .ptn file with list
%                     1st checks if var_list is a file, 2nd as variable list
%                     if var_list doesn exist it saves all variables
%
%        Example as list:
%                    var_list = ['ncc     '; 'variable1'; 'variable2'];
%        Example as file:
%                    var_list = 'variable.ptn';
%                      where variable.ptn file is:
%                         % comment line
%                         ncc
%                         varable1
%                         % comment line not copied
%                         variable2
%                
% OUTPUT:
%
%        list_save= list of variables saved in file SAVE_FILE
%
%        Creates file with name defined by SAVE_FILE [save.mat] containing
%        all variable defined in VAR_LIST. Those that dont exist aren't saved
%        Saves ALL VARIABLES if VAR_LIST doesnt exist.
%
%
% NOTE: To read Save set use:
%        load(save_file);
%        load save.mat;
%


% Jim Leuer GA 7-8-03
% taken from d3d_sim_save.m
% ============================================================================

 if exist('save_file')~= 1
    disp('%Note: save_mat saving variables to save.mat')
    save_file= 'save.mat';
 end

 if exist('var_list')~= 1 
    disp(['%Note: save_mat saving all variables in environment to: ',save_file])
    dum= whos;
    list_save= {};
    [list_save{1:length(dum),1}] = deal(dum.name);
    list_save= char(list_save);
    var_list=[];
  elseif exist(var_list(1,:))==2 % its a file 
    list_save= read_point_file(var_list(1,:));
  else % assume its a variable list
    list_save= var_list;
  end

  if isempty(list_save)
    disp(['%ERROR save_mat: no variable saved ',save_file])
    return 
  end
   
  dum= whos;
  list_all= {};
  [list_all{1:length(dum),1}] = deal(dum.name);
  list_all= char(list_all);
  id= strsmatch(list_save,list_all,'exact'); 
  idd= find(id>=1); % non-existant symbols are 0, good symbols >=1
  list_save= list_all(id(idd),:); % good symbols  

  if ~isempty(list_save)
  
    str= ['save ',save_file];
    for ii=1:length(list_save(:,1))
     str= [str ' ' deblank(list_save(ii,:))];
    end
  
    eval(str);
    disp(['% save_file just saved ',int2str(length(list_save(:,1))),...
          ' variables in file: ',save_file])
  else
    disp(['%ERROR save_mat Couldnt find any variables to save ',var_list])
  end 

  return

% ===============================================================
% test:
  save_mat

  var1= 1;
  var2= 2;
  list= ['var1';'var2'];
  var_list= 'var.ptn';
  num= write_ptn_file(list,var_list)
  save_mat

  var_list=list
  save_mat
  
  
  
