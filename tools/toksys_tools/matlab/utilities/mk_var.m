  function numvar= mk_var(variable, value);

% mk_var sets variable = value if it does not exist in the calling workspace
%        same as:
%                if exist('variable') variable = value; end
%
% SYNTAX:
%         numvar= mk_var(variable,value);
%
% INPUT:
%         variable= name of variable to make (IF IT DOESNT ALREADY EXIST)
%         value=    value of variable Default: [];
%            
% OUTPUT:
%         makes "variable=value" in calling workspace if it doesnt exist
%
%         numvar = number of variables made (0 if no variables made -1 if error)
%            
% LIMITATIONS:
%         variable can be an string array containing a list of vairables
%         however value then must be a single variable, or an array with
%         same row length as variable. If not it sets all variables  to value.
 
% Jim Leuer 12-04-02 Leuer@fusion.gat.com

% ============================
% check input arguments

  if nargin<=0
    disp('%ERROR mk_var needs at least one input variable')
    help mk_var;
    numvar=-1;
    return
  end
  if nargin==1
     value= []*ones(length(variable(:,1)),1);
  end

  if ~isstr(variable)
    disp('%ERROR mk_var: 1st argument must be a string or string array')
    output=-1;
    return
  end

  numvar= 0;

  if length(variable(:,1)) == 1
    str= ['exist(' '''' deblank(variable(1,:)) ''')'];
    chk = evalin('caller',str);
    if chk~=1
      assignin('caller',deblank(variable(1,:)), value)
      numvar= 1;
    end
  elseif length(variable(:,1))==  length(value(:,1))
   for ii= 1:length(variable(:,1))
    str= ['exist(' '''' deblank(variable(ii,:)) ''')'];
    chk = evalin('caller',str);
    if chk~=1
      assignin('caller',deblank(variable(ii,:)), value(ii,:))
      numvar= numvar+1;
    end    
   end
  else
    disp('%CAUTION mk_var: length variable(:,1) and value(:,1) not same')
    for ii= 1:length(variable(:,1))
      str= ['exist(' '''' deblank(variable(ii,:)) ''')'];
      chk = evalin('caller',str);
      if chk~=1
       assignin('caller',deblank(variable(ii,:)), value)
       numvar= numvar+1;
      end
    end   
  end

  return
  
