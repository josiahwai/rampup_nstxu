  function string_ot= remove_space(string, position)
 % 
% SYNTAX: 
%         string_ot= remove_space(string_in, position);
%
% PURPOSE: Removes spaces in string (like deblank but allows different location removal)
%
% INPUT: [default]
%       string_in=   input string (single row string NOT string array)
%       position=    0 all spaces, 1=front, [2=END], 3=front&end {OPTIONAL)
%
% OUTPUT:
%       string_ot= output string with spaces removed
%
%  RESTRICTIONS: Only tested on 2-d arrays (not on vectors or scalars)
 
%  WRITTEN BY:  jim leuer 4-19-02
% jal 12/19/2006 modified to handle null inputs and single space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)remove_space.m	1.3 02/24/10

   string_ot = string;
   if(iscell(string))
      if(max(size(string))>1)
         wait('ERROR remove_space: can only handle single cell element')
         return;
      end
      return_cell=1;
      string_in = string{:};
   else
      return_cell=0;
      string_in = string;
   end

  if nargin <=0
    help remove_space
    return
  end

  string_ot= string_in;

  if isempty(string_in) % do nothing if empty
     return
  end

  if isstr(string_in) == 0
      disp(['%ERROR remove_space.m: input string_in is not a string: ', ...
      num2str(string_in)]);
      return; 
  end
 
  if nargin == 1
    position= 2;
  end
  
  [row,col]=size(string_in);
  
  if row ~= 1 | col <= 0
    disp('%CAUTION remove_space.m: input string_in should be [1,:]');
  else 

    id1= isspace(string_in); % find space, tab, linefeed, etc
    id2= double(string_in)==0; % uninitialized array are: double(' ')=0 
    id= ~(id1+id2); % this is all characters that are not zero
  
    if position==0
      string_ot= string_in(find(id));
    elseif position==1
      id= find(id);
      if ~isempty(id)
        string_ot= string_in(id(1):end);
      else
        string_ot= char([]);
      end
    elseif position==2
      id= find(id);
      if ~isempty(id)
        string_ot= string_in(1:id(end));
      else
        string_ot= char([]);
      end
    elseif position==3
      id= find(id);
      if ~isempty(id)
        string_ot= string_in(id(1):id(end));
      else
        string_ot= char([]);
      end
    else
      disp('%CAUTION remove_space.m: position must be 0,1,2,3');
    end
  end

  if(return_cell)
      temp = string_ot;
      string_ot = cell(1);
      string_ot{1} = temp;
  end

return

% old below
 
  for i=1:col
   im= col-i+1;
   if string_in(1,im) ~= ' ', break, end
   if im == 1, string_ot=[],  break, end
    string_ot= string_in(1,1:im-1);
  end   

return
