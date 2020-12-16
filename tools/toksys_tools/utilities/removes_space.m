  function string_ot= removes_space(string_in, position)

% removes spaces in sring (like deblank but allows different location removal)
%
% SYNTAX: 
%         string_ot= removes_space(string_in, position);
%
% INPUT: [default]
%       string_in=   input string
%       position=    0 all spaces, 1=front, [2=END], 3=front&end {OPTIONAL)
%
%       note: this recursively runs remove_space on each line of string_in
%
% OUTPUT:
%       string_ot= output string with spaces removed
%
% Note: This works on columns

% New algorithm for more flexibility
% jim leuer 4-19-02


  if nargin <=0
    help remove_space
    return
  elseif nargin == 1
    if isstr(string_in) == 0
      disp(['%ERROR remove_space.m: input string_in is not a string: ', ...
      num2str(string_in)]);
      return; 
    end
    position= 2;
  end
  
  string_ot=string_in;

  if isempty(string_in)
    disp('%CAUTION remove_space.m: input string_in is empty');
    string_ot=string_in;
    return
  end
  
  [row,col]=size(string_in);

  string_ot= [];
  for ii=1:row
    str= string_in(ii,:);
    id1= isspace(str); % find space, tab, linefeed, etc
    id2= double(str)==0; % uninitialized array are: double(' ')=0 
    id= ~(id1+id2); % this is all characters that are not zero
  
    if position==0
      str= str(find(id));
    elseif position==1
      id= find(id);
      str= str(id(1):end);
    elseif position==2
      id= find(id);
      str= str(1:id(end));
    elseif position==3
      id= find(id);
      str= str(id(1):id(end));
    else
      disp('%CAUTION remove_space.m: position must be 0,1,2,3');
    end
    string_ot= strvcat(string_ot,str);
  end

  return

