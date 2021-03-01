  function output= whoss(varargin);

% WHOSS is like whos except it prints out actual variable value instead of
%       size information. It prints out first couple values of an array
%
% SYNTAX: whoss (structure and usage similar to whos)
%         whoss  a* *name* ...
%         out= whoss(var1,var2,...)
%
%         doesnt work for structures

%            

 
% Jim Leuer 8-13-02 Leuer@fusion.gat.com

% put together input arguments to whos
  if isempty(varargin)
    w= 'whos';
  else
    w= 'whos(';
    for ii= 1:length(varargin)
      w = [w,'''',char(varargin(ii)),''','];
    end
    w= [w(1:end-1),')'];
  end

  ww = evalin('caller',w);

% Check if any variables are found
  if isempty(ww)
     error('%Caution - whoss: reports NO variables Found')
     return
  end

  vmax= 20;  % maximum number of array variable to display
  cmax= 60;  % maximum number of characters to display
  lmax= 80;  % maximum length of a line
    
  nm=  char([]);
  sz=  char([]);
  cl=  char([]);
  val= char([]);

  for ii=1:length(ww)
    nam= char(ww(ii).name);
    siz= ww(ii).size;
    byt= ww(ii).bytes;
    cls= ww(ii).class;
        
    nm=  strvcat(nm,nam);
    sz=  strvcat(sz,[' ',int2str(siz(1)),'x',int2str(siz(2))]);
    cl=  strvcat(cl,[' ',cls,' ']);

    if byt ~= 0
     switch cls
      case 'double'
        len= min(prod(siz),vmax);
        str= [nam,'(1:',int2str(len),')'];
        value= evalin('caller',str);
	value= reshape(value,1,prod(size(value)));
        val=  strvcat(val,num2str(value));	
      case 'char'
        len= min(siz(1,2),cmax);
        str= [nam,'(1,1:',int2str(len),')'];
        value= evalin('caller',str);
        val=  strvcat(val,value);	
      case 'cell'
        str= [nam,'(1);'];
        value= char(evalin('caller',str));
        len= min(length(value),cmax);
        value= value(1:len);
        val=  strvcat(val,value);	
      case 'struct'
        str= ['fieldnames(',nam,');'];
        value= evalin('caller',str);
        len= min(length(value(1,:)),cmax);
        value= char(value(1,1:len));
        val=  strvcat(val,value);	
      otherwise
        val=  strvcat(val,'unknown');
     end
    else   
     val=  strvcat(val,'empty');
    end

  end
     
% output

  if nargout<=0
    dum= [nm,sz,cl,val];
    dum= dum(:,1:min([lmax, size(dum,2)]));
    disp(dum);
%    disp([nm,sz,cl,val]);
  else
    output= struct('name',nm,'size',sz,'class',cl,'val',val);
  end
  
  return
% one fix in char, remove ";" for work on PC jal. See old below
%        str= [nam,'(1,1:',int2str(len),');'];
% jal26/mar/2007  add lmax
