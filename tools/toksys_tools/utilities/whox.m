function ww = whox(varargin)

% WHOX List current variables as WHOS but also expand structures
%    to see all fields inside
%    WHOX by itself will list all variables with expanded structures
%    WHOX with one or more arguments will list those variables that
%    match at least one of the arguments. If an argument contains 
%    points, each 'field' of the argument must match the field of
%    the variable.
%
%    Examples:
%    whox *.*         % List only fields of structures
%    whox *.*.*       % List only fields of fields of structures
%    whox *.*.*.EC*   % Find fields that begin with EC 3 levels down
%
%    See also WHOS.

%  VERSION @(#)whox.m	1.2 11/30/13
%
%  WRITTEN BY:  Anders Welander ON 2013-11-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% All variables from whos returned in ww
ww = evalin('caller','whos');


% Expand ww all the way down every structure
n = 0;
while n < length(ww)
  n = n+1;
  nam = char(ww(n).name);
  if strcmp(ww(n).class,'struct')
    nam1 = strrep(nam,'.','(1).');
    try
      fs = evalin('caller',['fieldnames(' nam1 ')']);
      nfs = length(fs);
      ww((n+1:end)+nfs) = ww(n+1:end);
      for j = 1:nfs
	fnam = char(fs{j});
	ww(n+j).name = [nam '.' fnam];
	try
	  dum = evalin('caller',['getfield(' nam1 ',''' fnam ''');']);
	  ww(n+j) = whos('dum');
	  ww(n+j).name = [nam '.' fnam];
	catch
	end
      end
    catch
    end
  end
end


% Include names that are returned by 'whos' with arguments applied
if length(varargin) > 0
  w = 'whos(';
  for i = 1:length(varargin)
    w = [w,'''',char(varargin(i)),''','];
  end
  w = [w(1:end-1),')'];
  w1 = evalin('caller',w);
  ii = zeros(1,n); % Flags for selected elements of ww
  % Select those names in ww that can be found in w1
  j = 1;
  for i = 1:length(w1)
    while ~strcmp(ww(j).name,w1(i).name) & j < n
      j = j+1;
    end
    if strcmp(ww(j).name,w1(i).name)
      ii(j) = 1; % Was returned by whos call with arguments, so return
    end
  end
else % No arguments in call, so all of ww should be selected
  ii = ones(1,n);
end


% Include names that match arguments with points and wildcards
for m = 1:n
  if ii(n) == 0 % Test for names that aren't already selected
    nam = ww(m).name;
    i1 = [0, findstr(nam,'.'), length(nam)+1];
    for j = 1:length(varargin)
      flag = 1; % Keep if all fields in arg match those in nam
      arg = char(varargin(j));
      i2 = [0, findstr(arg,'.'), length(arg)+1];
      if length(i1) < length(i2)
	flag = 0;
      else
	for j2 = 1:length(i2)-1
          nam1 = nam(i1(j2)+1:i1(j2+1)-1);
          arg1 = regexptranslate('wildcard',arg(i2(j2)+1:i2(j2+1)-1));
	  k1 = regexp(nam1,arg1,'start');
	  k2 = regexp(nam1,arg1,'end');
	  if isempty(k1) | k1(1) ~= 1 | k2(1) ~= length(nam1)
	    flag = 0;
	  end
	end
      end
      if flag
	ii(m) = 1;
      end
    end
  end
end


% Perform the down-selection
ww = ww(ii==1);


if nargout == 0 % Report to the screen if nargout == 0
  if isempty(ww)
    clear ww
    return
  end
  sp = '';
  nm = '';
  s1 = '';
  s2 = '';
  cl = '';
  by = '';
  at = '';
  for i = 1:length(ww)
    nam = char(ww(i).name);
    siz = ww(i).size;
    byt = ww(i).bytes;
    cls = ww(i).class;
    att = '';
    if ww(i).global
      if isempty(att)
	att = [att 'global'];
      else
	att = [att ', global'];
      end
    end
    if ww(i).persistent
      if isempty(att)
	att = [att 'persistent'];
      else
	att = [att ', persistent'];
      end
    end
    if ww(i).sparse
      if isempty(att)
	att = [att 'sparse'];
      else
	att = [att ', sparse'];
      end
    end
    if ww(i).complex
      if isempty(att)
	att = [att 'complex'];
      else
	att = [att ', complex'];
      end
    end
    sp = char(sp,'  ');
    nm = char(nm,nam);
    s1 = char(s1,sprintf('%12s',num2str(siz(1))));
    ss = '';
    for j = 2:length(siz)
      ss = [ss,'x',num2str(siz(j))];
    end
    s2 = char(s2,ss);
    cl = char(cl,[' ',cls,' ']);
    by = char(by,[' ',sprintf('%17d',byt),' ']);
    at = char(at,[' ',att,' ']);
  end
  while sum(s1(:,1) == 32) == size(s1,1)
    s1 = s1(:,2:end);
  end
  sp = char(' ',sp);
  nm = char('Name',nm);
  s1 = char(' ',s1);
  s1(1,end) = 'S';
  s2 = char('ize',s2);
  by = char(' ',by);
  by(1,end-5:end) = 'Bytes ';
  cl = char(' Class',cl);
  at = char(' Attributes',at);
  disp([sp,nm,sp,sp,sp,s1,s2,sp,by,cl,sp,at]);
  disp(' ')
  clear ww
end
