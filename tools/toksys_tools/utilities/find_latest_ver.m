  function [name_latest,version_n]= find_latest_ver(name,iopt)

% find_latest_ver.m finds latest version of a routine
% assumes name is in format: xxxxxx54.m where xxxxxx is program name and
% 54 is version. It does a directory on xxxxxx*.m and finds largest number
% also will work on xxxxx.m* finding xxxxx.m45. 
%
% SYNTAX:
%  [name_latest,version_n]= find_latest_ver(name);
%
% INPUT:
%         name	= string of file name to do directory on
%                 (ex: program; program.m; program*.m)
%         iopt  = 0= return array of all matches not just latest
%                 1= return full name [default]
%                 2= eliminate extention .* 
%                 3= eliminate path */
%                 4= eliminate extention and path
%
%         if ipt negative (ie. -1) then returns all file names not just latest
%
% OUTPUT:
%         name_latest = latest version of file
%                       (ex: program106.m)
%         version_n   = version number (ex. 106)
%
% Caution: name without a number is considered version -1
%          name0 is considered version 0;
%          name01 is same as name1
%          this doesnt care what names wild card * gets only looks for 5555.
%
% NEW Search: If it doesnt find routine in current directory it looks for file
%             without the wild card in the path using which(). If it finds the
%             the file then it does a directory in this "BASE" area 

% Jim Leuer 12-98
% jal 4-23-02 include search of path using which if not in local directory 

  if nargin <= 0
     disp('%Caution: find_latest_version takes a single string argument')
     help find_latest_ver
     name_latest= []; version_n= -99;
     return
  end

  if nargin <= 1
     iopt= 1;
  end

  if ~isstr(name)
     disp('%Caution: find_latest_version needs STRING argument')
     help find_latest_ver
     name_latest= []; version_n= -99;
     return
  end

% parse name

  nname= name;
% eliminate below to make more general
%  if isempty(findstr(name,'*')) % if it contains * then use name as is
%     id= findstr(name,'.');     % if it contains a . add *.
%     if isempty(id)
%        nname= [name,'*'];   % add * to end if no . present
%     else           
%        nname= [name(1:(id(1)-1)),'*',name(id(1):end)];  
%     end
%  end

  d= dir(nname);      % structure of directory
  dc= struct2cell(d); % cell array

  if isempty(dc)
% try looking in path using which 
    id = findstr(nname,'*'); % get rid of wild card *
    idd= ones(size(nname));
    idd(id)=0;
    idd= find(idd);
    nnname= nname(idd);
    dnam= which(nnname); % look in path for name
    if isunix
       ch='/';
    else
       ch='\';
    end
    if isempty(dnam)
     id=[];
    else
     id=  find(dnam==ch);
    end
    if isempty(id)
      dr=[];
    else
     dr= dnam(1:id(end));
    end
    d= dir([dr nname]); % do directory with path name from which
    dc= struct2cell(d); % cell array

    if isempty(dc)
     disp(['%Caution: couldnt find any files using: dir(',nname,')'])
     name_latest= []; version_n= -99;
     return
    end
  end

  ds= char(dc(1,:));  % column string of directories 

% Find Version Number: version number for each entry: num
  num= zeros(length(ds(:,1)),1);
  for ii= 1:length(ds(:,1));
    dss= ds(ii,:);
    id= findstr(dss,'.');   % find '.' in name
    if isempty(id);
       id= length(dss) + 1;
    end
    ide= id - 1;
    dn= char(dss(1:ide));
    idd= find(dn>=48 & dn<=57);
    if isempty(idd)
       num(ii)= -1; %no number present so this is -1
    elseif idd(end) ~= ide
       num(ii)= -1; % no number before . or end of string
    else
      ids= ide;
      for iii= 1:length(idd)
        im= length(idd)-iii+1;
        idg= ide-iii+1;
        if idg ~= idd(im) break, end
        ids= idd(im);
      end
      num(ii)= str2num(dn(ids:ide));
    end
  end

% return list of all matching files if iopt=0

  if iopt == 0
   version_n= num;
   name_latest= ds;
   return
  else
   [version_n,idv]= max(num);
   name_latest= ds(idv,:);
   name_latest= remove_space(name_latest); % remove space on right hand side
  end

  if iopt <= 1 return; end

  if iopt == 2 | iopt == 4
    id= findstr(name_latest,'.');
    if ~isempty(id)
      id= max(1,id(end)-1);
      name_latest=  name_latest(1:id);
    end
  end

  if iopt == 3 | iopt == 4
    id= findstr(name_latest,'/');
    if isempty(id)
       return
    end
    id= min(length(name_latest),id(end)+1);
    name_latest=  name_latest(id:end);
  end

  return

