 function tokamak= tok_from_pwd(tokamak)
%
% PURPOSE: determines "tokamak" from your "present working directory" pwd
%          sets global variable TOKAMAK 
%          runs startup_"tok".m routine for particular machine
%
% SYNTAX:  
%          tokamak= tok_from_pwd;
%                   tok_from_pwd
%          tokamak= tok_from_pwd(tokamak); % if you know tokamak to work on
%          tokamak= tok_from_pwd([]);      % return tokamak but no run of startup
% INPUT:   
%          tokamak= OPTIONAL name of tokamak is known [default= looks in pwd]
%
% OUTPUT:
%          tokamak= tokamak name based on search of pwd
%          runs appropriate startup routine as long as tokamak not empty
%
% EXAMPLE:
%          pwd= '/home/leuer/tokamaks/d3d/d3d_sim'
%          tokamak= tok_from_pwd;
%  ==>     tokamak= 'd3d'; TOKAMAK= 'd3d'; Routine runs: ...startup_d3d.m
%

% Jim Leuer 20Oct2008 from my startup.m
% ==============================================================================

  if nargin >=1
     tokin= tokamak;
  else
     tokin= 'none';
  end

  gatools_root= getenv('GATOOLS_ROOT');
  if isempty(gatools_root)
     disp('% ERROR: tok_from_pwd environmental variable GATOOLS_ROOT doesnt exist')
     if exist('tokamak')~=1 tokamak= char([]); end
     return
  end
% ==============================================================================
% establish tokamak based on comparison of all tokamaks with name found in pwd
% ==============================================================================
% Note: Below only works if you have "no startup.m" file in your local directory

  if exist('tokamak')~=1 | isempty(tokamak) % check to see if already exist
    nam= dir([gatools_root '/tokamaks/']);
    if ~isempty(nam)
      if isfield(nam,'name')
        nam= char(nam.name);
        id= strmatch('.',nam); idd= ones(length(nam(:,1)),1); idd(id)= 0;
        tokamaks= nam(find(idd),:);
        clear id idd
      end
    end 
% Debug
%    pwd, 
%    tokamaks
    [itext,itok]= strsfind(pwd,tokamaks); % find if any tokamaks exist in pwd
    if ~isempty(itok)
     tokamak= deblank(tokamaks(itok(1),:));
    else
     disp(['% Caution: tok_from_pwd couldnt find tokamak in pwd: ' pwd])
     tokamak= char([]);     
    end
  end

% ==============================================================================
% if variable 'tokamak' exists => read approprate startup file and set TOKAMAK
% ==============================================================================

 if ~isempty(tokin)
  if exist('tokamak')==1 & ~isempty(tokamak) % tokamak is variable in environment
     disp(['% Setting global variable: TOKAMAK = ' tokamak])
     global TOKAMAK
     TOKAMAK= tokamak;
     dum= [gatools_root '/tokamaks/' tokamak '/' tokamak '_startup.m'];
     dum2=[gatools_root '/startups/'  tokamak '_startup.m'];
     if exist(dum)==2
        disp(['% Running: ' dum])
        run(dum)
     elseif exist(dum2)==2
        disp(['% Running: ' dum])
        run(dum2)        
     end
  end
 end % if ~isempty('tokin')
 return
