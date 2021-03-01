function [eq, neq, eqraw] = read_eq(shot,times,source,tokamak,server)
%
%  USAGE: eq = read_eq(shot,times,source,tokamak)
%
%  PURPOSE: Read equilibrium data from many possible sources and make structure
%           according to toksys convention. Implemented sources are:
%           1. mds
%           2. A directory with efit g-files
%           3. A directory with corsica flat files
%
%  INPUTS:  shot is the shot number (no surprise there)
%           times: times to fetch [sec]
%             If empty or 'all' then all equilibria are fetched
%             If single time specified, closest available time fetched
%             If two times are specified they are tmin and tmax
%             If >2 times are specified, only those exact times are fetched
%           source: mds tree or directory with g files or directory with corsica flatfiles
%             If source begins with a '/' it is a directory name, otherwise mds tree
%             If server is 'corsica', source should contain flatfiles by create_corsica_flat_files
%           tokamak: name of the tokamak, default 'DIII-D'
%	    server: optional mds server name. Default is standard mds server for tokamak.
%             Use 'local' for data from GA test server (currently soas).
%
%  OUTPUTS:  eq = equilibrium data on toksys form, that can be used with
%              cc_efit_to_tok and build_tokamak_system
%            neq = number of equilibria returned
%            eqraw = 'raw' eq data on original form before conversion to toksys convention
%
%  RESTRICTIONS: Needs modification to work with Corsica files for tokamaks other than ITER
%
%  METHOD:
	
%  WRITTEN BY:  Anders Welander  ON	6/15/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if nargin < 4
    tokamak='DIII-D';
    disp(['read_eq: tokamak was set to the default: ' tokamak])
  end
  if nargin < 3
    source = 'EFIT01';
    disp(['read_eq: source was set to the default: ' source])
  end
  if nargin < 2
    times='ALL';
    disp(['read_eq: times was set to the default: ' times])
  end
  if nargin < 1
    shot = 106652;
    disp(['WARNING read_eq: No inputs were specified! shot was set to a default of : ' num2str(shot)])
    beep, pause(0.25), beep, pause(0.25), beep, pause(0.25), beep
  end
  eqraw = [];
  
  
  % CORSICA FLAT FILES
  if strcmp(upper(tokamak),'CORSICA') | (exist('server','var') & strcmp(upper(server),'CORSICA'))
    neq = 0;
    if source(end)~='/', source(end+1)='/'; end
    [s,w] = unix(['ls -l ' source '*/cc.flat']);
    w(end+1)=10;
    k=1;
    while k<length(w)-1
      while k<length(w)-1 & w(k)~='/', k=k+1; end
      k1 = k;
      while k<length(w)-1 & w(k)~=10, k=k+1; end
      dirname = strrep(w(k1:k-1),'cc.flat','');
      neq = neq+1;
      if ~exist('eq','var')
        q=read_corsica_flat_files(dirname);
        eq.gdata=q;
      else
        q=read_corsica_flat_files(dirname);
        eq.gdata(end+1)=q;
      end
    end
    
  % A&G FILES FROM A DIRECTORY
  elseif source(1) == '/'
    if source(end) ~= '/', source(end+1) = '/'; end
    [s,w] = unix(['ls -l ' source]);
    shotstr = num2str(shot);
    k = findstr(w,['g' shotstr '.']);
    for j=1:6, shotstr = ['0' shotstr]; k=[k findstr(w,['g' shotstr '.'])]; end % Find even if shot begins with zeros
    gtimestrs = '';
    gtimes = [];
    for j=1:length(k)
      j1 = min(findstr(w(k(j):end),'.'));
      j2 = min(find(w(k(j)+1:end)==10));
      s = w(k(j)+j1:k(j)+j2);
      if ~isempty(str2num(strrep(s,'_','.')))
        gtimestrs = strvcat(gtimestrs,s);
        gtimes(end+1) = str2num(strrep(s,'_','.')); %gtimes in ms
      end
    end
    if exist('gtimes') ~= 1
      disp(['WARNING read_eq: No gfile with shot# ' num2str(shot) ' found in ' source])
      eq = [];
      return
    end
    [gtimes, k] = sort(gtimes);
    gtimestrs = gtimestrs(k,:);
    if length(times) == 0
      times = gtimes;
      timesstr = gtimestrs;
    elseif length(times) == 1
      [dum, k] = min(abs(times-gtimes*1e-3));
      times = gtimes(k); timesstr = gtimestrs(k,:);
    elseif length(times) == 2
      k = find(gtimes*1e-3>=min(times) & gtimes*1e-3<=max(times));
      times = gtimes(k); 
      timesstr = gtimestrs(k,:);
    else
      k = [];
      for j=1:length(gtimes)
	if find(times == gtimes(j)*1e-3)
	  k(end+1) = j;
	end
      end
      times = gtimes(k);  %now times in ms
      timesstr = gtimestrs(k,:);
    end
    eq.shotnum = shot;
    eq.time = times(:)/1000;
    shotstr = num2str(shot); while length(shotstr)<6, shotstr = ['0' shotstr]; end
    neq = length(times);
    for j = 1:neq
      eq.gdata(j) = read_gfile_tok([source 'g' shotstr '.' deblank(timesstr(j,:))],tokamak);
    end
    for j = 1:neq
      eq.gdata(j).time = times(j)/1000;
    end
    if nargout > 2
      eqraw = eq;
    end
    
  % MDS
  else
    if length(times) == 1 | length(times) == 2 % return the 1 nearest time or a time range
      if exist('server','var')
        [eq, neq, eqraw] = read_mds_eqdsk(shot, sort(times), source, tokamak, server);
      else
        [eq, neq, eqraw] = read_mds_eqdsk(shot, sort(times), source, tokamak);
      end
    elseif isempty(times) | strcmp(upper(times),'ALL') % return all times
      if exist('server','var')
        [eq, neq, eqraw] = read_mds_eqdsk(shot, 'ALL', source, tokamak, server);
      else
        [eq, neq, eqraw] = read_mds_eqdsk(shot, 'ALL', source, tokamak);
      end
    else % return specific times
      if exist('server','var')
        [eqdata, neq, eqraw] = read_mds_eqdsk(shot,[-1e9 1e9],source,tokamak,server);
      else
        [eqdata, neq, eqraw] = read_mds_eqdsk(shot,[-1e9 1e9],source,tokamak);
      end
      k = [];
      for j=1:length(eqdata.time)
	if length(find(times == eqdata.time(j))) | isempty(times)
	  k(end+1) = j;
	end
      end
      eq.shotnum = eqdata.shotnum;
      eq.time = eqdata.time(k);
      eq.gdata = eqdata.gdata(k);
      eq.adata = eqdata.adata(k);
      eq.descriptions = eqdata.descriptions;
      neq = length(eq.gdata);
    end
  end
  
