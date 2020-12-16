function write_gfile(eq)
%  USAGE:  write_gfile(eq)
%
%  PURPOSE: write gfiles
%
%  INPUTS:  eq: equilibrium on the toksys format (that is returned by read_mds_eqdsk)
%
%  OUTPUTS: gfiles to disk
%
%  RESTRICTIONS: 
%

%  VERSION @(#)write_gfile.m	1.3 09/17/12
%
%  WRITTEN BY:  Anders Welander  ON	6/1/11
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  xlim = 0; ylim = 0; % fixes matlab bug
  
  struct_to_ws(eq);
  
  if ~exist('pasmat','var')
    pasmat = cpasma; % pasmat is perhaps measured plasma current
  end

  if ~exist('ecurrt','var')
    ecurrt = [];
  end

  tims = round(1000000*time)/1000;
  ss = num2str(shotnum); while length(ss)<6, ss = ['0' ss]; end
  tt = num2str(fix(tims)); while length(tt)<5, tt = ['0' tt]; end
  if mod(tims,1), tt = [tt '_' num2str(mod(tims,1))]; end
  filename = ['g' ss '.' tt];
  
  fid = fopen(filename,'w');
  
  fprintf(fid,'  EFITD    ');
  fprintf(fid,datestr(now,23));
  
  str = num2str(shotnum);
  while length(str)<6, str = [' ' str]; end
  str = ['    #' str];
  fprintf(fid,str);
  
  str = [num2str(tims) 'ms'];
  while length(str)<8, str = [' ' str]; end  
  str = [str blanks(8)];
  fprintf(fid,str);
  
  fprintf(fid,'%4d%4d%4d\n',3,nw,nh);
  
  if ~exist('xdim','var')
    xdim = rg(end)-rg(1);
  end
  if ~exist('zdim','var')
    zdim = zg(end)-zg(1);
  end
  if ~exist('zmid','var')
    zmid = (zg(end)+zg(1))/2;
  end
  fprintf(fid,'% 11.9E',xdim);
  fprintf(fid,'% 11.9E',zdim);
  fprintf(fid,'% 11.9E',rzero);
  fprintf(fid,'% 11.9E',rg(1));
  fprintf(fid,'% 11.9E\n',zmid);

  fprintf(fid,'% 11.9E',rmaxis);
  fprintf(fid,'% 11.9E',zmaxis);
  fprintf(fid,'% 11.9E',ssimag);
  fprintf(fid,'% 11.9E',ssibry);
  fprintf(fid,'% 11.9E\n',bzero); % bcentr = bzero?
  
  xdum = 0; % introducing xdum...
  
  fprintf(fid,'% 11.9E',cpasma);
  fprintf(fid,'% 11.9E',ssimag);
  fprintf(fid,'% 11.9E',xdum);
  fprintf(fid,'% 11.9E',rmaxis);
  fprintf(fid,'% 11.9E\n',xdum);
  
  fprintf(fid,'% 11.9E',zmaxis);
  fprintf(fid,'% 11.9E',xdum);
  fprintf(fid,'% 11.9E',ssibry);
  fprintf(fid,'% 11.9E',xdum);
  fprintf(fid,'% 11.9E\n',xdum);
  
  for j = 1:nw
    fprintf(fid,'% 11.9E',fpol(j));
    if mod(j,5) == 0, fprintf(fid,'\n'); end
  end
  if mod(j,5), fprintf(fid,'\n'); end
  for j = 1:nw
    fprintf(fid,'% 11.9E',pres(j));
    if mod(j,5) == 0, fprintf(fid,'\n'); end
  end
  if mod(j,5), fprintf(fid,'\n'); end
  for j = 1:nw
    fprintf(fid,'% 11.9E',-sign(pasmat)*ffprim(j));
    if mod(j,5) == 0, fprintf(fid,'\n'); end
  end
  if mod(j,5), fprintf(fid,'\n'); end
  for j = 1:nw
    fprintf(fid,'% 11.9E',-sign(pasmat)*pprime(j));
    if mod(j,5) == 0, fprintf(fid,'\n'); end
  end
  if mod(j,5), fprintf(fid,'\n'); end
  for j = 1:nh
    for k = 1:nw
      fprintf(fid,'% 11.9E',psirz(k,j));
      if mod((j-1)*nw+k,5) == 0, fprintf(fid,'\n'); end
    end
  end
  if mod(j,5), fprintf(fid,'\n'); end
  for j = 1:nw
    fprintf(fid,'% 11.9E',qpsi(j));
    if mod(j,5) == 0, fprintf(fid,'\n'); end
  end
  if mod(j,5), fprintf(fid,'\n'); end
  
  if ~exist('limitr','var')
    limitr = length(xlim);
  end  
  fprintf(fid,'% 5d',nbbbs);
  fprintf(fid,'% 5d\n',limitr);
  for j = 1:nbbbs
    fprintf(fid,'% 11.9E',rbbbs(j));
    if mod(2*j-1,5) == 0, fprintf(fid,'\n'); end
    fprintf(fid,'% 11.9E',zbbbs(j));
    if mod(2*j,5) == 0, fprintf(fid,'\n'); end
  end
  if mod(2*nbbbs,5)~=0
    fprintf(fid,'\n');
  end
  for j = 1:limitr
    fprintf(fid,'% 11.9E',xlim(j));
    if mod(2*j-1,5) == 0, fprintf(fid,'\n'); end
    fprintf(fid,'% 11.9E',ylim(j));
    if mod(2*j,5) == 0, fprintf(fid,'\n'); end
  end
  
  % Rotation information
  if ~exist('kvtor','var')
    kvtor = 0;
  end  
  if ~exist('rvtor','var')
    rvtor = 1.7;
  end  
  if ~exist('nmass','var')
    nmass = 0;
  end  
  if mod(2*limitr,5)~=0
    fprintf(fid,'\n');
  end

  fprintf(fid,'% 5d',kvtor);
  fprintf(fid,'% 11.9E',rvtor);
  fprintf(fid,'% 5d\n',nmass);
  
  % write out rotation information
  if kvtor>0
    for j = 1:nw
      fprintf(fid,'% 11.9E',pressw(j));
      if mod(j,5) == 0, fprintf(fid,'\n'); end
    end
    fprintf(fid,'\n');
    for j = 1:nw
      fprintf(fid,'% 11.9E',-sign(pasmat)*pwprim(j));
      if mod(j,5) == 0, fprintf(fid,'\n'); end
    end
    fprintf(fid,'\n');
  end
  
  % write out ion mass density profile if available
  if nmass>0
    for j = 1:nw
      fprintf(fid,'% 11.9E',dmion(j));
      if mod(j,5) == 0, fprintf(fid,'\n'); end
    end    
    fprintf(fid,'\n');
  end  
  
  if ~exist('rhovn','var')
    rhovn(1:nw) = 0;
  end  
  for j = 1:nw
    fprintf(fid,'% 11.9E',rhovn(j));
    if mod(j,5) == 0, fprintf(fid,'\n'); end
  end
  if mod(j,5), fprintf(fid,'\n'); end
  
  if ~exist('keecur','var')
    keecur = 0;
  end
  fprintf(fid,'% 5d\n',keecur);
  if keecur>0
    for j = 1:nw
      fprintf(fid,'% 11.9E',-sign(pasmat)*epoten(j));
      if mod(j,5) == 0, fprintf(fid,'\n'); end
    end
    if mod(j,5), fprintf(fid,'\n'); end
  end
  
  if exist('jphi','var') % if true then iplcout was 1
    fprintf(fid,'% 5d',nw);
    fprintf(fid,'% 5d',nh);
    fprintf(fid,'% 5d',shotnum);
    fprintf(fid,'% 5d\n',tims);

    fprintf(fid,'% 11.9E',rg(1));
    fprintf(fid,'% 11.9E',rg(end));
    fprintf(fid,'% 11.9E',zg(1));
    fprintf(fid,'% 11.9E\n',zg(end));
    
    if ~exist('nfcoil','var')
      nfcoil = length(brsp);
    end
    for j = 1:nfcoil
      fprintf(fid,'% 11.9E',brsp(j));
      if mod(j,5) == 0, fprintf(fid,'\n'); end
    end
    if mod(j,5), fprintf(fid,'\n'); end
    
    if ~exist('nesum','var')
      nesum = length(ecurrt);
    end
    for j = 1:nesum
      fprintf(fid,'% 11.9E',ecurrt(j));
      if mod(j,5) == 0, fprintf(fid,'\n'); end
    end
    if mod(j,5), fprintf(fid,'\n'); end
    
    for j = 1:nw*nh
      fprintf(fid,'% 11.9E',pcurrt(j)); % Need to check r,z order
      if mod(j,5) == 0, fprintf(fid,'\n'); end
    end
    if mod(j,5), fprintf(fid,'\n'); end
  end
  
  out1{ 1,1} = 'ishot';
  out1{ 1,2} = 'shot';
  out1{ 2,1} = 'itime';
  out1{ 2,2} = 'tims';
  out1{ 3,1} = 'betap0';
  out1{ 3,2} = 'betap';
  out1{ 4,1} = 'btor';
  out1{ 4,2} = 'bzero';
  out1{ 5,1} = 'rcentr';
  out1{ 5,2} = 'rzero';
  out1{ 6,1} = 'rbdry';
  out1{ 6,2} = 'rbbbs';
  out1{ 7,1} = 'zbdry';
  out1{ 7,2} = 'zbbbs';
  out1{ 8,1} = 'nbdry';
  out1{ 8,2} = 'nbbbs';
  out1{ 9,1} = 'mw';
  out1{ 9,2} = 'nw';
  out1{ 9,3} = 'nr';
  out1{10,1} = 'mh';
  out1{10,2} = 'nh';
  out1{10,3} = 'nz';
  out1{11,1} = 'psirz';
  out1{12,1} = 'xlim';
  out1{13,1} = 'ylim';
  out1{14,1} = 'limitr';
  out1{15,1} = 'brsp';
  out1{15,2} = 'cc';
  out1{15,3} = 'ic';

  fprintf(fid,[' &OUT1' 10]);
  for i = 1:size(out1,1)
    S = upper(out1{i,1});
    for j = 1:size(out1,2)
      s = out1{i,j};
      if ~isempty(s) & exist(s,'var')
        dum = eval(s);
	n = min(size(dum));
	for k = 1:numel(dum)
	  if k == 1
	    V = [' ' S ' = '];
	  else
	    V = 32+zeros(1,length(S)+3);
	  end
	  if dum(k) > 0
	    x = ' ';
	  else
	    x = '';
	  end
	  m = num2str(dum(k));
	  if isa(dum,'logical')
	    if dum(k)
	      m = 'T';
	    else
	      m = 'F';
	    end
	  end
	  if k/n == round(k/n) % Add new line at the end
	    if (k-1)/n == round((k-1)/n) % Add V at beginning of line
	      fprintf(fid,[V x m ',' 10]);
	    else
	      fprintf(fid,[x m ', ' 10]);
	    end
	  else
	    if (k-1)/n == round((k-1)/n) % Add V at beginning of line
	      fprintf(fid,[V x m ',']);
	    else
	      fprintf(fid,[x m ', ']);
	    end
	  end
	end
	break % Don't look for alternatives
      end
    end
  end
  fprintf(fid,[' /' 10]);
  
  fclose(fid);
