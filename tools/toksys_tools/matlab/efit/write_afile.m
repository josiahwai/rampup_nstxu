function write_afile(equil,eqversion,aday)
%  USAGE:  write_afile(equil,eqversion,aday)
%
%  PURPOSE: write afiles
%
%  INPUTS:  equil: equilibrium on the toksys format (that is returned by read_mds_eqdsk)
%           eqversion: (optional) source of equilibrium (default is TOKSYS)
%           aday: (optional) the date the equilibrium was created
%
%  OUTPUTS: afiles to disk
%
%  RESTRICTIONS: 
%

%  VERSION @(#)write_afile.m	1.3 10/09/14
%
%  WRITTEN BY:  Anders Welander  ON	6/1/11
%
%  MODIFICATION HISTORY:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  xlim = 0; ylim = 0; % fixes matlab bug
  mu0 = 4e-7*pi;
  
  struct_to_ws(equil);
  if exist('adata','var') & isstruct(adata)
    struct_to_ws(adata);
  end
  if ~exist('shotnum','var')
    shotnum = 0;
  end
  if ~exist('dr','var') & exist('dw','var')
    dr = dw;
  end
  if ~exist('dz','var') & exist('dh','var')
    dz = dh;
  end
  
  if ~exist('uday','var')
    uday = datestr(date,'dd-mmm-yyyy');
    uday = uday([1:7 10:11]);
  end
  while length(uday)<10, uday = [' ' uday]; end
  if ~exist('mfvers','var')
    mfvers = datestr(date,23);
  end
  
  if ~exist('tims','var')
    if exist('time','var')
      tims = round(1000000*time)/1000;
    else
      tims = 0;
    end
  end
  ss = num2str(shotnum); while length(ss)<6, ss = ['0' ss]; end
  tt = num2str(fix(tims)); while length(tt)<5, tt = ['0' tt]; end
  if mod(tims,1), tt = [tt '_' num2str(mod(tims,1))]; end
  filename = ['a' ss '.' tt];
  
  fid = fopen(filename,'w');
  
% 1 -------------------------- Start writing file
  fprintf(fid,'%s %s\n',uday,mfvers);
  
% 2 --------------------------------------------
  if ~exist('ktime','var'), ktime = 1; end % ktime is total number of time slices 
  if length(num2str(shotnum))>5
    fprintf(fid,'%7d',shotnum);
  else
    fprintf(fid,'%6d',shotnum);
  end
  fprintf(fid,'%16d\n',ktime);

% 3 --------------------------------------------
  fprintf(fid,' % 11.9E\n',tims); % time in ms

% 4 --------------------------------------------
  if ~exist('jflag','var'), jflag = 1; end
  if ~exist('lflag','var'), lflag = 0; end % lflag is an efit error flag
  if ~exist('limloc','var'), limloc = 'SNB'; end % Type of shape: IN , OUT, TOP, BOT, DN , SNT, SNB, MAR
  if ~exist('mco2v','var'), mco2v = 3; end % Number of co2v items
  if ~exist('mco2r','var'), mco2r = 2; end % Number of co2r items
  if ~exist('fwtqa','var'), fwtqa = 0; end % Fit weight for assumed q on axis
  if ~exist('qvfit','var'), qvfit = 0; end % A parameter for q-profile model used by fitting logic
  if fwtqa > 0 & qvfit > 0
    qmflag='FIX';
  else
    qmflag='CLC';
  end
  if ~exist('n1old','var'), n1old = 39; end
  if ~exist('n1new','var'), n1new = 40; end
  fprintf(fid,'*%8.3f',tims);
  fprintf(fid,'%14d',jflag);
  fprintf(fid,'%16d',lflag);
  fprintf(fid,'%4s',limloc);
  fprintf(fid,'%4d%4d',mco2v,mco2r);
  fprintf(fid,'%4s',qmflag);
  fprintf(fid,'%6d',n1old);
  fprintf(fid,'%5d',n1new);
  fprintf(fid,'\n ');
  
% 5 --------------------------------------------
  if ~exist('rcentr','var')
    rcentr = rzero;
  end
  if ~exist('tsaisq','var'), tsaisq = 0; end % saisq for 1 or several time slices
  if ~exist('rcencm','var'), rcencm = rzero*100; end % rcentr*100 for 1 or several time slices
  if ~exist('bcentr','var'), bcentr = bzero; end % Bt at rcentr for 1 or several time slices
  if ~exist('pasmat','var'), pasmat = cpasma; end % Measured? plasma current for 1 or several time slices
  skriv(fid,tsaisq);
  skriv(fid,rcencm);
  skriv(fid,bcentr);
  skriv(fid,pasmat);
  fprintf(fid,'\n ');

% 6 --------------------------------------------
  if ~exist('rout','var'), rout = 100*(max(rbbbs(1:nbbbs))+min(rbbbs(1:nbbbs)))/2; end
  if ~exist('zout','var'), zout = 100*(max(zbbbs(1:nbbbs))+min(zbbbs(1:nbbbs)))/2; end
  if ~exist('aout','var'), aout = 100*(max(rbbbs(1:nbbbs))-min(rbbbs(1:nbbbs)))/2; end
  skriv(fid,cpasma);
  skriv(fid,rout);
  skriv(fid,zout);
  skriv(fid,aout);
  fprintf(fid,'\n ');

% 7 --------------------------------------------
  if ~exist('eout','var')
    eout = 100*(max(zbbbs(1:nbbbs))-min(zbbbs(1:nbbbs)))/2/aout;
  end
  if ~exist('doutu','var')
    [ztop, k] = max(zbbbs(1:nbbbs));
    rtop = 100*rbbbs(k);
    doutu = (rout-rtop)/aout;
  end
  if ~exist('doutl','var')
    [zbot, k] = min(zbbbs(1:nbbbs));
    rbot = 100*rbbbs(k);
    doutl = (rout-rbot)/aout;
  end
  if ~exist('vout','var')
    vout = 0;
    for j = 1:nh
      for k = 1:nw
        if jphi(j,k)~=0
	  vout = vout+2*pi*rg(k)*dr*dz*1e6;
        end
      end
    end
  end
  skriv(fid,eout);
  skriv(fid,doutu);
  skriv(fid,doutl);
  skriv(fid,vout);
  fprintf(fid,'\n ');

% 8 --------------------------------------------
  if ~exist('rcurrt','var')
    rcurrt = 0;
    for j = 1:nh
      for k = 1:nw
	rcurrt = rcurrt+jphi(j,k)*rg(k)*100;
      end
    end
    rcurrt = rcurrt/sum(sum(jphi));
  end
  if ~exist('zcurrt','var')
    zcurrt = 0;
    for j = 1:nh
      for k = 1:nw
	zcurrt = zcurrt+jphi(j,k)*zg(j)*100;
      end
    end
    zcurrt = zcurrt/sum(sum(jphi));
  end
  aspect = rout/aout;
  if ~exist('qsta','var')
    qsta   = abs(rcentr*bcentr/mu0/cpasma*(eout^2+1)/2/aspect^2);
  end
  if ~exist('betat','var'),  betat  = 0; end
  skriv(fid,rcurrt); % R of current centroid? in cm
  skriv(fid,zcurrt); % Z of current centroid? in cm
  skriv(fid,qsta);
  skriv(fid,betat);
  fprintf(fid,'\n ');

% 9 --------------------------------------------
  if ~exist('betap','var'), betap = 0; end
  if ~exist('ali','var')
    if exist('li','var')
      ali = li;
    else
      ali = 0;
    end
  end
  if ~exist('oleft','var') % Closest distance to limiter
    oleft   = 1e10;
  end
  if ~exist('oright','var'),  oright  = 1e10; end
  skriv(fid,betap);
  skriv(fid,ali);
  skriv(fid,oleft);
  skriv(fid,oright);
  fprintf(fid,'\n ');

% 10 --------------------------------------------
  if ~exist('otop','var'), otop = 1e10; end
  if ~exist('obott','var'), obott = 1e10; end
  if ~exist('qpsib','var')
    psibar = linspace(0,1,nw);
    qpsib = spline(psibar(1:nw-1),qpsi(1:nw-1),0.95);
  end
  if ~exist('vertn','var'),  vertn  = 0; end
  skriv(fid,otop);
  skriv(fid,obott);
  skriv(fid,qpsib);
  skriv(fid,vertn);
  fprintf(fid,'\n ');
  
% 11 --------------------------------------------
  if ~exist('rco2v','var'), rco2v = zeros(1,mco2v); end
  if ~exist('dco2v','var'), dco2v = zeros(1,mco2v); end
  for j = 1:mco2v
    skriv(fid,rco2v(j));
  end
  fprintf(fid,'\n ');

% 12 --------------------------------------------
  for j = 1:mco2v
    skriv(fid,dco2v(j));
  end
  fprintf(fid,'\n ');
  
% 13 --------------------------------------------
  if ~exist('rco2r','var'), rco2r = zeros(1,mco2r); end
  if ~exist('dco2r','var'), dco2r = zeros(1,mco2r); end
  for j = 1:mco2r
    skriv(fid,rco2r(j));
  end
  fprintf(fid,'\n ');

% 14 --------------------------------------------
  for j = 1:mco2r
    skriv(fid,dco2r(j));
  end
  fprintf(fid,'\n ');

% 15 --------------------------------------------
  if ~exist('shearb','var'), shearb = 0; end
  if ~exist('bpolav','var'), bpolav = 0; end
  if ~exist('s1','var'),   s1   = 0; end
  if ~exist('s2','var'),  s2  = 0; end
  skriv(fid,shearb);
  skriv(fid,bpolav);
  skriv(fid,s1);
  skriv(fid,s2);
  fprintf(fid,'\n ');

% 16 --------------------------------------------
  if ~exist('s3','var'), s3 = 0; end
  if ~exist('qout','var'), qout = 0; end
  if ~exist('olefs','var'),   olefs   = 0; end
  if ~exist('orighs','var'),  orighs  = 0; end
  skriv(fid,s3);
  skriv(fid,qout);
  skriv(fid,olefs);
  skriv(fid,orighs);
  fprintf(fid,'\n ');

% 17 --------------------------------------------
  if ~exist('otops','var'), otops = 0; end
  if ~exist('sibdry','var'), sibdry = 0; end
  if ~exist('areao','var'), % Area of plasma in cm^2
    areao   = 0;
    for j = 1:nh
      for k = 1:nw
        if jphi(j,k) ~= 0
          areao = areao+dr*dz*1e4;
	end
      end
    end
  end
  if ~exist('wplasm','var'),  wplasm  = 0; end
  skriv(fid,otops);
  skriv(fid,sibdry);
  skriv(fid,areao);
  skriv(fid,wplasm);
  fprintf(fid,'\n ');

% 18 --------------------------------------------
  if ~exist('terror','var'), terror = 0; end
  if ~exist('elongm','var'), elongm = 0; end
  if ~exist('qqmagx','var'), qqmagx = qpsi(1); end
  if ~exist('cdflux','var'), cdflux  = 0; end
  skriv(fid,terror);
  skriv(fid,elongm);
  skriv(fid,qqmagx);
  skriv(fid,cdflux);
  fprintf(fid,'\n ');

% 19 --------------------------------------------
  if ~exist('alpha','var'), alpha = 0; end
  if ~exist('rttt','var'), rttt = 0; end
  if ~exist('psiref','var'),   psiref   = 0; end
  if ~exist('xndnt','var'),  xndnt  = 0; end
  skriv(fid,alpha);
  skriv(fid,rttt);
  skriv(fid,psiref);
  skriv(fid,xndnt);
  fprintf(fid,'\n ');

% 20 --------------------------------------------
  if ~exist('rseps1','var'), rseps1 = 0; end
  if ~exist('zseps1','var'), zseps1 = 0; end
  if ~exist('rseps2','var'), rseps2 = 0; end
  if ~exist('zseps2','var'), zseps2 = 0; end
% THIS NEEDS TO BE FIXED to correctly handle empty objects.
  if isempty(rseps1), rseps1=0; end
  if isempty(zseps1), zseps1=0; end
  if isempty(rseps2), rseps2=0; end
  if isempty(zseps2), zseps2=0; end
  skriv(fid,rseps1);
  skriv(fid,zseps1);
  skriv(fid,rseps2);
  skriv(fid,zseps2);
  fprintf(fid,'\n ');

% 21 --------------------------------------------
  if ~exist('sepexp','var'), sepexp = 0; end
  if ~exist('obots','var'), obots = 0; end
  if ~exist('btaxp','var'), btaxp = 0; end
  if ~exist('btaxv','var'), btaxv = 0; end
  skriv(fid,sepexp);
  skriv(fid,obots);
  skriv(fid,btaxp);
  skriv(fid,btaxv);
  fprintf(fid,'\n ');

% 22 --------------------------------------------
  if ~exist('aaq1','var'), aaq1 = 0; end
  if ~exist('aaq2','var'), aaq2 = 0; end
  if ~exist('aaq3','var'), aaq3 = 0; end
  if ~exist('seplim','var'), seplim = 0; end
  skriv(fid,aaq1);
  skriv(fid,aaq2);
  skriv(fid,aaq3);
  skriv(fid,seplim);
  fprintf(fid,'\n ');

% 23 --------------------------------------------
  if ~exist('rmagx','var'), rmagx = 100*rmaxis; end
  if ~exist('zmagx','var'), zmagx = 100*zmaxis; end
  if ~exist('simagx','var'), simagx = psimag/2/pi; end
  if ~exist('taumhd','var'), taumhd = 0; end
  skriv(fid,rmagx);
  skriv(fid,zmagx);
  skriv(fid,simagx);
  skriv(fid,taumhd);
  fprintf(fid,'\n ');

% 24 --------------------------------------------
  if ~exist('betapd','var'), betapd = 0; end
  if ~exist('betatd','var'), betatd = 0; end
  if ~exist('wplasmd','var'), wplasmd = 0; end
  if ~exist('fluxx','var'), fluxx = 0; end
  skriv(fid,betapd);
  skriv(fid,betatd);
  skriv(fid,wplasmd);
  skriv(fid,fluxx);
  fprintf(fid,'\n ');

% 25 --------------------------------------------
  if ~exist('vloopt','var'), vloopt = 0; end
  if ~exist('taudia','var'), taudia = 0; end
  if ~exist('qmerci','var'), qmerci = 0; end
  if ~exist('tavem','var'), tavem = 0; end
  skriv(fid,vloopt);
  skriv(fid,taudia);
  skriv(fid,qmerci);
  skriv(fid,tavem);
  fprintf(fid,'\n ');

% 26 --------------------------------------------
  if ~exist('nsilop','var'), nsilop = 1; end
  if ~exist('magpri','var'), magpri = 1; end
  if ~exist('nfcoil','var'), nfcoil = 1; end
  if ~exist('nesum','var')
    if exist('ecurrt','var')
      nesum = length(ecurrt);
      if nesum == 0
        nesum = 1;
	ecurrt = 0;
      end
    else
      nesum = 1;
    end
  end
  fprintf(fid,'%5i',nsilop);
  fprintf(fid,'%5i',magpri);
  fprintf(fid,'%5i',nfcoil);
  fprintf(fid,'%5i\n ',nesum);

% 27 --------------------------------------------
  if ~exist('csilop','var'), csilop = zeros(nsilop,1); end
  if ~exist('cmpr2','var'), cmpr2 = zeros(magpri,1); end
  if ~exist('ccbrsp','var'), ccbrsp = zeros(nfcoil,1); end
  if ~exist('ecurrt','var'), ecurrt = zeros(nesum,1); end
  for j = 1:nsilop+magpri
    if j<=nsilop
      skriv(fid,csilop(j));
    else
      skriv(fid,cmpr2(j-nsilop));
    end
    if mod(j,4) == 0
      fprintf(fid,'\n ');
    end
  end
  if mod(j,4) ~= 0
    fprintf(fid,'\n ');
  end

% 27+ --------------------------------------------
  for j = 1:nfcoil
    skriv(fid,ccbrsp(j));
    if mod(j,4) == 0
      fprintf(fid,'\n ');
    end
  end
  if mod(j,4) ~= 0
    fprintf(fid,'\n ');
  end

% 28+ --------------------------------------------
  for j = 1:nesum
    skriv(fid,ecurrt(j));
    if mod(j,4) == 0
      fprintf(fid,'\n ');
    end
  end
  if mod(j,4) ~= 0
    fprintf(fid,'\n ');
  end

% 29+ --------------------------------------------
  if ~exist('pbinj','var'), pbinj = 0; end
  if ~exist('rvsin','var'), rvsin = 0; end
  if ~exist('zvsin','var'), zvsin = 0; end
  if ~exist('rvsout','var'), rvsout = 0; end
  skriv(fid,pbinj);
  skriv(fid,rvsin);
  skriv(fid,zvsin);
  skriv(fid,rvsout);
  fprintf(fid,'\n ');
  
  if ~exist('zvsout','var'), zvsout = 0; end
  if ~exist('vsurfa','var'), vsurfa = 0; end
  if ~exist('wpdot','var'), wpdot = 0; end
  if ~exist('wbdot','var'), wbdot = 0; end
  skriv(fid,zvsout);
  skriv(fid,vsurfa);
  skriv(fid,wpdot);
  skriv(fid,wbdot);
  fprintf(fid,'\n ');
  
  if ~exist('slantu','var'), slantu = 0; end
  if ~exist('slantl','var'), slantl = 0; end
  if ~exist('zuperts','var'), zuperts = 0; end
  if ~exist('chipre','var'), chipre = 0; end
  skriv(fid,slantu);
  skriv(fid,slantl);
  skriv(fid,zuperts);
  skriv(fid,chipre);
  fprintf(fid,'\n ');
  
  if ~exist('cjor95','var'), cjor95 = 0; end
  if ~exist('pp95','var'), pp95 = 0; end
  if ~exist('ssep','var'), ssep = 0; end
  if ~exist('yyy2','var'), yyy2 = 0; end
  skriv(fid,cjor95);
  skriv(fid,pp95);
  skriv(fid,ssep);
  skriv(fid,yyy2);
  fprintf(fid,'\n ');
  
  if ~exist('xnnc','var'), xnnc = 0; end
  if ~exist('cprof','var'), cprof = 0; end
  if ~exist('oring','var'), oring = 0; end
  if ~exist('cjor0','var'), cjor0 = 0; end
  skriv(fid,xnnc);
  skriv(fid,cprof);
  skriv(fid,oring);
  skriv(fid,cjor0);
  fprintf(fid,'\n ');
  
  if ~exist('fexpan','var'), fexpan = 0; end
  if ~exist('qqmin','var'), qqmin = min(qpsi); end
  if ~exist('chigamt','var'), chigamt = 0; end
  if ~exist('ssi01','var'), ssi01 = 0; end
  skriv(fid,fexpan);
  skriv(fid,qqmin);
  skriv(fid,chigamt);
  skriv(fid,ssi01);
  fprintf(fid,'\n ');
  
  if ~exist('fexpvs','var'), fexpvs = 0; end
  if ~exist('sepnose','var'), sepnose = 0; end
  if ~exist('ssi95','var'), ssi95 = 0; end
  if ~exist('rqqmin','var'), rqqmin = 0; end
  skriv(fid,fexpvs);
  skriv(fid,sepnose);
  skriv(fid,ssi95);
  skriv(fid,rqqmin);
  fprintf(fid,'\n ');
  
  if ~exist('cjor99','var'), cjor99 = 0; end
  if ~exist('cjlave','var'), cjlave = 0; end
  if ~exist('rmidin','var') % min(rbbbs) at zbbbs=0
    k = find(zbbbs(1:nbbbs-2).*zbbbs(2:nbbbs-1)<=0);
    [dum, kk] = sort(rbbbs(k)); k = k(kk(1));
    rmidin = interp1(zbbbs(k+[0 1]),rbbbs(k+[0 1]),0);
  end
  if ~exist('rmidout','var')
    k = find(zbbbs(1:nbbbs-2).*zbbbs(2:nbbbs-1)<=0);
    [dum, kk] = sort(rbbbs(k)); k = k(kk(end));
    rmidout = interp1(zbbbs(k+[0 1]),rbbbs(k+[0 1]),0);
  end
  skriv(fid,cjor99);
  skriv(fid,cjlave);
  skriv(fid,rmidin);
  skriv(fid,rmidout);
  fprintf(fid,'\n ');
  
  if ~exist('psurfa','var'), psurfa = 0; end
  if ~exist('peak','var'), peak = 0; end
  if ~exist('dminux','var'), dminux = 0; end
  if ~exist('dminlx','var'), dminlx = 0; end
  skriv(fid,psurfa);
  skriv(fid,peak);
  skriv(fid,dminux);
  skriv(fid,dminlx);
  fprintf(fid,'\n ');
  
  if ~exist('dolubaf','var'), dolubaf = 0; end
  if ~exist('dolubafm','var'), dolubafm = 0; end
  if ~exist('diludom','var'), diludom = 0; end
  if ~exist('diludomm','var'), diludomm = 0; end
  skriv(fid,dolubaf);
  skriv(fid,dolubafm);
  skriv(fid,diludom);
  skriv(fid,diludomm);
  fprintf(fid,'\n ');
  
  if ~exist('ratsol','var'), ratsol = 0; end
  if ~exist('rvsiu','var'), rvsiu = 0; end
  if ~exist('zvsiu','var'), zvsiu = 0; end
  if ~exist('rvsid','var'), rvsid = 0; end
  skriv(fid,ratsol);
  skriv(fid,rvsiu);
  skriv(fid,zvsiu);
  skriv(fid,rvsid);
  fprintf(fid,'\n ');
  
  if ~exist('zvsid','var'), zvsid = 0; end
  if ~exist('rvsou','var'), rvsou = 0; end
  if ~exist('zvsou','var'), zvsou = 0; end
  if ~exist('rvsod','var'), rvsod = 0; end
  skriv(fid,zvsid);
  skriv(fid,rvsou);
  skriv(fid,zvsou);
  skriv(fid,rvsod);
  fprintf(fid,'\n ');
  
  if ~exist('zvsod','var'), zvsod = 0; end
  if ~exist('condno','var'), condno = 0; end
  if ~exist('psin32','var'), psin32 = 0; end
  if ~exist('psin21','var'), psin21 = 0; end
  skriv(fid,zvsod);
  skriv(fid,condno);
  skriv(fid,psin32);
  skriv(fid,psin21);
  fprintf(fid,'\n ');

  if ~exist('rq32in','var'), rq32in = 0; end
  if ~exist('rq21top','var'), rq21top = 0; end
  if ~exist('chilibt','var'), chilibt = 0; end
  if ~exist('ali3','var'), ali3 = 0; end
  skriv(fid,rq32in);
  skriv(fid,rq21top);
  skriv(fid,chilibt);
  skriv(fid,ali3);
  fprintf(fid,'\n ');

  fprintf(fid,'\n ');

  if ~exist('header','var'), header = blanks(42); end
  if ~exist('fit_type','var'), fit_type = blanks(3); end
  fprintf(fid,' %42s',header);
  fprintf(fid,' %3s',fit_type);
    
  fclose(fid);
  return
  
  % function to print to file EXACTLY like Fortran's write(fid,1040) does
  function skriv(fid,val)
  if val >= 0
    fprintf(fid,' %11.9E',val);
  else
    fprintf(fid,'%11.9E',val);
  end
  return
  
  % function copied from efit code
  function dismin = dslant(x,y,np,xmin,xmax,ymin,ymax,x1,y1,x2,y2)
  dismin=inf;
  delx=x2-x1;
  dely=y2-y1;
  dels=sqrt(delx^2+dely^2);
  nn=dels/0.002;
  nn=max(5,nn);
  delx=delx/(nn-1);
  dely=dely/(nn-1);
  for j=1:nn
    xw=x1+delx *(j-1);
    yw=y1+dely *(j-1);
    for m=1:np
      if x(m) > xmin & x(m)<xmax & y(m) > ymin & y(m) < ymax
        disw=sqrt((xw-x(m))^2+(yw-y(m))^2);
        dismin=min(dismin,disw);
      end
    end
  end
  dismin=dismin*100;
  return

