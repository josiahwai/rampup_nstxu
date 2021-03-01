
% read_afile.m: reads efit ASCII a file with form: A_shot#.time#
%  produced from read_a.for in [.efit] area
%  efit must be run with eqdsk output file set to ASCII (normal mode is binary)
%  READ OF A______._____ FILE PRODUCED BY EFIT
%  READS FILE: a085352.03000 by default
%  WARNING: # OF E-COIL SEGMENTS HARDWIRED BELOW!!!!!!
%
% ---------------------------------------------
%  Input: (Internally default to D3D values)
%       filename= 'a085352.03000 ' = a file name
%       nsilop=     41               = number of psi loops
%       magpri=     60               = number of mag probes
%       nfcoil=     18               = number of fcoils
%       nesum=       6               = number of e-coil circuits(A,B, for d3d)
%
%  Output:
%        (all variables in a file with names based on EFIT convention)
%
%       aaq1           eccurt         pasmat         taudia         
%       aaq2           filename       pbinj          taumhd         
%       aaq3           elongm         psiref         tavem          
%       ali            eout           qmflag         terror         
%       alpha          fluxx          qout           time           
%       aout           header         qpsib          tsaisq         
%       areao          ishot          qqmagx         uday           
%       bcentr         jflag          qsta           vertn          
%       betap          ktime          rcencm         vloopt         
%       betapd         lflag          rco2r          vout           
%       betat          limloc         rco2v          vsurfa         
%       betatd         magpri         rcurrt         wbdot          
%       bpolav         mco2r          rmagx          wpdot          
%       btaxp          mco2v          rout           wplasm         
%       btaxv          mfver1         rseps1         wplasmd        
%       ccbrsp         mfver2         rseps2         xndnt          
%       cdflux         nesum          rttt           xxx            
%       chipre         nfcoil         rvsin          xxxxx          
%       cmerci         nsilop         rvsout         xxxxxx         
%       cmpr2          obots          s1             zcurrt         
%       cpasma         obott          s2             zmagx          
%       csilop         olefs          s3             zout           
%       dco2r          oleft          sepexp         zseps1         
%       dco2v          orighs         seplim         zseps2         
%       diamag         oright         shearb         zvsin          
%       doutl          otop           sibdry         zvsout         
%       doutu          otops          simagx         

% ---------------------------------------------
% Jim Leuer, General Atomics, 7-27-95
% ---------------------------------------------
% Modifications:
%   9/12/97  Added read of nsilop, magpri, nfcoil, nesum
%  12/04/00  Removed clear of ireadok from end so ireadok passes to
%		calling routine for read assessment  DAH
%

% DEFAULTS: (D3D)

%   clear ireadok
   format compact;

% NO LONGER NEEDED:
%   if (exist('nsilop') ~= 1), nsilop= 41, end
%   if (exist('magpri') ~= 1), magpri= 60, end
%   if (exist('nfcoil') ~= 1), nfcoil= 18, end
%   if (exist('nesum')  ~= 1), nesum =  6, end
   
% -------------------------- Open File: filename
   if ~isstr(filename)
      disp([' %ERROR read_afile: filename must be a string: ',filename]);
   end

   if (exist(filename) ~= 2)
      disp([' %ERROR read_afile: Can''t find file: ',filename]);
   end

  fid= fopen(filename, 'r');

   if fid == (-1)
      disp([' %ERROR read_afile: Couldn''t open file: ',filename]);
      ireadok= 0;
      return
   end

   ireadok= 1;

%Derive shotname from file:
   shotname = filename(end-11:end-6); %str w/shotname,assume shot.time format

% 1 -------------------------- Start Read of File

% READ (neqdsk,1055) uday,mfver1,MFVER2; 1055 format (1x,a10,2a5)

   line=   fgetl(fid);
    uday=   line(1,1:11);
    mfver1= line(1,12:16);
    mfver2= line(1,17:min(length(line),21));

% 2 --------------------------------------------
% ? watch out not implemented:C        if (ishot.le.99999) then
% write (neqdsk,1053) ishot,ktime; 1053 format (1x,i6,11x,i5)

   dum= fscanf(fid,'%f',2);
   ishot=     dum(1);
   ktime=     dum(2);

   line=   fgetl(fid); % needed after fscanf to read end of line charac.

% 3 --------------------------------------------
% READ (neqdsk,1060) time;
   line=   fgetl(fid);
    [dum,count]= sscanf(line,'%f');
    time= dum(1);

% 4 --------------------------------------------
% READ (neqdsk,1060) time,jflag,lflag,limloc,mco2v,mco2r,qmflag
% 1060 format (1h ,f7.2,10x,i5,11x,i5,1x,a3,1x,i3,1x,i3,1x,a3)
   line=   fgetl(fid);
    limloc= line(41:43);
    qmflag= line(53:55);
    [dum,count]= sscanf(line(2:40),'%f');
     time1= dum(1);
     jflag= dum(2);
     lflag= dum(3);
    [dum,count]= sscanf(line(44:52),'%f');
     mco2v= dum(1);
     mco2r= dum(2);

% 5 --------------------------------------------
% READ (neqdsk,1040) tsaisq ,rcencm,bcentr ,pasmat 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    tsaisq= dum(1);
    rcencm= dum(2);
    bcentr= dum(3);
    pasmat= dum(4);

% 6 --------------------------------------------
% READ (neqdsk,1040) cpasma ,rout ,zout ,aout 
 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    cpasma= dum(1);
    rout= dum(2);
    zout= dum(3);
    aout= dum(4);

% 7 --------------------------------------------
% READ (neqdsk,1040) eout ,doutu ,doutl ,vout 
  
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    eout= dum(1);
    doutu= dum(2);
    doutl= dum(3);
    vout= dum(4);

% 8 --------------------------------------------
% READ (neqdsk,1040) rcurrt ,zcurrt ,qsta ,betat 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    rcurrt= dum(1);
    zcurrt= dum(2);
    qsta= dum(3);
    betat= dum(4);

% 9 --------------------------------------------
% READ (neqdsk,1040) betap ,ali ,oleft ,oright 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    betap= dum(1);
    ali= dum(2);
    oleft= dum(3);
    oright= dum(4);

% 10 --------------------------------------------
% READ (neqdsk,1040) otop ,obott ,qpsib ,vertn 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    otop= dum(1);
    obott= dum(2);
    qpsib= dum(3);
    vertn= dum(4);

% 11,12,13,14 --------------------------------------------
%   
   if mco2r > 4 | mco2v > 4 
     disp('%ERROR READ_AFILE: not set up to read more than 4 co2 items')
     return
   end

   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    for i=1:count
       rco2v(i)= dum(i);
    end

   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    for i=1:count
       dco2v(i)= dum(i);
    end

   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    for i=1:count
       rco2r(i)= dum(i);
    end

   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    for i=1:count
       dco2r(i)= dum(i);
    end

% 15 --------------------------------------------
% READ (neqdsk,1040) shearb ,bpolav ,s1 ,s2 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    shearb= dum(1);
    bpolav= dum(2);
    s1= dum(3);
    s2= dum(4);

% 16 --------------------------------------------
% READ (neqdsk,1040) s3 ,qout ,olefs ,orighs 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    s3= dum(1);
    qout= dum(2);
    olefs= dum(3);
    orighs= dum(4);

% 17 --------------------------------------------
% READ (neqdsk,1040) otops ,sibdry ,areao ,wplasm 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    otops= dum(1);
    sibdry= dum(2);
    areao= dum(3);
    wplasm= dum(4);

% 18 --------------------------------------------
% READ (neqdsk,1040) terror ,elongm ,qqmagx ,cdflux 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    terror= dum(1);
    elongm= dum(2);
    qqmagx= dum(3);
    cdflux= dum(4);

% 19 --------------------------------------------
% READ (neqdsk,1040) alpha ,rttt ,psiref ,xndnt 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    alpha= dum(1);
    rttt= dum(2);
    psiref= dum(3);
    xndnt= dum(4);

% 20 --------------------------------------------
% READ (neqdsk,1040) rseps1,zseps1,rseps2,zseps2
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    rseps1= dum(1);
    zseps1= dum(2);
    rseps2= dum(3);
    zseps2= dum(4);

% 21 --------------------------------------------
% READ (neqdsk,1040) sepexp ,obots ,btaxp ,btaxv 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    sepexp= dum(1);
    obots= dum(2);
    btaxp= dum(3);
    btaxv= dum(4);

% 22 --------------------------------------------
% READ (neqdsk,1040) aaq1 ,aaq2 ,aaq3 ,seplim 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    aaq1= dum(1);
    aaq2= dum(2);
    aaq3= dum(3);
    seplim= dum(4);

% 23 --------------------------------------------
% READ (neqdsk,1040) rmagx ,zmagx ,simagx ,taumhd 
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    rmagx= dum(1);
    zmagx= dum(2);
    simagx= dum(3);
    taumhd= dum(4);

% 24 --------------------------------------------
% READ (neqdsk,1040) betapd ,betatd ,wplasmd ,fluxx
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    betapd= dum(1);
    betatd= dum(2);
    wplasmd= dum(3);
    fluxx= dum(4);
    diamag=fluxx /1.0e-03;

% 25 --------------------------------------------
% READ (neqdsk,1040) vloopt ,taudia ,cmerci ,tavem
   line=   fgetl(fid);
   [dum,count]= sscanf(line,'%f');
    vloopt= dum(1);
    taudia= dum(2);
    cmerci= dum(3);
    tavem= dum(4);

% 25a --------------------------------------------
%  READ (neqdsk,****) nsilop, magpri, nfcoil, nesum
    nsilop = fscanf(fid,'%i',1);
    magpri = fscanf(fid,'%i',1);
    nfcoil = fscanf(fid,'%i',1);
    nesum = fscanf(fid,'%i',1);

% 26 --------------------------------------------
%  READ (neqdsk,1040) (csilop(k),k=1,nsilop),(cmpr2(k),k=1,magpri)
    csilop= fscanf(fid,'%f',nsilop);
    cmpr2=  fscanf(fid,'%f',magpri);

% 27 --------------------------------------------
%  READ (neqdsk,1040) (ccbrsp(k),k=1,nfcoil)
    ccbrsp= fscanf(fid,'%f',nfcoil);

% 28 --------------------------------------------
%  READ (neqdsk,1040) (eccurt(k),k=1,nesum)
    eccurt= fscanf(fid,'%f',nesum);

% 29 --------------------------------------------
%  READ (neqdsk,1040) pbinj ,rvsin ,zvsin ,rvsout 
    dum= fscanf(fid,'%f',4);
     pbinj= dum(1);
     rvsin= dum(2);
     zvsin= dum(3);
     rvsout= dum(4);

% 30 --------------------------------------------
%  READ (neqdsk,1040) zvsout ,vsurfa ,wpdot ,wbdot 
    dum= fscanf(fid,'%f',4);
     zvsout= dum(1);
     vsurfa= dum(2);
     wpdot= dum(3);
     wbdot= dum(4);

% 31 --------------------------------------------
%  READ (neqdsk,1040) slantu,slantl,zuperts,chipre
    dum= fscanf(fid,'%f',4);
     slantu= dum(1);
     slantl= dum(2);
     zuperts= dum(3);
     chipre= dum(4);

% 32 --------------------------------------------
%  READ (neqdsk,1040) cjor95,pp95,ssep,yyy2
    dum= fscanf(fid,'%f',4);
     cjor95= dum(1);
     pp95= dum(2);
     ssep= dum(3);
     yyy2= dum(4);

% 33 --------------------------------------------
%  READ (neqdsk,1040) xnnc,cprof,oring,cjor0
    dum= fscanf(fid,'%f',4);
     xnnc = dum(1);
     cprof = dum(2);
     oring = dum(3);
     cjor0 = dum(4);

% 34 --------------------------------------------
%  READ (neqdsk,1040) xdum,xdum1,xdum2.xdum3
    dum= fscanf(fid,'%f',4);
    xxxxx = dum(1);
    xxx = dum(2);
    xxx = dum(3);
    xxx = dum(4);

% 35 --------------------------------------------
%   READ (neqdsk,1040) xxxxxx,xxx,xxx,xxx
    dum= fscanf(fid,'%f',4);
    xxxxxx = dum(1);
    xxx = dum(2);
    xxx = dum(3);
    xxx = dum(4);

    line=   fgetl(fid); % need to read end of line after fscanf & before fgetl

% 36 --------------------------------------------
%   READ (neqdsk,1042) header; 1042 format (1x,a42)
    line=   fgetl(fid);
    header= line(2:length(line));
% ------------------------end of a_file read

    fclose(fid);
    
    clear ans count fid i dum line xdum time1; 

    disp(['  %READ_AFILE just read efit A_file: ',filename])

%    return

% xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% junk below


