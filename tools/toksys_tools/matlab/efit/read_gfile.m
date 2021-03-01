 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  USAGE:  >>read_gfile
%
%  PURPOSE: Script to read efit ASCII g-eqdsk file w/ flag iecurr=2 set
%
%  INPUTS: [default]
%     filename = ['g083750.04838']; % g0-eqdsk filename - Can contain directory
%     nfcoil = [18]; % # of F-coils; must be specified for non-D3D (Optional)
%			 defaults to D3D values of nfcoil=18
%     nesum =  [6];  % # of E-coils; must be specified for non-D3D (Optional)
%			 defaults to D3D value of nesum=6 (for > 6seg case)
%     old_efit=[0];  % For iecurr=2 several different flavors of efit output
%                        exists. (Optional) For old versions of efit (i.e. NSTX)
%                        you want old_efit=1 since they use older efit version.
%
%  OUTPUTS:
%   (Produces large set of data objects from g0-file in Matlab environment.)
%    ireadok = flag to report good read of gfile (0=bad, 1=good)
%    jphi=current density on grid  MA/m^2
%    psizr = true total flux on grid in Wb
%    psimag = axis flux in true Wb
%    psibnd = boundary flux in true Wb  (psibry also defined same)
%    cc = E/F coil currents in MA-turns 
%    cc2 = E/F coil currents in MA-turns (for 2-segment E-coil)
%          CAUTION: cc2 has min(2,nesum) E-coil elements (could be 1 vs std: 2)
%
%
%  RESTRICTIONS:
%     g0 file must be produced by EFIT run with flag iecurr=2 in order to
%      write current density distribution on grid.
%   Must set nesum to correct number of e-coil segments at top of file.
%   May overwrite some data objects in D3D environment from load_d3denv.
%
%  METHOD:  
%  Follows write/read format of EFIT file weqdskx.for, reads usual g-file data,
%      reads coil current and jphi data (if iecurr was set =2 in EFIT run),
%      defines some flux variables and converts some variables to D3D Matlab
%      data environment standards.

%  WRITTEN BY:  Jim Leuer ON	
%
%  MODIFICATION HISTORY:
%     DAH  Added new variable definitions and unit converted quantities.
%     DAH  8/8/97  Modified to read new (~>5/2/97) g0-file format *only*
%     DAH  1/13/00 Modified to have nfcoil, nesum defaulted like read_gfile1
%		Now read_gfile *should* be usable for all non-D3D devices 
%		as well...
%     JAL 11-15-02 Fixed so it reads old efit iecurr=2 format for old_efit=1
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prelims:
 if exist('nfcoil')~=1, nfcoil=18; end  %default to D3D values...
 if exist('nesum')~=1, nesum=6; end
 if exist('old_efit')~=1, old_efit=0; end % change in structure below iecurr=2

 %Former harwired values:
   %nfcoil = 18;
   %%nesum = 2;   %for pre-Ecoil breakage 2-segment E-coil.
   %nesum = 6;    %for post-Ecoil breakage 6-segment E-coil.

% -------------------------- Open File: filename

   if exist('filename')~=1
      filename= 'g083750.04838'
   end
   
   if ~isstr(filename)
      disp([' %ERROR read_gfile: Filename must be a string',filename]);
   end

  fid= fopen(filename, 'r');

   if fid == (-1)
      disp([' %ERROR read_gfile: Couldn''t open file',filename]);
      ireadok= 0;
      return
   end

   ireadok= 1;

%Derive shotname from file:
   shotname = filename(end-11:end-6); %str w/shotname,assume shot.time format

% -------------------------- Start Read of File

  disp('First line:')
%      read (neqdsk,2000,err=30000) (case(i),i=1,6),imfit,nw,nh
  sline= fgetl(fid);
  disp(sline)
  ecase= sline(1,1:48);
  [dum,count]= sscanf(sline(1,49:length(sline)),'%d');
  imfit= dum(1);
  nw=    dum(2);
  nh=    dum(3);

%     read (neqdsk,2020) xdim,zdim,rzero,rgrid(1),zmid
  dum= fscanf(fid,'%f',5);
  xdim=     dum(1);
  zdim=     dum(2);
  rzero=    dum(3);
  rgrid1= dum(4);
  zmid=     dum(5);

%      read (neqdsk,2020) rmaxis,zmaxis,ssimag,ssibry,bzero
  dum= fscanf(fid,'%f',5);
  rmaxis= dum(1);
  zmaxis= dum(2);
  ssimag= dum(3);
  ssibry= dum(4);
  bzero=  dum(5);

%      write (neqdsk,2020) cpasma(jtime),ssimag,xdum,rmaxis,xdum
  dum= fscanf(fid,'%f',5);
  cpasma= dum(1);
  ssimag= dum(2);
  xdum=   dum(3);
  rmaxis= dum(4);
  xdum=   dum(5);

%      write (neqdsk,2020) zmaxis,xdum,ssibry,xdum,xdum
  dum= fscanf(fid,'%f',5);
  zmaxis= dum(1);
  xdum=   dum(2);
  ssibry= dum(3);
  xdum=   dum(4);
  xdum=   dum(5);

%      write (neqdsk,2020) (fpol(i),i=1,nw)
  fpol= fscanf(fid,'%f',nw);

%      write (neqdsk,2020) (pres(i),i=1,nw)
  pres= fscanf(fid,'%f',nw);

% ffprim and pprime are input with sign correct for EFIT
% eqdsk is written out as a negative if plasma current is positive

%      write (neqdsk,2020) (workk(i),i=1,nw)= -ffprim(i)
  ffprim= -sign(cpasma)*fscanf(fid,'%f',nw);

%      write (neqdsk,2020) (workk(i),i=1,nw)= -pprime
  pprime= -sign(cpasma)*fscanf(fid,'%f',nw);

%      write (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
  psirz= fscanf(fid,'%f',[nw,nh]);
%  note: nw= plasma width is stored in columns; nh= height in rows

%      write (neqdsk,2020) (qpsi(i),i=1,nw)
  qpsi= fscanf(fid,'%f',nw);

%      write (neqdsk,2022) nbbbs,limitr
  dum= fscanf(fid,'%f',2);
  nbbbs= dum(1);
  limitr= dum(2);

%      write (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
  dum= fscanf(fid,'%f',[2,nbbbs]);
  rbbbs= dum(1,:)';
  zbbbs= dum(2,:)';
  
%      write (neqdsk,2020) (xlim(i),ylim(i),i=1,limitr)
  dum= fscanf(fid,'%f',[2,limitr]);
  xlim= dum(1,:)';
  ylim= dum(2,:)';


% ------------------------end of original eqdsk read

%    fclose(fid);
%    clear dum sline xdum; 
%    return

%-------------------
% Create useful fluxes (even for iecurr not = 2)
% Convert flux objects to REAL units:
  psizr = -psirz'*2*pi;
  psimag = -ssimag*2*pi;
  psibry = -ssibry*2*pi;
  psibnd = -ssibry*2*pi;


% ------------------ 
% Continued read for read_eqdsk_iec2.m below:

%      read next line and display to see where we are
%  sline= fgetl(fid);
%  sline
%  sline= fgetl(fid);
%  sline

%      read (neqdsk,2024,err=30900,end=30900) kvtor,rvtor,nmass
  kvtor= 0;
  dum= fscanf(fid,'%f',3);
  if (isempty(dum))    fclose(fid); return, end
    kvtor= dum(1);
    rvtor= dum(2);
    nmass= dum(3);

    if (kvtor > 0)
%      read (neqdsk,2020) (presw(i),i=1,mw)
       presw= fscanf(fid,'%f',nw);
%      read (neqdsk,2020) (preswp(i),i=1,mw)
       preswp= fscanf(fid,'%f',nw);
    end

% Read Ion Mass Density and associated Data (if present):
   if nmass>0
%    read (neqdsk,2020) (dmion(i),i=1,nw)
     dmion = fscanf(fid,'%f',nw);
   end

%% Loop to find ishot (if exists):
%   while 1
%     dum = fscanf(fid,'%f',1);
%     ***WORK HERE***
%   end
 
% Read New (sometime > 5/2/97) data in g-file:
  if old_efit==0
    rhovn = fscanf(fid,'%f',nw);    %always written now in g0-files!!
    keecur = fscanf(fid,'%f',1);
    if keecur>0
      workk=fscanf(fid,'%f',nw);
    end        
  end
   
% Read next line of (run-on) integers (and display to see where we are):
%  disp('Line after kvtor, etc...:')
    dum = fscanf(fid,'%f',4);
    disp('ISHOT data line to check:')
    disp(dum)
  if (isempty(dum))    
    fclose(fid); 
    disp('No IECURR=2 data.')
    return
  end
  disp('IECURR=2 data being read...')

%  read (neqdsk) rgrid(1),rgrid(nw),zgrid(1),zgrid(nh)
     dum = fscanf(fid,'%f',4);
     rgefit=linspace(dum(1),dum(2),nw)';
     zgefit=linspace(dum(3),dum(4),nh)';
     drefit = (dum(2)-dum(1))/(nw-1);
     dzefit = (dum(4)-dum(3))/(nh-1);
%     disp('Grid extreme R,Z values are :')
%     disp(dum)

% read (neqdsk) (brsp(i),i=1,nfcoil) Lang says this is Amp-Turns (fitted)
    brsp = fscanf(fid,'%f',nfcoil);

% read (neqdsk) (ecurrt(i),i=1,nesum)
    ecurrt = fscanf(fid,'%f',nesum);

% read (neqdsk) (pcurrt(i),i=1,nwnh)
    [pcurrt,count] = fscanf(fid,'%f',nw*nh);   %current in grid elements (A)
    pcurrt = reshape(pcurrt,nh,nw);

% Convert pcurrt to nh x nw array and scale to MA/m^2:
   %jphi=zeros(nh,nw);
   jphi = pcurrt*1.e-6/(dzefit*drefit);

% Construct coil current vector:
   dum= min(2,nesum); % if nesum=1 then use ONLY 1 ECOIL current
   cc2 = [ecurrt(1:dum);brsp]*1.e-6; %cc-vector with 2 e-coil segments (MA-t)
   cc = [ecurrt;brsp]*1.e-6;   %cc-vector with all e-coil segs(2 or 6) (MA)

  clear dum
  
  fclose(fid);

%  return
    

        
