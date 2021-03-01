function gfile_data= ...
         read_gfile_tok(filename,tokamak,nfcoil,nesum,nves,old_efit,cc_file)
 %
%  USAGE:
%   gfile_data = ...
%        read_gfile_tok(filename,tokamak,nfcoil,nesum,nves,old_efit,cc_file)
%
%   gfile_data = read_gfile_tok(filename); % default DIII-D read
%
%  PURPOSE: Read efit ASCII g-eqdsk file for specific tokamak.
%
%  INPUTS: <default>
%
%     filename= g0-eqdsk filename (printed with iecurr=2 or iplcout=1)
%               Note: filename= gShot.Time name is parced to get shot number 
%     tokamak = one of 'DIII-D','d3d','NSTX','EAST','KSTAR','ITER','CTF'
%               PEGASUS, (Not CASE sensitive so DIII-D <=> diii-d)
%               if not one of above then need nfcoil,nesum,...
%               <'DIII-D'>)
%     nfcoil =  number of F-coils <default specific to each machine - below>
%     nesum =   number of E-coils <default specific to each machine - below>
%     nves =    number of vessel elements included in efit data <defaults below>
%     old_efit= For iecurr=2 several different flavors of efit output
%               exists. (Optional, default=0)  For old versions of efit
%               (e.g. NSTX) use old_efit=1 since it uses older efit version.
%     cc_file=  coil current file name [cc=load(cc_file)] or coil current vector
%               use when coil currents not available from efit (or to override)
%		(units = MA-turns)
%
%     Current Machine Defaults: (Default tokamak is DIII-D)
%	Tokamak	nfcoil	nesum	nves   old_efit Size(cc)	Note
%	DIII-D	18	6	0	0	18+6=24		6-seg E-coil
%	NSTX	11	1	30	1	41+1=42		1E+11F+30VV
%               after shot 115265, April 2 2005 set old_efit=0 or shot#>=115265
%	NSTX	17	1	35	1	1+17+35=42	1E+17F+35VV
%               NSTX internally switches on shot # determined from gfile name
%                    ex: /g120423.00533 => shot= 120423 => old_efit=0
%	KSTAR	14	1	0	0	14+4=18		no E; 4IC's
%	 sab10	14	1	60+24	0	14+4=18		12 VV Groups
%        rtefit 18      1       60+10	0	18 ic included	33x33 14VCefit
%	EAST	12	1	0	0	14		no E; IC's?
%	CTF	14	1	24	0	14		NO E;
%	FDF	14	0	24	0	14		NO E; FDF2009?
%               22      0       24      0       22              FDF2011
%       PEGASUS 17      1       0       0       17
%       ITER    12      1       113     0       12              NO E 
%       HL2M    16      1       24      0       16              1E 
%	other	nfcoil	nesum	0	0	nfcoil+nesum	need coil #'s	
%
%      Note: Including nfcoil, nesum ... in argument list with a particular
%            machine overrides the machine default numbers
%
%      Note: filename= "..../gShot.Time" format is parced to get shot number
%            if needed to differentiate between machine configurations (i.e. NSTX)
%
%  OUTPUTS:
%   gfile_data = structure containing data read directly from g-file and 
%		 variables derived from this original data.
%
%  RESTRICTIONS:
%     g0 file must be produced by EFIT run with flag iecurr=2 or iplcout=1
%         in order to write current density distribution on grid and coil cc's
 
%  METHOD:  
%  Follows write/read format of EFIT file weqdskx.for, reads usual g-file data,
%      reads coil current and jphi data (if iecurr was set =2 in EFIT run),
%      defines some flux variables and converts some variables to D3D Matlab
%      data environment standards.

%  WRITTEN BY:  Jim Leuer ON	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)read_gfile_tok.m	1.23 11/02/11

 if nargin==0
    help read_gfile_tok
    return
 end

 if nargin==1
    tokamak='DIII-D'
 end
 if nargin <= 2
    nfcoil = [];
 end
 if nargin <= 3
    nesum = [];
 end
 if nargin <= 4
    nves = [];
 end
 if nargin <= 5
    old_efit= [];
 end
 if nargin <= 6
    cc_file= [];
 end

 jphi=[]; % needed to make sure this is a variable not a function
 gfile_data  = [];	% return empty data if error in processing

% ===========================================================================
  switch upper(tokamak)

% ===========================================================================
  case {'DIII-D','D3D','DIIID'}
% ===========================================================================

    if isempty(old_efit)
       old_efit= 0;
    end
    if isempty(nfcoil)   nfcoil=  18; end
    if isempty(nesum)    nesum=    6; end

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nesum = nesum;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);


% ===========================================================================
  case {'EAST'}
% ===========================================================================

    if isempty(old_efit)
       old_efit= 0;
    end
   if isempty(nfcoil)   nfcoil=  13; end
    if isempty(nesum)    nesum=    0; end

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nesum = nesum;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);


% ===========================================================================
  case {'KSTAR'}
% ===========================================================================

    if isempty(old_efit)
       old_efit= 0;
    end
    if isempty(nfcoil)   nfcoil=  14; end
    if isempty(nesum)    nesum=    1; end

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end

% Modify read based on gfile header ecase date which is from EFITD fortran
    ecase= gdata.ecase;
    if ~isempty(findstr('01/23/2002',ecase)) % this is SABBAGH EFIT
      dum= sscanf(ecase(27:38),'%f');
      if ~isempty(dum)
      	shot= dum(1);
         if length(dum)>=2
           time= dum(2);
         end
         if shot<=3000 % ? guess where is shot transition from 2009 to 2010?
%           Need to reread since has F and VV as active coils
            nfcoil= 14;
            nves=    56;
            nesum=    1;
           [gdata,ireadok] = read_gfile_func(filename,nfcoil+nves,nesum,old_efit);
           if(isempty(gdata))
              wait(['ERROR read_gfile_tok: unable to read ' filename])
           return;
    end
         else
%           Need to reread since has F and VV as active coils
            nfcoil= 14;
            nves=    60;
            nesum=    1;
           [gdata,ireadok] = read_gfile_func(filename,nfcoil+nves,nesum,old_efit);
           if(isempty(gdata))
              wait(['ERROR read_gfile_tok: unable to read ' filename])
           return;
    end
         end % if shot
       end % if ~isempty(dum)
    end

    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nesum = nesum;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);


% ===========================================================================
  case {'NSTX'}
% ===========================================================================
%   CHECK FOR OLD OR NEW NSTX CONFIGURATION: Two EFIT builds are supported:
%             data     config_name            shot       old_efit
%        old: ~2002    02072002Av1.0         <~115100	  1
%        new: apr05,   04202005Av1.0/        >115265      0
% Note: some shots (100's) below 115265 sabaugh was switching different efits
%       and nether version will work. these were done during Apr 05 startup
    if isempty(old_efit)
       [shot,tms,tree,server] = shot_from_gfile(filename); % shot,tms,tree/dir,ser
 	 if ~isempty(shot)
	    if shot < 100000  shot > 1000000
	       disp(['CAUTION: read_gfile_tok Doesnt understand shot# from: ', ...
	            filename])
	       wait(' ? CONTINUE using NEW Ver: apr05 04202005Av1.0/ shot>115265') 
	       old_efit=0;
	    elseif shot < 115265 % Transition to New EFIT at this shot (gates)
	       old_efit=1;
	    else
	       old_efit=0;
	    end
         else
	    disp(['CAUTION: read_gfile_tok No shot Number from filename', ...
	            filename])
	    wait(' ? CONTINUE using NEW Ver: apr05 04202005Av1.0/ shot>115265') 
            old_efit= 0; % this corresponds to New: apr05 efit
         end
    end
         
  if old_efit==0 % NEW Apr 05 EFIT Config: 04202005Av1.0
    
    if isempty(nfcoil),   nfcoil=17; end  % default NSTX values
    if isempty(nves),     nves=35; end    % default NSTX values
    if isempty(nesum),    nesum=1; end    % default NSTX values

    [gdata,ireadok] = read_gfile_func(filename,nfcoil+nves,nesum,0);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    gfile_data = gdata;

    gfile_data.config_name= '04202005Av1.0';
    if isfield(gdata,'nh') &  isfield(gdata,'nw') 
        gfile_data.config_name= [gfile_data.config_name ...
	                         '_' int2str(gdata.nw) int2str(gdata.nh)];
    end
  
    disp(['read_gfile_tok Reading New NSTX: config_name= ' ...
          gfile_data.config_name])

    id= find(filename=='/');
    if ~isempty(id) & filename(id(end)+1)=='k'
      disp('%!!! read_gfile_tok NSTX  CORRECTING WRONG JSOLVER pcurrt & jphi !!!')
      disp('%in read_gfile_tok NSTX NOT TESTED !!!')
      gfile_data.pcurrt= -gdata.pcurrt*dr*dz;
      gfile_data.jphi= -gdata.jphi*dr*dz;
      disp(['% Plasma Current from jphi= ',num2str(sum(jphi(:))*dr*dz)])
      disp(['% Plasma Current cpasma= ',num2str(gdata.cpasma)])
    end

  else % if old_efit==1 OLD Feb 2002 EFIT BELOW  Config_name= '02072002Av1.0';

    if isempty(nfcoil),   nfcoil=11; end    % default NSTX values
    if isempty(nves),     nves=30; end    % default NSTX values
    if isempty(nesum),    nesum=1; end       % default NSTX values

    [gdata,ireadok] = read_gfile_func(filename,nfcoil+nves,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    gfile_data = gdata;

    gfile_data.config_name= '02072002Av1.0';
    disp(['read_gfile_tok Reading Old NSTX: config_name= ' gfile_data.config_name])

    id= find(filename=='/');
    if ~isempty(id) & filename(id(end)+1)=='k'
      disp('%!!! read_gfile_tok NSTX  CORRECTING WRONG JSOLVER pcurrt & jphi !!!')
      disp('%in read_gfile_tok NSTX NOT TESTED !!!')
      gfile_data.pcurrt= -gdata.pcurrt*dr*dz;
      gfile_data.jphi= -gdata.jphi*dr*dz;
      disp(['% Plasma Current from jphi= ',num2str(sum(jphi(:))*dr*dz)])
      disp(['% Plasma Current cpasma= ',num2str(gdata.cpasma)])
    end

  end % if old_efit 

  options.old_efit = old_efit;
  options.nfcoil = nfcoil;
  options.nves = nves;
  options.nesum = nesum;
  gfile_data = std_efit_units(gfile_data,upper(tokamak),options);
  
% ===========================================================================
  case {'NSTXU'}
% ===========================================================================
%  20121225: MJL Modified to specify fcturn(:)=1
%             data     config_name            shot       old_efit
%      only: 12nov,   12212012v0.0        >115265      0
%
    
    if isempty(nfcoil),   nfcoil=14; end  % default NSTXU values
    if isempty(nves),     nves=40; end    % default NSTXU values
    if isempty(nesum),    nesum=1; end    % default NSTXU values
    if isempty(old_efit), old_efit= 0; end

    [gdata,ireadok] = read_gfile_func(filename,nfcoil+nves,nesum,0);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    gfile_data = gdata;
    
    if isfield(gdata,'nh') &  isfield(gdata,'nw') 
        gfile_data.config_name= [int2str(gdata.nw) int2str(gdata.nh)];
    end

    id= find(filename=='/');
    if ~isempty(id) & filename(id(end)+1)=='k'
      disp('%!!! read_gfile_tok NSTX  CORRECTING WRONG JSOLVER pcurrt & jphi !!!')
      disp('%in read_gfile_tok NSTX NOT TESTED !!!')
      gfile_data.pcurrt= -gdata.pcurrt*dr*dz;
      gfile_data.jphi= -gdata.jphi*dr*dz;
      disp(['% Plasma Current from jphi= ',num2str(sum(jphi(:))*dr*dz)])
      disp(['% Plasma Current cpasma= ',num2str(gdata.cpasma)])
    end

  options.old_efit = old_efit;
  options.nfcoil = nfcoil;
  options.nves = nves;
  options.nesum = nesum;
  gfile_data = std_efit_units(gfile_data,upper(tokamak),options);
  
%Workaround for cc vector until we can get cc's in geqdsk from LRDFIT
  if(~isfield(gfile_data,'cc')), gfile_data.cc=zeros(1+nfcoil,1); end

% ===========================================================================
  case {'ITER'}
% ===========================================================================

    if isempty(old_efit)
       old_efit= 0;
    end
    nfcoil = 12;
    nesum = 0;

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nesum = nesum;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);

% ===========================================================================
  case {'CTF'}
% ===========================================================================

    if isempty(old_efit)
       old_efit= 0;
    end
    if exist('filename')~=1 | isempty(filename)
     disp('% NOTE: No filename GIVEN: Defaulting to gfile below:')
     filename='/u/leuer/efit/ctf/run/g010010.00008'
    end    % default
    if isempty(nfcoil),   nfcoil=  14; end
    if isempty(nesum),    nesum=    1; end
    if isempty(old_efit), old_efit=0; end % read new

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    gfile_data = gdata;
    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nesum = nesum;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);

% ===========================================================================
  case {'FDF'}
% ===========================================================================

    if isempty(old_efit)
       old_efit= 0;
    end
    if exist('filename')~=1 | isempty(filename)
     disp('% NOTE: No filename GIVEN: Defaulting to gfile below:')
     filename='/m/GAtools/tokamaks/fdf/define/g020000.00215_FDF2_Cur'
    end    % default

% read first line to see whcih efit is there
    fid= fopen(filename, 'r');
    if fid == -1
       wait([' %ERROR read_gfile_tok FDF filename :',filename])
    end
    ecase= fgetl(fid);

    if ~isempty(findstr('08/14/2006',ecase))  % original 2009 EFIT
%      filename='/m/GAtools/tokamaks/fdf/define/g020000.00215_FDF2_Cur'
       icase = 0;
       nfcoil=  14
       nesum=    1;
       old_efit=0;
       turnfc0= 1;
    end

    if ~isempty(findstr('03/28/2011',ecase))  % new FDF 2011 version
%      filename= '/usc-data/m/GAtools/tokamaks/fdf/efit/g000138.00650'; 
       icase = 1;
       nfcoil=  22
       nesum=    1;
       old_efit=0;
       turnfc0= 50; % 50turns (multiplyier of experimental current )
    end

    if isempty(nfcoil),   nfcoil=  14; end
    if isempty(nesum),    nesum=    1; end
    if isempty(old_efit), old_efit=0; end % read new

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    gfile_data = gdata;
    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nves = nves;
    options.nesum = nesum;
    options.turnfc0 = turnfc0;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);

% ===========================================================================
  case {'CFETR'}
% ===========================================================================

    if isempty(old_efit)
       old_efit= 0;
    end
    nfcoil = 14;
    nesum = 0;

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nesum = nesum;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);

% ===========================================================================
  case {'PEGASUS'}
% ===========================================================================

    if isempty(old_efit)
       old_efit= 0;
    end
    if isempty(nfcoil)   nfcoil=  17; end
    if isempty(nesum)    nesum=    1; end

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    gfile_data = gdata;
    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nesum = nesum;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);

% ===========================================================================
  case {'SST'}
% ===========================================================================

    if isempty(old_efit)
       old_efit= 0;
    end
    if isempty(nfcoil)   nfcoil=  14; end
    if isempty(nesum)    nesum=    1; end

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    gfile_data = gdata;
    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nesum = nesum;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);

% ===========================================================================
  case {'HL2M'}
% ===========================================================================
% jal 2014jul2 EFIT from LIJX at Swip
%   disp('111111')
    if isempty(old_efit)
       old_efit= 0;
    end
    if isempty(nfcoil)   nfcoil=  16; end
    if isempty(nesum)    nesum=    1; end

    [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
    if(isempty(gdata))
       wait(['ERROR read_gfile_tok: unable to read ' filename])
       return;
    end
    options.old_efit = old_efit;
    options.nfcoil = nfcoil;
    options.nesum = nesum;
    gfile_data = std_efit_units(gdata,upper(tokamak),options);

% check currents in Fcoils= -3300, -3000, ... looks good ecurrt=0

% ===========================================================================
  otherwise
% ===========================================================================

  disp(['%ERROR: read_gfile_tok does not recognize tokamak:' tokamak])

   [gdata,ireadok] = read_gfile_func(filename,nfcoil,nesum,old_efit);
   if(isempty(gdata))
      wait(['ERROR read_gfile_tok: unable to read ' filename])
      return;
   end
   gfile_data = gdata;

end  % switch


% ===========================================================================
% Overwrite coil data if cc_file not empty
% ===========================================================================
  if ~isempty(cc_file)
     if isnumeric(cc_file)
        cc = cc_file;
        gfile_data.cc = cc;
        gfile_data.cc2 = cc;
     elseif exist(cc_file)==2
        cc = load(cc_file);
       gfile_data.cc = cc;
       gfile_data.cc2 = cc;
     else
        disp(['%ERROR: read_gfile_tok could not set cc based on cc_file'])
     end
  end

  gfile_data.gdef.filename = 'read_gfile_tok input "filename"';

  return

% ========================
% TESTING
% ========================
% kstar sab09 test jal01jun2010
  addpath /home/leuer/matlab/efit/
  addpath /home/leuer/tokamaks/kstar/efit/sab/
  clear
  filename ='/home/leuer/tokamaks/kstar/efit/sab/gfiles/g001938.00109'
      nfcoil= 14;
      nves=    56;
      nesum=    1;
 [gdata,ireadok] = read_gfile_func(filename,nfcoil+nves,nesum,0)

 filename ='/home/leuer/tokamaks/kstar/efit/sab/gfiles/g001938.00109'
 gdata= read_gfile_tok(filename,'kstar')

 filename ='/home/leuer/tokamaks/kstar/efit/sab/gfiles/g001938.00109'
 efit_gfile= filename;
 build_kstar_sys

% Old nstx test:
  clear
  figure(1), clf
  filename ='/home/leuer/tokamaks/nstx/efit/gfiles/g110843.00294'
   tokamak='nstx'
%  gfile_data = read_gfile_tok(filename,tokamak)
  gfile_data = read_gfile_tok(filename,tokamak)
  struct_to_ws(gfile_data);
  clear gfile_data
  save old
  figure(1), clf
  [c,h]= contour(rg,zg,pcurrt);
  clabel(c,h);
   hold on
   plot(rbbbs,zbbbs)
   title(filename)
  axis image

% New April 05 NSTX EFIT
  clear
  figure(2), clf
%  filename ='/home/leuer/tokamaks/nstx/efit/gfiles/g120434.00617_iecurr2'
%  filename ='/home/leuer/tokamaks/nstx/efit/gfiles/g120423.00553' % neg Ip
%  dir_name='/home/leuer/tokamaks/nstx/efit/gfiles/';
%  dir(dir_name)
%  filename ='/home/leuer/tokamaks/nstx/efit/gfiles/g122645.00781'; % Good Ref
  filename ='/m/GAtools/tokamaks/nstx/efit/gfiles/g122645.00781'; % Good Ref
  tokamak='nstx'
%  gfile_data = read_gfile_tok(filename,tokamak)
  gfile_data = read_gfile_tok(filename,tokamak)
  struct_to_ws(gfile_data);
  clear gfile_data
  save new
  figure(2), clf
  [c,h]= contour(rg,zg,pcurrt);
  clabel(c,h);
   hold on
   plot(rbbbs,zbbbs)
   title(filename)
  axis image

  diff_mat('old.mat','new.mat')


% kstar test:
  filename= '/u/leuer/efit/kstar/run/g010002.01010';
  gfile_data = read_gfile_tok(filename,'kstar'); 
  gfile_data = read_gfile_tok(filename,'kstar',[],[],[],[],[]); 
  gfile_data.pcurrt
  
%  gfile_data = read_gfile_tok([],'CTF'); % default read
 
%  figure(1)
%  clf 
%  [c,h]= contour(gfile_data.rg,gfile_data.zg,gfile_data.jphi,10,'r');
%  hold on
%  axis equal
%  plot(gfile_data.rbbbs(1:gfile_data.nbbbs),...
%       gfile_data.zbbbs(1:gfile_data.nbbbs),'k')
  

% 1) brsp is Amps (or Amp-turns), depending on FCTURN in mhdin.dat. For 
% Kstar FCTURN=14*1 so brsp is in amp-turns just like DIII-D. If FCTURN 
% has actual number of turns in coil then objects are in AMP's
% 
% 2) cc is derived from brsp. Daves original definition of cc is 
% MA-turns and the 1e+6 Multiplier always needed to convert brsp to 
% cc[MA-turns]. FCTURN may also be needed in conversion if EFIT was 
% generated with FCTURN~=1 so we always generate cc in Amp-turn. This 
% has been the standard in all read-gfile stuff.
% 
% 3) Here is description EFIT mhdin.dat definitions of variables as I 
% believe them
% 
% FCTURN= turns multiplier in generating EFIT Greens Functions.
% TURNFC= turns multiplier of input current for fitting
% FCID=   Coil grouping vector.
%          (ex: [1 2 2 3] => 3 groups with 2nd group containing seried coil
% 
% 4) Here are variables for EAST MHDIN.dat
% 
%   FCID=1., 2., 3., 4.,  4., 5., 6.,
%        7., 8., 9., 10., 10. 11.,  12.
%   FCTURN=1,1,1,0.193548,0.806452,1,1,1,1,1,0.193548,0.806452,1,1
%   TURNFC= 3*140.0, 248.0,  60.0, 32.0,
%           3*140.0, 248.0,  60.0, 32.0
% Note that the FCTURN & FCID make a 1 turn system
% 
% 5) Here are variables for  KSTAR MHDIN.dat
%   FCID=1., 2., 3., 4., 5., 6.,
%        7., 8., 9., 10., 11., 12., 13., 14.
%   FCTURN=14*1., ECTURN=1*1.
%   TURNFC= 180, 144, 72, 108, 224, 128, 72,
%           180, 144, 72, 108, 224, 128, 72
% 
% 
% 6) mhdin.dat has ALL information we need and should be used for all 
% this input. It can even get most dimensions since fcid, fcturn and 
% turnfc read in by namelist will have dimensions associated with them. 
% At present it can be hardwired into read_gfile. Later we should use 
% my read_mhdin.m to do all work so it is transparent. This information 
% and the proj_turn function should allow us to do any conversion we 
% desire.

% 
%
%  MODIFICATION HISTORY:
%   5/15/06 jal implemented cc_file per RD and new NSTX
%   Robert Deranian ON 4/3/06 - Added PEGASUS case
%   jal 4/11/06 fix empty brsp ecurrt
% pegasis mod number of turns 6/8/06 jal/walker
