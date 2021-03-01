  function  eq= mds_eq(shot, tree, gam, mk_var, server)
%
%  SYNTAX:
%         eq= mds_eq(shot, tree, gam, mk_var);    % defaults to DIII-D server
%             mds_eq(shot, 'EFIT01, 'gam', 1, 'NSTX'); % makes vars in workspace 
%
%  PURPOSE: Get EFIT equilibrium (GEQDSK AEQDSK MEQDSK) from mdsplus database. 
%
%  INPUT: <default>
%    shot   = shot number
%    tree   = which efit tree to use <'EFIT01'>
%    gam    = ['g']; get only g file, a= only a file, m= mfile, <gam>
%    mk_var = [1]; Will make efit variables in calling workspace, <0> makes none
%                      2 will make lower case, 3 will make upper case variable
%                      note default variables in MDS are upper case
%                  NOTE: see eq_to_env for way to make tstart<=t<tend arrays
%    server = which data base to use: defaults to [DIII-D] also NSTX
%
%  OUTPUT:
%    eq      = structure containing EFIT eqdsk variables
%              See: eq.gnames .anames .mnames for g,a,m file names
%                   All EFIT variables are as in MDSPLUS => Upper Case
%              Ex:  eq.BCENTR, eq.BDRY, eq.PSIRZ, 
%              Use: plot(eq.GTIME, eq.CPASMA),
%                  contour(eq.R, eq.Z, eq.PSIRZ(:,:,round(length(eq.GTIME)/2))')
%
%    Extra items added to structure (all lower case)
%      eq.id   = sting array of important data identifyer enf
%      eq.gnames= sting array of all gfile names collected (eq to see all)
%      eq.anames= sting array of all afile names collected (eq to see all)
%      eq.mnames= sting array of all mfile names collected (eq to see all)
%      eq.shot = shot number
%      eq.tree = mds tree used to get data (ex. EFIT01 EFIT02 (=MSD));
%
%    if mk_var>0 it makes variables in base workspace
%                 BCENTR, BDRY, CPASMA, PSIRZ, ...
% 
%  NOTE: 1)Function returns all TIME data for shot in large arrays.
%          Example:
%          eq.BCENTER(255,1)
%          eq.GTIME(255,1); % Time data vector for all .GNAME arrays 255 times
%          eq.ATIME(251,1); % Time data vector for all .ANAME arrays 251 times
%          eq.MTIME or eq.TIME % Note old "m" uses TIME new "m" uses MTIME
%          eq.PPRIME(65,255),PSIRZ(65,65,225) 
%
%          Similar arrays (without eq. strucutre) are generated in workspace 
%          if mk_var>=1 : BCENTER(255,1), GTIME(255,1)  or LC if mk_var==2
%
%        2) Time is in vector eq.GTIME and is in ms. All other units are as
%           they come from EFIT. psi(Vs/r) I(A) R(m) ...
%
% CAUTION: All time vectors are last index of vector of array. 
% CAUTION: GTIME ATIME & (MTIME or TIME) are not necesaarly the same. 
% CAUTION: Current Bad Variable in d3dMDSPLUS ARE: CASE, HEADER, ZGRID (use Z)
% CAUTION: THIS ROUTINE GETS ALL MDS EQDSK DATA FOR SHOT AND EQ CAN BE HUGE
% CAUTION: Afile read 1st then Gfile: variables in common: BCENTR, CPASMA, ...
%
%  SEE ALSO: eq_time_lim eq_ga_env
%

%  WRITTEN BY:  Jim Leuer    ON      6/11/03
%  USES:   eq_mod
% To see MDS structure on HYDRA run traverser
% should work for JET data but not tested
% ==========================================================================
% 
% Defaults
  wait('%CAUTION: mds_eq is OBSOLETE use eq_mds instead')
  return;

  if nargin==4
     server='DIII-D';
  elseif nargin==3
     server='DIII-D';
     mk_var= 0;
  elseif nargin==2
     server='DIII-D';
     mk_var= 0;
     gam= 'gam';
  elseif nargin==1
     server='DIII-D';
     mk_var= 0;
     gam= 'gam';
     tree= 'EFIT01';
  elseif nargin==0
     disp('% mds_eq needs at least a "shot" argument')
     help mds_eq
     return
  end

  if isempty(server)  server= 'DIII-D'; end
  if isempty(tree)    tree= 'EFIT01'; end
  if isempty(gam)     gam= 'gam';     end
  if isempty(mk_var)  mk_var= 0;      end

% -----------------------------------------------
% Open and check conneciton to MDSPLUS data base:
% -----------------------------------------------
  tic
  ier= 0;
% [shoto,status]=mdsopen('atlas.gat.com::EFIT01',shot)
  if strcmp(upper(server),'DIII-D')
    mdsconnect('atlas.gat.com');
  elseif strcmp(upper(server),'NSTX')
    mdsconnect('europa.pppl.gov:8501');
  else
    disp(['%ERROR: Couldnt recognize parameter server =',server])
  end
  [shoto,status]=mdsopen(tree,shot);
  if ~mod(status,2)
    ier=1;
    fprintf(['%ERROR mds_eq: unable to open ' tree ' for shot '  ...
        int2str(shot) '\n'])
    eq=[];
    status=mdsdisconnect;
    return;
  end

% add identifier string to structure:
  id= str2mat('mds_eq ', int2str(shot), tree, date);
  eq.id= id;
  eq.shot= shot;
  eq.tree= tree;
  eq.creator= 'mds_eq';

% ===============================================================
% START READ OF AEQDSK DATA
% ===============================================================

 if ~isempty(strfind(gam,'a')) 

  eq.anamesmds= mdsvalue('getnci(".RESULTS.AEQDSK:*","FULLPATH")');

  if isempty(eq.anamesmds)
    disp('%ERROR mds_eq cant read AEQDSK Names, returning')
    status=mdsdisconnect;
    return
  end

  fullnames= char(eq.anamesmds);
  id= strfind(fullnames(1,:),'.RESULTS.AEQDSK:');
  eq.anames=  fullnames(:,id+16:end);

% loop over all variables
  for ii= 1:length(eq.anames);
%    ii=0; ii=ii+1;
     fullnam= char(eq.anamesmds(ii));
     nam= deblank(eq.anames(ii,:));
     dat=  mdsvalue(fullnam); % 
     str= ['eq.' nam '= dat;'];
     eval(str);
     if mk_var
       switch mk_var
         case 1
           assignin('base',nam,dat);
         case 2
           assignin('base',lower(nam),dat);
         case 3
           assignin('base',upper(nam),dat);
	 end
     end
     if toc > 10
       disp(['%mds_eq still reading at ' int2str(ii) ' of ',...
               int2str(length(eq.anames)),' at variable: ', nam]);
       tic
     end
   end
 end

% ===============================================================
% START READ OF MEQDSK DATA which is in measurments area time in eq.TIME
% ===============================================================

 if ~isempty(strfind(gam,'m')) 

  eq.mnamesmds= mdsvalue('getnci(".MEASUREMENTS:*","FULLPATH")');

  if ~isempty(eq.mnamesmds)
   fullnames= char(eq.mnamesmds);
   id= strfind(fullnames(1,:),'.MEASUREMENTS:');
   eq.mnames=  fullnames(:,id+14:end);

% loop over all variables
   for ii= 1:length(eq.mnames);
%    ii=0; ii=ii+1;
     fullnam= char(eq.mnamesmds(ii));
     nam= deblank(eq.mnames(ii,:));
     dat=  mdsvalue(fullnam); % 
     str= ['eq.' nam '= dat;'];
     eval(str);
     if mk_var
       switch mk_var
         case 1
           assignin('base',nam,dat);
         case 2
           assignin('base',lower(nam),dat);
         case 3
           assignin('base',upper(nam),dat);
	 end
     end
     if toc > 10
       disp(['%mds_eq still reading at ' int2str(ii) ' of ',...
               int2str(length(eq.mnames)),' at variable: ', nam]);
       tic
     end
   end
  else
   disp('%mds_eq COULD NOT FIND ANY EFIT Measurement data: mnamesmds=[]');
  end
  
 end

% ===============================================================
% START READ OF GEQDSK DATA
% ===============================================================

 if ~isempty(strfind(gam,'g'))
 
  eq.gnamesmds= mdsvalue('getnci(".RESULTS.GEQDSK:*","FULLPATH")');

  if ~isempty(eq.gnamesmds)
   fullnames= char(eq.gnamesmds);
   id= strfind(fullnames(1,:),'.RESULTS.GEQDSK:');
   eq.gnames=  fullnames(:,id+16:end);

% loop over all variables
   for ii= 1:length(eq.gnames);
%    ii=0; ii=ii+1;
     fullnam= char(eq.gnamesmds(ii));
     nam= deblank(eq.gnames(ii,:));
     dat=  mdsvalue(fullnam); % 
     if isfield(eq,nam);
        rmfield(eq,nam);
	disp(['%Note: mds_eq is overwriting structure field: ',nam]);
     end
     str= ['eq.' nam '= dat;'];
     eval(str);
     if mk_var
       switch mk_var
         case 1
           assignin('base',nam,dat);
         case 2
           assignin('base',lower(nam),dat);
         case 3
           assignin('base',upper(nam),dat);
	 end
     end
     if toc > 10
       disp(['%mds_eq still reading at ' int2str(ii) ' of ',...
               int2str(length(eq.gnames)),' at variable: ', nam]);
       tic
     end
   end
  else
   disp('%mds_eq COULD NOT FIND ANY EFIT GEQDSK Data: gnamesmds=[]');
  end

 end

% ========================================================
     
  status=mdsdisconnect;

  eq= eq_mod(eq); % fix time so end of vectors
  
  return
  
% Testing SEE test_mds_eq
% (WATCH OUT 114504 has problems use 98549 Ferron High Performance)

%  eqg= mds_eq(114504, 'EFIT01', 'g', 0);
%  eqa= mds_eq(114504, 'EFIT01', 'a', 0);
%  eqm=  mds_eq(98549, 'EFIT01', 'm', 0);

  shot=98549; tree='EFIT01'; gam='gam'; mk_var= 0;
  filename='/u/leuer/efit/diiid/s98549/g098549.04000' %
  read_gfile
  eq=  mds_eq(shot, tree, 'gam', 0);
  eq = mds_eq(shot);
  eq=  mds_eq(shot, tree, 'gm', 0 );
  eq= eq_time_lim(eq,4,4); % make only 4s time
  eq_ga_env(eq,[],[],'g');

%  mds_eq(114504, 'EFIT01',[], 1);
%  mds_eq(114504, 'EFIT02',[], 2);

%   eqg= mds_eq(114504, 'EFIT01', 'g', 0);
%   id= round(length(eqg.GTIME)/2);
%   contour(eqg.R, eqg.Z, eqg.PSIRZ(:,:,id)'); hold on, axis equal
%   plot(eqg.RBBBS(1:eqg.NBBBS(id),id),eqg.ZBBBS(1:eqg.NBBBS(id),id),'k')

% ========================================================
% Check NSTX:
  shot=110843; tree='EFIT01'; gam='gam'; mk_var= 0; server='NSTX'
  clear eq
  eq=  mds_eq(shot, tree, 'gam', 0, server);
  plot(eq.ATIME,eq.IPMHD)

