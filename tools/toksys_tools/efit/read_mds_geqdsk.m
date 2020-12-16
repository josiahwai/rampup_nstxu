function [gdata,neq,eq]= read_mds_geqdsk(shotnum,times,efit_source,tokamak)
 %
%  SYNTAX: [gdata,neq,eq] = read_mds_geqdsk(shotnum,times,efit_source,tokamak)
%
%  PURPOSE:  Read G=file data from mdsplus, return data like read_gfile_tok.m.
%*************************************************************************
%     OBSOLETE FUNCTION - please use read_mds_eqdsk
%*************************************************************************
%  INPUT:
%	shotnum = shot number
%	times   = equilibrium time(s), one of:
%			if single float#: closest time available in mds (ms)
%			if pair [t1,t2]: times between t1 and t2 (ms)
%			if 'all': return all times in tree
%	efit_source = source tree for EFIT (e.g. 'EFIT01', 'EFITRT1')
%	tokamak = name of tokamak
%
%  OUTPUT:
%	gdata = array of structures, each containing info from 1 geqdsk
%       (single time causes structure returned to be like read_gfile_tok; 
%		others return extended structures)
%	neq   = number of eqdsk's returned
%	eq    = entire tree containing efit equilibrium data
% 
%  RESTRICTIONS: Only handles D3D, EAST, and NSTX(EFITRT only) right now.  Need
%	to add code for the other devices supported by read_gfile_tok.m.

%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	10/8/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)read_mds_geqdsk.m	1.6 01/22/10

wait('read_mds_geqdsk is OBSOLETE - please replace with read_mds_eqdsk');
return;

 if nargin<3
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

gdata=[]; neq=[]; eq=[];

if(isstr(times))
   if(~strcmp(upper(times),'ALL'))
      wait('ERROR read_mds_geqdsk: invalid value for times')
      return;
   end
end

  jphi=[]; % needed to make sure this is a variable not a function

switch upper(tokamak)

% ===========================================================================
  case {'DIII-D','D3D','DIIID'}
% ===========================================================================

    nfcoil=  18;
    nesum=    6;

%    eq = eq_mds(shotnum,efit_source,'DIII-D',-1);
%    geq = eq.results.geqdsk;
    [gdata,neq,eq] = read_mds_g_func(shotnum,efit_source,'DIII-D');
%    meq = eq.measurements;

% include d3d build parameters (see build area => must be same)
    gdata.fcturn = ones(1,18);
    gdata.turnfc = [58 58 58 58 58 55 55 58 55 58 58 58 58 58 55 55 58 55];
    gdata.fcid = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
%   NOTE: EFIT IS CONSTRUCTED with ecid having 1,2,3,4,5,6 numbers. Even are 1's
%         odds hook to 2
    gdata.ecid = [ ...                     % must be same as ecdata(5,:)
                        1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 ...
                        1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 ...
                        1 2 2 1 1 2 2 1 1 1 1 1 2 2 2 2 2 2 2 1 ...
                        1 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 ...
                        2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 ...
                        2 2 1 1 2 2 1 1 2 2 2 2 2 1 1 1 1 1 1 1 ...
                        2 2];
    gdata.ecturn = ones(122,1);    % required if ecid is not empty


   nee = length(unique(gdata.ecid));
   for k=1:neq
    if ~isempty(gdata.data(k).brsp) & ~isempty(gdata.data(k).ecurrt)

% Construct E-coil Number of turns to make output MA-turns:
    dum= min(2,nesum); % if nesum=1 then use ONLY 1 ECOIL current

%  find reduced ecoil set data 
      clear neturn
      for ii=1:nee
         idx = find(gdata.ecid==ii);
         neturn(ii,1) = sum(gdata.ecturn(idx));
      end

 %cc-vector with 2 e-coil segments (MA-t):
    cc2 = [gdata.data(k).ecurrt(1:dum).*neturn(1:dum);gdata.data(k).brsp]*1.e-6;
%    cc2 = [diag(neturn(1:dum))*meq.cecurr(1:dum,:);meq.ccbrsp]*1.e-6;

% CAUTION: we are not using the 6segment e-coil anymore so this is crippled
%cc-vector with all e-coil segs(2 or 6) (MA-t):
%    cc = [gdata.ecurrt;gdata.brsp]*1.e-6;
    cc= cc2;

    gdata.data(k).cc     = cc;
    gdata.data(k).cc2    = cc2;
 
    gdata.gdef.cc = 'E/F coil currents in MA-turns ';
    gdata.gdef.cc2 = 'E/F coil currents in MA-turns (for 2-segment E-coil). CAUTION: cc2 has min(2,nesum) E-coil elements (could be 1 vs std: 2)';
    end
   end

% ===========================================================================
  case {'east','EAST'}
% ===========================================================================

    nfcoil=  12;
    nesum=    1;

    [gdata,neq,eq] = read_mds_g_func(shotnum,efit_source,'thor');

    ncadd=4;    % fix cc so it has extra control coils
    for k=1:length(gdata.data)
       if ~isempty(gdata.data(k).brsp)
          cc = [gdata.data(k).brsp*1e-6;zeros(ncadd,1)];
          cc2= cc;

          gdata.data(k).cc  = cc;
          gdata.data(k).cc2 = cc2;
       end
    end
    gdata.gdef.cc = 'PF coil currents in MA-turns ';
    gdata.gdef.cc2 = 'PF coil currents in MA-turns ';

% include east build parameters 

    gdata.fcturn = [1 1 1 44/(44+204) 204/(44+204) 1 1 ...
                         1 1 1 44/(44+204) 204/(44+204) 1 1]';
    gdata.turnfc = [140 140 140 44 204 60 32 140 140 140 44 204 60 32]';
    gdata.fcid = [1 2 3 4 4 5 6 7 8 9 10 10 11 12];
    gdata.ecid = [];

    gdata.data(k).cc     = cc;
    gdata.data(k).cc2    = cc2;
 
    gdata.gdef.cc = 'E/F coil currents in MA-turns ';
    gdata.gdef.cc2 = 'E/F coil currents in MA-turns (for 2-segment E-coil). CAUTION: cc2 has min(2,nesum) E-coil elements (could be 1 vs std: 2)';

% ===========================================================================
  case {'nstx','NSTX'}
% ===========================================================================

%   CHECK FOR OLD OR NEW NSTX CONFIGURATION: Two EFIT builds are supported:
%             data     config_name            shot       old_efit
%        old: ~2002    02072002Av1.0         <~115100	  1
%        new: apr05,   04202005Av1.0/        >115265      0
% Note: some shots (100's) below 115265 sabaugh was switching different efits
%       and nether version will work. these were done during Apr 05 startup

    if shotnum < 100000
       disp('ERROR: read_gfile_tok doesnt know shots less than 100k')
       old_efit=-1; % this is error output
       gdata=[];
       return
    elseif shotnum < 115265 % Transition to New EFIT at this shot (gates)
       old_efit=1;
    else
       old_efit=0;
    end

% -------------------------------------------------------
    if old_efit==0 % NEW Apr 05 EFIT Config: 04202005Av1.0

      nfcoil=17;  % default NSTX values
      nves=35;    % default NSTX values
      nesum=1;    % default NSTX values

      [gdata,neq,eq]=read_mds_g_func(shotnum,efit_source,'nstx');

      gdata.time = gdata.time * 1e3;

      gdata.config_name= '04202005Av1.0';
      disp(['read_gfile_tok Reading New NSTX: config_name= ' ...
          gdata.config_name])

% Coil Note except for changes to pf1aul,
%      below I add a 2nd vector which adds from the Old to the New(Apr05)
      gdata.fcid= [[1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 11] ...
                        [12 13 14 15 16 17]];
      gdata.turnfc=[[20 14 14 15 15 9 8 12 12 12 12 9 8 15 15 14 14 20 32]...
                         [1 1 1 1 48 48]];

      gdata.fcturn = ...
     [[1 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1 1]...
      [1 1 1 1 1 1]];

% vv has some differences internal so new from "make"
      gdata.vvid = [ ...
        1 2 3 4 4 5 6 6 7 7 7 8 8 8 8 8 8 8 8 9 9 9 9 9 10 10 11 11 11 12 ...
        13 14 15 16 16 16 17 17 18 18 18 18 18 19 19 19 19 19 20 20 21 ...
        22 22 23 24 25 26 26 27 27 28 29 30 31 32 33 34 35];

      gdata.vvfrac = [ ...
      1.0000    1.0000    1.0000    0.5000    0.5000    1.0000    0.5000 ...
      0.5000    0.3333333 0.3333333 0.3333333 0.1250    0.1250    0.1250 ...
      0.1250    0.1250    0.1250    0.1250    0.1250    0.2000    0.2000 ...
      0.2000    0.2000    0.2000    0.5000    0.5000    0.3333333 0.3333333 ...
      0.3333333 1.0000    1.0000    1.0000    1.0000    0.3333333 0.3333333 ...
      0.3333333 0.5000    0.5000    0.2000    0.2000    0.2000    0.2000 ...
      0.2000    0.2000    0.2000    0.2000    0.2000    0.2000    0.5000 ...
      0.5000    1.0000    0.5000    0.5000    1.0000    1.0000    1.0000 ...
      0.5000    0.5000    0.5000    0.5000    1.0000    1.0000    1.0000 ...
      1.0000    1.0000    1.0000    1.0000    1.0000];

% ecoil - no change old to new
      gdata.ecid = ones(240,1);
      gdata.ecturn = [ ...
       4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
       4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
       4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 ...
     4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ];

% Changing cc and cc2 so they are amp-turns
% brsp is terminal current (A) NOT Amp-turns because NSTX uses groups and turns
% for each group. 

      for k=1:neq
        if ~isempty(gdata.data(k).brsp) & ~isempty(gdata.data(k).ecurrt)
          gdata.data(k).vc = gdata.data(k).brsp(nfcoil+1:nfcoil+nves)*1.e-6;
          cc_terminal=[gdata.data(k).ecurrt; gdata.data(k).brsp(1:nfcoil)]*1.e-6;
 
          turns = zeros(2,1);	%force to be column
          for j=1:nesum
             turns(j) = sum(gdata.ecturn);
          end
          for j=1:nfcoil
             idx = find(gdata.fcid==j);
             turns(j+nesum) = sum(gdata.turnfc(idx));
          end

          gdata.data(k).cc     = cc_terminal.*turns;
          gdata.data(k).cc2    = gdata.data(k).cc;
        end
      end
      gdata.gdef.cc =  'Fcoil currents in MA-turns ';
      gdata.gdef.cc2 = 'Fcoil currents in MA-turns ';
      gdata.gdef.vc =  'vessel currents in MA';
  
% ----------------------------------------------------------------------------
    else % if old_efit==1 OLD Feb 2002 EFIT BELOW  Config_name= '02072002Av1.0';
% ----------------------------------------------------------------------------
      nfcoil=11; 	% default NSTX values
      nves=30; 		% default NSTX values
      nesum=1; 		% default NSTX values
  
      [gdata,neq,eq]=read_mds_g_func(shotnum,efit_source,'nstx');

      gdata.time = gdata.time * 1e3;

      gdata.config_name= '02072002Av1.0';
      disp(['read_gfile_tok Reading Old NSTX: config_name= ' gdata.config_name])

      gdata.fcid = [1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 11];
      gdata.turnfc = [48 14 14 15 15 9 8 12 12 12 12 9 8 15 15 14 14 48 32];
      gdata.fcturn = ...
       [1 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1 1];

      gdata.vvid = [ ...
	1 2 3 4 4 5 6 6 6 6 6 6 6 7 7 7 7 8 8 8 9 9 9 10 10 11 11 12 12 12 ...
	13 13 13 14 14 14 14 15 15 15 15 15 15 15 16 17 17 18 19 20 21 21 ...
	22 22 23 24 25 26 27 28 29 30];
      gdata.vvfrac = [ ...
   	1 1 1 0.5 0.5 1 ...
	0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 ...
	0.25 0.25 0.25 0.25 ...
	0.33333 0.33333 0.33333 0.33333 0.33333 0.33333 ...
	0.5 0.5 0.5 0.5 ...
	0.33333 0.33333 0.33333 0.33333 0.33333 0.33333 ...
	0.25 0.25 0.25 0.25 ...
	0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 ...
	1 0.5 0.5 1 1 1 0.5 0.5 0.5 0.5 1 1 1 1 1 1 1 1];

      gdata.ecid = ones(240,1);
      gdata.ecturn = [ ...
     4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 ...
     4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 4.0875 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     4.033333 4.033333 4.033333 4.033333 4.033333 4.033333 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.995833 3.995833 3.995833 3.995833 3.995833 3.995833 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ...
     3.933333 3.933333 3.933333 3.933333 3.933333 3.933333 ];

% Changing cc and cc2 so they are amp-turns
% brsp is terminal current (A) NOT Amp-turns because NSTX uses groups and turns
% for each group. 

      for k=1:neq
        if ~isempty(gdata.data(k).brsp) & ~isempty(gdata.data(k).ecurrt)
          gdata.data(k).vc = gdata.data(k).brsp(nfcoil+1:nfcoil+nves)*1.e-6;
          cc_terminal=[gdata.data(k).ecurrt; gdata.data(k).brsp(1:nfcoil)]*1.e-6;

          turns = zeros(2,1);	%force to be column
          for j=1:nesum
             turns(j) = sum(gdata.ecturn);
          end
          for j=1:nfcoil
             idx = find(gdata.fcid==j);
             turns(j+nesum) = sum(gdata.turnfc(idx));
          end

          gdata.cc     = cc_terminal.*turns;
          gdata.cc2    = gdata.cc;

        end
      end
      gdata.gdef.cc =  'Fcoil currents in MA-turns ';
      gdata.gdef.cc2 = 'Fcoil currents in MA-turns ';
      gdata.gdef.vc =  'vessel currents in MA';

   end % if old_efit 
% --------------------------------------------------

  otherwise
    fprintf('ERROR read_mds_geqdsk: unsupported device %s\n',tokamak)
    return;

end

% ===========================================================================
% DERIVED DATA
% ===========================================================================

for k=1:neq
% Derive efit grid:

  rg= linspace(gdata.data(k).rgrid1, gdata.data(k).rgrid1+gdata.data(k).xdim, gdata.data(k).nw)';
  zg= linspace(gdata.data(k).zmid-gdata.data(k).zdim/2,gdata.data(k).zmid+gdata.data(k).zdim/2,gdata.data(k).nh)';
  dr= (rg(end)-rg(1))/(gdata.data(k).nw-1);
  dz= (zg(end)-zg(1))/(gdata.data(k).nh-1);
   
% Create useful fluxes (even for iecurr not = 2)
% Convert flux objects to REAL units:

  psizr  = -gdata.data(k).psirz'*2*pi*sign(gdata.data(k).cpasma);
  psimag = -gdata.data(k).ssimag*2*pi*sign(gdata.data(k).cpasma);
  psibry = -gdata.data(k).ssibry*2*pi*sign(gdata.data(k).cpasma);
  psibnd = -gdata.data(k).ssibry*2*pi*sign(gdata.data(k).cpasma);

% Convert pcurrt to nh x nw array and scale to MA/m^2:

 if ~isempty(gdata.data(k).pcurrt) 
  jphi = gdata.data(k).pcurrt*1.e-6/(dz*dr);
  gdata.data(k).jphi   = jphi;
  gdata.data(k).gdef.jphi = 'current density on grid  MA/m^2';
 end
% Construct output object:

  gdata.data(k).rg= rg;
  gdata.data(k).zg= zg;
  gdata.data(k).dr= dr;
  gdata.data(k).dz= dz;
  gdata.data(k).psizr  = psizr;
  gdata.data(k).psimag = psimag;
  gdata.data(k).psibry = psibry;
  gdata.data(k).psibnd = psibnd;

  gdata.data(k).gdef.rg=  ' ';
  gdata.data(k).gdef.zg= ' ';
  gdata.data(k).gdef.dr= ' ';
  gdata.data(k).gdef.dz= ' ';
  gdata.data(k).gdef.psizr = 'true total flux on grid in Wb';
  gdata.data(k).gdef.psimag = 'axis flux in true Wb';
  gdata.data(k).gdef.psibry = ' ';
  gdata.data(k).gdef.psibnd = 'boundary flux in true Wb  (psibry also defined same)';


 gdata.data(k).fcturn = gdata.fcturn;
 gdata.data(k).turnfc = gdata.turnfc;
 gdata.data(k).fcid   = gdata.fcid;
 gdata.data(k).ecid   =  gdata.ecid;
 gdata.data(k).ecturn = gdata.ecturn;
 if(isfield(gdata,'vvid'))
    gdata.data(k).vvid = gdata.vvid;
 end
 if(isfield(gdata,'vvfrac'))
    gdata.data(k).vvfrac = gdata.vvfrac;
 end

end

if(isstr(times))	% must be 'all'
   return
elseif(length(times) > 1)
   idx = find(gdata.time >= times(1) & gdata.time <= times(2));
   gdata.data = gdata.data(idx);
   gdata.time = gdata.time(idx);
   neq = length(idx);
else
   [mm,mi] = min(abs(times-gdata.time));
   temp = gdata; clear gdata;
   gdata = temp.data(mi);
   gdata.gdef = temp.gdef;
   neq = 1;
end

