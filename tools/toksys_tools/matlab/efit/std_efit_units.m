function gfile_data = std_efit_units(gdata,tokamak,options)
 %
%  SYNTAX:  gfile_data = std_efit_units(gdata,tokamak)
%
%  PURPOSE:  Do device-dependent conversions of eqdsk data to a common
%	          set of units and conventions.
%
%  INPUT:
%	   gdata = structure containing g eqdsk data
%	   tokamak = name of device
%	   options = parameters needed to customize logic for a particular device
%
%  OUTPUT:
%	   gfile_data = structure containing g eqdsk data with standard units (i.e., MA-turns)
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	8/25/11
%	(extracted from read_gfile_tok.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)std_efit_units.m	1.16 07/12/12

  gfile_data = gdata;
  if(nargin>2)
     struct_to_ws(options);
  end
  if(nargin<3 | ~isfield(options,'old_efit'))
     old_efit = 0;
  end

  switch upper(tokamak)

% ===========================================================================
  case {'DIII-D','D3D','DIIID'}
% ===========================================================================
    if(~exist('nesum','var'))
       nesum = 6; 
    end

% include d3d build parameters (see build area => must be same)
    gfile_data.fcturn = ones(1,18);
    gfile_data.turnfc = [58 58 58 58 58 55 55 58 55 58 58 58 58 58 55 55 58 55];
    gfile_data.fcid = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
%   NOTE: EFIT IS CONSTRUCTED with ecid having 1,2,3,4,5,6 numbers. Even are 1's
%         odds hook to 2
    gfile_data.ecid = [ ...			% must be same as ecdata(5,:)
			1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 ...
			1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 ...
			1 2 2 1 1 2 2 1 1 1 1 1 2 2 2 2 2 2 2 1 ...
			1 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 ...
			2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 ...
			2 2 1 1 2 2 1 1 2 2 2 2 2 1 1 1 1 1 1 1 ...
			2 2];
    gfile_data.ecturn = ones(122,1);	% required if ecid is not empty
    gfile_data.gdef.ecturn = 'individual turns in ecoil(s)';

% Construct E-coil Number of turns to make output MA-turns:
    dum= min(2,nesum); % if nesum=1 then use ONLY 1 ECOIL current

  if ~isempty(gdata.brsp) & ~isempty(gdata.ecurrt)
% Construct E-coil Number of turns to make output MA-turns:
    dum= min(2,nesum); % if nesum=1 then use ONLY 1 ECOIL current
%  find reducted ecoil set data 
      nee = length(unique(gfile_data.ecid));
      clear neturn
      for ii=1:nee
         idx = find(gfile_data.ecid==ii);
         neturn(ii,1) = sum(gfile_data.ecturn(idx));
      end

 %cc-vector with 2 e-coil segments (MA-t):
    cc2 = [gdata.ecurrt(1:dum).*neturn(1:dum);gdata.brsp]*1.e-6;

% CAUTION: we are not using the 6segment e-coil anymore so this is crippled
%cc-vector with all e-coil segs(2 or 6) (MA-t):
% 6 ecoils: ECOILA    ECOILB    E567UP    E567DN    E89DN     E89UP
% (When only 2 ecoils, ECOILA = E567UP = E89UP, ECOILB = E567DN = E89DN.)
%    cc = [gdata.ecurrt;gdata.brsp]*1.e-6;
    cc= cc2;

    gfile_data.cc     = cc;
    gfile_data.cc2    = cc2;
 
    gfile_data.gdef.cc = 'E/F coil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
    gfile_data.gdef.cc2 = 'E/F coil currents in MA-turns (for 2-segment E-coil). CAUTION: cc2 has min(2,nesum) E-coil elements (could be 1 vs std: 2)';
  end

  if(isfield(gdata,'gtime') | isfield(gdata,'time'))
     if(~isfield(gdata,'gtime') & isfield(gdata,'time'))
        gdata.gtime = gdata.time;
     end
     gfile_data.time=gdata.gtime*1e-3;
     gfile_data.tms= gdata.gtime;
  end

% ===========================================================================
  case {'EAST'}
% ===========================================================================

% Construct coil current vector:
    ncadd=2; 	% fix cc so it has extra control coils
    ncadd=0; 	% fix cc so it has extra control coils
 if ~isempty(gdata.brsp)
    cc = [gdata.brsp*1e-6;zeros(ncadd,1)];
    gfile_data.cc     = cc;
    gfile_data.cc2    = cc;
    gfile_data.gdef.cc = 'PF coil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
    gfile_data.gdef.cc2 = 'PF coil currents in MA-turns ';
 end
%    gfile_data.fcturn = [1 1 1 0.17741935 0.82258065 1 1 ...
%                         1 1 1 0.17741935 0.82258065 1 1]';
    gfile_data.fcturn = [1 1 1 44/(44+204) 204/(44+204) 1 1 ...
                         1 1 1 44/(44+204) 204/(44+204) 1 1 1/2 -1/2]';
    gfile_data.turnfc = [140 140 140 44 204 60 32 140 140 140 44 204 60 32 2 2]';
    gfile_data.fcid = [1 2 3 4 4 5 6 7 8 9 10 10 11 12 13 13];
    gfile_data.ecid = [];

% EFIT order of coils is not the same as toksys order:
    gfile_data.idx_efit_to_tok = [1 8 2 9 3 10 4 11 5 12 6 13 7 14 15 16];

  if(isfield(gdata,'gtime') | isfield(gdata,'time'))
     if(~isfield(gdata,'gtime') & isfield(gdata,'time'))
        gdata.gtime = gdata.time;
     end
     gfile_data.time = gdata.gtime;
     gfile_data.tms = gdata.gtime*1e+3;
  end

% ===========================================================================
  case {'KSTAR'}
% ===========================================================================

% select case based on gfile header ecase date which is from EFITD fortran
    ecase= gdata.ecase;
    icase = 0;
    if ~isempty(findstr('03/16/2004',ecase))
      icase = 0;		% Leuer's Original KSTAR efit's
    elseif ~isempty(findstr('01/23/2002',ecase)) % this is SABBAGH EFIT
      dum= sscanf(ecase(27:38),'%f');
      if ~isempty(dum)
        shot= dum(1);
        if length(dum)>=2
          time= dum(2);
        end
        if shot<=3000 % ? guess where is shot transition from 2009 to 2010?
          icase = 1;		% Sabbagh efit 2008/2009
          nfcoil= 14;
          nves=    56;
          nesum=    1;
        else
          icase = 3;		% Sabbagh efit 2010
          nfcoil= 14;
          nves=    60;
          nesum=    1;
        end % if shot
      end % if ~isempty(dum)
    elseif (gdata.nh==33 & gdata.nw==33 & length(gdata.brsp)==32) % rtEFIT v7b (2011)
      icase = 5;
      nfcoil = 18;
      nves = 70; % Like sabbagh, vv is listed as Fcoil
      nesum = 1;
   elseif  (gdata.nh==33 & gdata.nw==33 & length(gdata.brsp)==42) % rtEFIT v7d2 (2012 active version w/ more passive els)
      icase= 6; 
      nfcoil = 18;
      nves = 180;
   elseif (gdata.nh==33 & gdata.nw==33 & length(gdata.brsp)==37) % rtefit_v7e (2012 active version, series coils and slightly shifted divertor indices)
      icase = 7;
      nfcoil = 13;
      nves = 180;
    end
      
% =============
    switch icase
% =============
       case 0	% Original Leuer EFIT
         ncadd=4; 	% fix cc so it has extra control coils
         if ~isempty(gdata.brsp)
           cc = [gdata.brsp*1e-6;zeros(ncadd,1)];
           cc2= cc;
           gfile_data.cc  = cc;
           gfile_data.cc2 = cc2;
           gfile_data.gdef.cc = 'PF coil currents in MA-turns ';
           gfile_data.gdef.cc2 = 'PF coil currents in MA-turns ';
         end
         gfile_data.fcturn = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
         gfile_data.turnfc = [180 144 72 108 208 128 72 180 144 72 108 208 128 72 6 4 6 4];
         gfile_data.fcid = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
         gfile_data.ecid = [];


% =============
       case 1	% Sabbagh efit 2009 => should use tok_sys.config_name= sab09
%       %  PF Order is: % f1u:7u,7l,6l,5l...1l, must add 4 zero current IC's IC1U..
         ncadd=4; 	% add in extra ic coils with zero current
         gfile_data.fcid=   1:(nfcoil+4); % add in 4 IC's	    
         gfile_data.turnfc= [180 144 72 108 208 128 72 72 128 208 108 72 144 180 6 4 6 4];
         gfile_data.fcturn=  ones(1,nfcoil+4);

       % VV has 56elements broken into 6 cur groups. We expand 6 actual currents to 56 elements	    
          %  load /home/leuer/tokamaks/kstar/make/sab/vvid.mat
          %  disp(int2str(vvid'))
	  %  load /home/leuer/tokamaks/kstar/make/sab/vvfrac.mat
          vvid= [...
             1  1  2  2  2  2  2  2  3  3  3  3  3  4  4  5  5  5  5  5  6  6  6  6  6  6  1  1 ...
             1  1  2  2  2  2  2  2  3  3  3  3  3  4  4  5  5  5  5  5  6  6  6  6  6  6  1  1]';
	  gfile_data.vvfrac= ones(nves,1);

         if ~isempty(gdata.brsp)
       %  EFIT brsp is terminal current (A) => Changing cc and cc2 are Mamp-turns
	    cc_terminal = [gdata.brsp(1:nfcoil)*1.e-6; zeros(ncadd,1)] ; % 1st 14 brsp => F-coils [A]
            cc= cc_terminal.*gfile_data.turnfc'; % convert to MA-turns

       % VV has 56elements broken into 6 cur groups. We expand 6 actual currents to 56 elements	    
	    proj= proj_turn(vvid,gfile_data.vvfrac);
            vc6_terminal=  gdata.brsp(nfcoil + (1:6)); % 6 currents after nfcoil are VV groups [A]
	    vc= proj'*vc6_terminal*1e-6; % expand 6 => 56 and change to MA-turn

            if any( abs(gdata.brsp((nfcoil+6+1):end)) > 1e-5*max(abs(gdata.brsp)) )
	       wait(' % CAUTION: some Sabbagh EFIT brsp(21:60) have substantial current???')
	    end

            gfile_data.cc= cc;
            gfile_data.cc2= cc;
            gfile_data.gdef.cc= 'PF coil currents in MA-turns; sab09; 14F & 4IC=0 ';
            gfile_data.gdef.cc2= 'PF coil currents in MA-turns ';
            gfile_data.vc= vc;
            gfile_data.gdef.vc =  'vessel currents in MA; sab09; 6 groups => 56elm; 1 turn per elm.';
         end
         gfile_data.vvid= 1:nves; % NOTE: we are expanding 6=>56 so keep all VV elements
         gfile_data.gdef.vvid =  'VV identifier; Note: sab09 has 6 VV  groups but we maintain 56elm';
         gfile_data.gdef.vvfrac =  'VV fraction used in Green generation; sab09 WILL CHANGE LATER ';

% =============
       case 2			% KI You efit 2009
          ncadd=4; 	% fix cc so it has extra control coils holding place of IC
          if ~isempty(gdata.brsp)
             gfile_data.cc = [gdata.brsp*1e-6;zeros(ncadd,1)];
             gfile_data.cc2 = gfile_data.cc;
          end

% =============
       case 3			% Sabbagh efit 2010
         %  PF Order is: % f1u:7u,7l,6l,5l...1l, must add 4 zero current IC's IC1U..
         ncadd=4; 	% add in extra ic coils with zero current
         %      Need to map Sab which is clockwise from inner midplane to toksysKSTAR which is CW top and CCW bottom

         id= [1 2 3 4 5 6 7 14 13 12 11 10 9 8]; % EFIT/Sab PF coil convert to sab10.mat

         gfile_data.fcid=   1:(nfcoil+4); % add in 4 IC's	   

         % NOTE: Some how the conventions for turnfc and fcturn seems opposite of EFIT!!!!!!
         % gfile order gfile_data.turnfc= [180 144 72 108 208 128 72 72 128 208 108 72 144 180 6 4 6 4]';
         gfile_data.turnfc= [180 144 72 108 208 128 72 180 144 72 108 208 128 72 6 4 6 4]';
         gfile_data.fcturn= ones(1,length(gfile_data.fcid));

         % VV see /home/leuer/tokamaks/kstar/efit/sab10/
         vvid0= [...
           1 1 2 2 2 3 3 4 4 5 5 5 5 6 6 7 7 8 8 8 8 9 9 10 10 11 11 11 12 12 ... % Outer VV
           1 1 2 2 2 3 3 4 4 5 5 5 5 6 6 7 7 8 8 8 8 9 9 10 10 11 11 11 12 12 ]'; % Inner VV
         vvid= [vvid0; (13:36)']; % add in 10 passive plates and 14 inner limiter = 24 elements
         gfile_data.vvid= vvid;

         vvfrac0= [...  
           0.2500    0.2500    0.1667    0.1667    0.1667    0.2500    0.2500 ... 
           0.2500    0.2500    0.1250    0.1250    0.1250    0.1250    0.2500 ... 
           0.2500    0.2500    0.2500    0.1250    0.1250    0.1250    0.1250 ... 
           0.2500    0.2500    0.2500    0.2500    0.1667    0.1667    0.1667 ... 
           0.2500    0.2500]'; % this is Outer VV 

         vvfrac0= [vvfrac0; vvfrac0]; % add in Inner VV which is same

         vvfrac= [vvfrac0; ones(24,1)]; % add in 10 PP + 14 IL

         if ~isempty(gdata.brsp)
           %  EFIT brsp in SABBAGH EFIT is terminal current (A) => Changing cc and cc2 are Mamp-turns

           id= [1 2 3 4 5 6 7 14 13 12 11 10 9 8]; % EFIT/Sab PF coil convert to sab10.mat
           cc_terminal = [gdata.brsp(id)*1.e-6; zeros(ncadd,1)] ; % 1st 14 brsp => F-coils [MA]
           cc= cc_terminal.*gfile_data.turnfc; % convert to MA => MA-turns

           % VV has 60elements broken into 12 cur groups. We expand 12 actual currents to 60 elements	    
           proj= proj_turn(vvid0,vvfrac0);
           vc6_terminal=  gdata.brsp(nfcoil + (1:12)); % 12 currents after nfcoil are VV groups [A]
           vc0= proj'*vc6_terminal*1e-6; % expand 12 => 60 and change to MA-turn
           vc= [vc0; zeros(24,1)]; % add in zeros for 10PP and 12IL
           vvid= (1:length(vvid))'; % since we have expaned up to all elements keep all
           gfile_data.gdef.vvid =     'VV id; We expand 12 groups in 60 elemnts +24 extra so 1:84';
           vvfrac= ones(size(vvfrac)); % since we have expaned up fractons are one
           gfile_data.gdef.vvid =  'VV id; We expand 12 groups in 60 elemnts +24 extra so 1:84';

           if any( abs(gdata.brsp((nfcoil+12+1):end)) > 1e-5*max(abs(gdata.brsp)) )
             wait(' % CAUTION: some Sabbagh EFIT brsp(15:74) have substantial current???')
           end

           gfile_data.cc= cc; 
           gfile_data.cc2= cc;
           gfile_data.gdef.cc= 'PF coil currents in MA-turns; sab09; 14F & 4IC=0 ';
           gfile_data.gdef.cc2= 'PF coil currents in MA-turns ';
           gfile_data.vc= vc; % example of vv is vc(1)= brsp(15)*1e-6*[vvturn(1) == 0.25]
           gfile_data.gdef.vc =  'vessel currents in MA; sab10; 12 groups => 60elm;';

           vvid= (1:length(vvid))'; % since we have expaned up to all elements keep all
           gfile_data.gdef.vvid =  'VV id; We expand 12 groups in 60 elemnts +24 extra so 1:84';

         end
         gfile_data.vvfrac= vvfrac;
         gfile_data.vvid= vvid;  % NOTE: we are expanding 12=>56+24 so keep all VV elements
         gfile_data.gdef.vvid =   'VV identifier; Note: sab10 has 12 VV  groups but we maintain 60elm';
         gfile_data.gdef.vvfrac = 'VV fraction used in Green generation; sab10 WILL CHANGE LATER ';

       case 4			% KI You efit 2010

       case 5  % 2011 rtEFIT
%         This EFIT seems to be different than SAB2010 in that:
%         1) PF ordering is like DIII-D so we do not have to reorder to get to TokSys (DIII-D) ordering
%         2) Inside Coils (IC) are included with the coils (14PF + 4IC = 18 currents)
%         3) Only the passive plates (10 elements) are added VV; BUT NOT Center Limiter (14 elements)
%         Below"if 0" allows VC ouput of EFIT (reduced) size elements or EFUNT (full size) elements

          ncadd = 0;
          id = [1:18]; % Fcoils
          % PF coils
          gfile_data.fcturn = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
          gfile_data.turnfc = [180 144 72 108 208 128 72 180 144 72 108 208 128 72 6 4 6 4];
          gfile_data.fcid = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
          gfile_data.ecid = [];
          gfile_data.fcturn= ones(1,length(gfile_data.fcid)); 
          if ~isempty(gdata.brsp)
             cc = [gdata.brsp(id)*1e-6;zeros(ncadd,1)];
             cc2= cc;
             gfile_data.cc  = cc;
             gfile_data.cc2 = cc2;
             gfile_data.gdef.cc = 'PF coil currents in MA-turns ';
             gfile_data.gdef.cc2 = 'PF coil currents in MA-turns ';
          end
          
          % Vessel = 60 VV -> 12 sections + 10 PP -> 2 sections
          vvid= [...
             1 1 2 2 2 3 3 4 4 5 5 5 5 6 6 7 7 8 8 8 8 9 9 10 10 11 11 11 12 12 ... % Outer VV
             1 1 2 2 2 3 3 4 4 5 5 5 5 6 6 7 7 8 8 8 8 9 9 10 10 11 11 11 12 12 ... % Inner VV
             13 13 13 14 14 14 13 13 14 14]'; % Passive Plates
          vvfrac0= [...  
            0.2500    0.2500    0.1667    0.1667    0.1667    0.2500    0.2500 ... 
            0.2500    0.2500    0.1250    0.1250    0.1250    0.1250    0.2500 ... 
            0.2500    0.2500    0.2500    0.1250    0.1250    0.1250    0.1250 ... 
            0.2500    0.2500    0.2500    0.2500    0.1667    0.1667    0.1667 ... 
            0.2500    0.2500]'; % this is Outer VV 
          vvfrac0= [vvfrac0; vvfrac0]; % add in Inner VV which is same
          vvfrac= [vvfrac0; 0.2*ones(10,1)]; % add in 10 PP elements
          vvfrac = [vvfrac; 
                 [0.250, 0.250, 0.250, 0.250, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200,...
                 0.200, 0.200, 0.200, 0.200, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050,...
                 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050,...
                 0.050, 0.050, 0.050, 0.050, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034,...
                 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034,...
                 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034, 0.034,...
                 0.034, 0.034, 0.034, 0.062, 0.062, 0.062, 0.062, 0.062, 0.062, 0.062,...
                 0.062, 0.062, 0.062, 0.062, 0.062, 0.062, 0.062, 0.062, 0.062, 0.040,...
                 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040,...
                 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040,...
                 0.040, 0.040, 0.040, 0.040]'];

          vvid(71:174) = 0; % These elements do not exist in 2011 rtefit_v7b 

          vc14_terminal=  gdata.brsp(nfcoil + (1:14)); % 14 currents after nfcoil are VV groups [A]
          vc = vc14_terminal*1e-6; % Keep EFIT size (14) vc=> 1:14

          gfile_data.vc= vc; % example of vv is vc(1)= brsp(15)*1e-6*[vvturn(1) == 0.25]
          gfile_data.vvid = vvid; % NOTE: we are expanding 14=>70 so keep all VV elements
          gfile_data.vvfrac = vvfrac; % Expanded out, so all ones
         
          gfile_data.gdef.vc =  'vessel currents in MA; 14 groups => 70elm;';
          gfile_data.gdef.vvid =   'VV identifier';
          gfile_data.gdef.vvfrac = 'VV fraction used in Green generation';
        
        case 6  %  rtEFIT v7d2 (2012 intermediate version)
          % Similar to previous rtEFIT, except includes inner limiter, cryo, and additional
          % anti-parallel circuit to model eddy currents in passive plate. This means
          % 180 vessel elements instead of 70. 
          ncadd = 0;
          id = [1:18]; % Fcoils
          % PF coils
          gfile_data.fcturn = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
          gfile_data.turnfc = [180 144 72 108 208 128 72 180 144 72 108 208 128 72 6 4 6 4];
          gfile_data.fcid = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
          gfile_data.ecid = [];
          gfile_data.fcturn= ones(1,length(gfile_data.fcid)); 
          if ~isempty(gdata.brsp)
             cc = [gdata.brsp(id)*1e-6;zeros(ncadd,1)];
             cc2= cc;
             gfile_data.cc  = cc;
             gfile_data.cc2 = cc2;
             gfile_data.gdef.cc = 'PF coil currents in MA-turns ';
             gfile_data.gdef.cc2 = 'PF coil currents in MA-turns ';
          end
          
          % Vessel = 180 VV -> 12 vessel sections + 2 inner lim + 2 diverter + 2 symmetrical
          % passive plates + 4 cryo sections + 2 antiparallel pp eddy circuits = 24 sections
          vvid= [...
              1  1  2  2  2  3  3  4  4  5  5  5  5  6  6  7  7  8  8  8 ...
              8  9  9 10 10 11 11 11 12 12  1  1  2  2  2  3  3  4  4  5 ...
              5  5  5  6  6  7  7  8  8  8  8  9  9  10 10 11 11 11 12 12 ...
              13 13 13 14 14 14 15 15 16 16 15 15 16 16 15 15 16 16 17 17 ...
              17 18 18 18 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 ...
              19 19 19 19 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 ...
              20 20 20 20 20 20 20 20 20 20 20 20 20 21 21 21 21 21 21 21 ...
              21 21 21 21 21 21 21 21 21 22 22 22 22 22 22 22 22 22 22 22 ...
              22 22 22 22 22 22 22 22 22 22 22 22 22 22 23 23 23 24 24 24]'; %rtefit_v7e
           vvfrac= [...  
                    0.2500  0.2500  0.1667  0.1667  0.1667  0.2500  0.2500  0.2500  0.2500  0.1250 ...
                    0.1250  0.1250  0.1250  0.2500  0.2500  0.2500  0.2500  0.1250  0.1250  0.1250 ...
                    0.1250  0.2500  0.2500  0.2500  0.2500  0.1667  0.1667  0.1667  0.2500  0.2500 ...
                    0.2500  0.2500  0.1667  0.1667  0.1667  0.2500  0.2500  0.2500  0.2500  0.1250 ...
                    0.1250  0.1250  0.1250  0.2500  0.2500  0.2500  0.2500  0.1250  0.1250  0.1250 ...
                    0.1250  0.2500  0.2500  0.2500  0.2500  0.1667  0.1667  0.1667  0.2500  0.2500 ...
                    0.3333  0.3333  0.3333  0.3333  0.3333  0.3333  0.1667  0.1667  0.1667  0.1667 ...
                    0.1667  0.1667  0.1667  0.1667  0.1667  0.1667  0.1667  0.1667  0.3333  0.3333 ...
                    0.3333  0.3333  0.3333  0.3333  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500 ...
                    0.0500  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500 ...
                    0.0500  0.0500  0.0500  0.0500  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345 ...
                    0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345 ...
                    0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345 ...
                    0.0345  0.0345  0.0345  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625 ...
                    0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0400 ...
                    0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400 ...
                    0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400 ...
                    0.0400  0.0400  0.0400  0.0400  0.500   0.500  -1.000  -1.000   0.500   0.500]';
          idxvv_efit_to_tok = [1:length(vvfrac)-6]; % toksys model does not have anti-series PP

          vc =  gdata.brsp(nfcoil + (1:24))*1e-6; % 24 currents after nfcoil are VV groups [MA-turn]
          
          gfile_data.vc= vc; % example of vv is vc(1)= brsp(15)*1e-6*[vvturn(1) == 0.25]
          gfile_data.vvid = vvid; 
          gfile_data.vvfrac = vvfrac; 
          gfile_data.idxvv_efit_to_tok = idxvv_efit_to_tok;
         
          gfile_data.gdef.vc =  'vessel currents in MA; 24 groups => 180 elm;';
          gfile_data.gdef.vvid =   'VV identifier';
          gfile_data.gdef.vvfrac = 'VV fraction used in Green generation';

        case 7  %  rtEFIT v7e ()
           % Similar to rtEFIT_v7d, but with series PF coils and slightly shifted divertor indices 
           ncadd = 0;
           id = [1:13]; % Fcoils
           % PF coils
           gfile_data.fcid   = [ 1  2  3  4  5  6  7  1  2  8  9 10 11  7 12 13  12  13];
           gfile_data.fcturn = [.5 .5  1  1  1  1 .5 .5 .5  1  1  1  1 .5 .5 .5 -.5  .5]';
           gfile_data.turnfc = [360 288 72 108 208 128 144 72 108 208 128 12 8];

           gfile_data.ecid = [];
           if ~isempty(gdata.brsp)
              cc = [gdata.brsp(id)*1e-6;zeros(ncadd,1)];
              cc2= cc;
              gfile_data.cc  = cc;
              gfile_data.cc2 = cc2;
              gfile_data.gdef.cc = 'PF coil currents in MA-turns ';
              gfile_data.gdef.cc2 = 'PF coil currents in MA-turns ';
           end

           % Vessel = 180 VV -> 12 vessel sections + 2 inner lim + 2 diverter + 2 symmetrical
           % passive plates + 4 cryo sections + 2 antiparallel pp eddy circuits = 24 sections
           vvid= [...
              1  1  2  2  2  3  3  4  4  5  5  5  5  6  6  7  7  8  8  8 ...
              8  9  9 10 10 11 11 11 12 12  1  1  2  2  2  3  3  4  4  5 ...
              5  5  5  6  6  7  7  8  8  8  8  9  9  10 10 11 11 11 12 12 ...
              13 13 13 14 14 14 15 15 16 16 15 15 16 16 15 15 16 16 17 17 ...
              17 18 18 18 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 19 ...
              19 19 19 19 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 ...
              20 20 20 20 20 20 20 20 20 20 20 20 20 21 21 21 21 21 21 21 ...
              21 21 21 21 21 21 21 21 21 22 22 22 22 22 22 22 22 22 22 22 ...
              22 22 22 22 22 22 22 22 22 22 22 22 22 22 23 23 23 24 24 24]'; %rtefit_v7e
           vvfrac= [...  
                    0.2500  0.2500  0.1667  0.1667  0.1667  0.2500  0.2500  0.2500  0.2500  0.1250 ...
                    0.1250  0.1250  0.1250  0.2500  0.2500  0.2500  0.2500  0.1250  0.1250  0.1250 ...
                    0.1250  0.2500  0.2500  0.2500  0.2500  0.1667  0.1667  0.1667  0.2500  0.2500 ...
                    0.2500  0.2500  0.1667  0.1667  0.1667  0.2500  0.2500  0.2500  0.2500  0.1250 ...
                    0.1250  0.1250  0.1250  0.2500  0.2500  0.2500  0.2500  0.1250  0.1250  0.1250 ...
                    0.1250  0.2500  0.2500  0.2500  0.2500  0.1667  0.1667  0.1667  0.2500  0.2500 ...
                    0.3333  0.3333  0.3333  0.3333  0.3333  0.3333  0.1667  0.1667  0.1667  0.1667 ...
                    0.1667  0.1667  0.1667  0.1667  0.1667  0.1667  0.1667  0.1667  0.3333  0.3333 ...
                    0.3333  0.3333  0.3333  0.3333  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500 ...
                    0.0500  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500  0.0500 ...
                    0.0500  0.0500  0.0500  0.0500  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345 ...
                    0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345 ...
                    0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345  0.0345 ...
                    0.0345  0.0345  0.0345  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625 ...
                    0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0625  0.0400 ...
                    0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400 ...
                    0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400  0.0400 ...
                    0.0400  0.0400  0.0400  0.0400  0.500   0.500  -1.000  -1.000   0.500   0.500]';
 
           idxvv_efit_to_tok = [1:length(vvfrac)-6]; % toksys model does not have anti-series PP

           vc =  gdata.brsp(nfcoil + (1:24))*1e-6; % 24 currents after nfcoil are VV groups [MA-turn]

           gfile_data.vc= vc; % example of vv is vc(1)= brsp(15)*1e-6*[vvturn(1) == 0.25]
           gfile_data.vvid = vvid; 
           gfile_data.vvfrac = vvfrac; 
           gfile_data.idxvv_efit_to_tok = idxvv_efit_to_tok;

           gfile_data.gdef.vc =  'vessel currents in MA; 24 groups => 180 elm;';
           gfile_data.gdef.vvid =   'VV identifier';
           gfile_data.gdef.vvfrac = 'VV fraction used in Green generation';

       otherwise
          wait(['ERROR read_gfile_tok: unsupported KSTAR efit format ' gfile_data.ecase])
          return;
    end

    gfile_data.ecid = [];
    if(isfield(gdata,'gtime') | isfield(gdata,'time'))
       if(~isfield(gdata,'gtime') & isfield(gdata,'time'))
          gdata.gtime = gdata.time;
       end
       gfile_data.time = gdata.gtime;
       gfile_data.tms = gdata.gtime*1e+3;
    end

% ===========================================================================
  case {'NSTXU'}
% ===========================================================================

  if(~exist('nfcoil','var'))
     nfcoil = [];
  end
  if(~exist('nesum','var'))
     nesum = [];
  end
  if(~exist('nves','var'))
     nves = [];
  end

  icase = 'rtefit_v7a';
  
  switch icase
     
    case 'lrdfit'   %used for initial model development
     
     if isempty(nfcoil),   nfcoil=14; end  % default NSTXU values
     if isempty(nves),     nves=12; end     % default NSTXU values
     if isempty(nesum),    nesum=1; end    % default NSTXU values
     nec=6;	   %number of ohmic coil elements
     
     % Believe turnfc should be tok_data_struct.fcnturn  
     if(~isfield(gfile_data,'turnfc'))
       %This turnfc consistent with Menard Nucl. Fusion 52 (2012) 083015
       gfile_data.turnfc=[64 32 20 14 14 7 8 7 8  4  5  8  4  5 ...
 			   8 12 12 12 12 7 8 7 8 14 14 20 32 64];
       %This turnfc tweaked so build_tokamak_system doesn't complain
       %about PF3U,L coil currents being different (line 562)
       %Change is to increase turns to 32 instead of 30
       %Eg turnfc=8 8 8 8 instaed of 7 8 7 8
       % (there must be a better way to do this)
       %gfile_data.turnfc=[64 32 20 14 14 8 8 8 8  4  5  8  4  5 ...
       %		   8 12 12 12 12 8 8 8 8 14 14 20 32 64];
     end
     
     % Group 28 VF coils into 14 circuits
     if(~isfield(gfile_data,'fcid'))
       gfile_data.fcid=[1 2 3 4 4 5 5 5 5 6 6 6 7 7 7 8 8 9 9 10 10 10 10 11 11 12 13 14];
     end
     % Specify fraction of current in each coil for the specified fcid
     % Hence, if you do k = find(fcid==2) % where 2 is an example
     % then sum(fcturn(k)) % should be 1
     if(~isfield(gfile_data,'fcturn'))
 	 for k=1:nfcoil
 	   idfc = gfile_data.fcid(k);
 	   tmp = sum(gfile_data.turnfc(gfile_data.fcid==idfc)); % get total turns according to fcid
 	   gfile_data.fcturn(k)=gfile_data.turnfc(k)/tmp;    % apply fcturn definition
 	 end	   
     end

   % vv 
     gfile_data.vvid = [1:nves];
     gfile_data.vvfrac = ones(nves,1);
   
   % ec
   % To make the object have a consistent interpretation, this needs to be the
   % fraction of current carried by each defined group????
   %	gfile_data.ecturn = gfile_data.ecturn/sum(gfile_data.ecturn);
     if(~isfield(gfile_data,'ecid'))
       gfile_data.ecid=ones(59,1);
     end  
     if(~isfield(gfile_data,'ecturn'))
       gfile_data.ecturn=repmat(15,[59 1]);
     end   
     gfile_data.gdef.ecturn = 'individual turns in ecoil(s)';
   
   % Changing cc and cc2 so they are amp-turns
   % brsp is terminal current (A) NOT Amp-turns because NSTX uses groups and turns
   % for each group. 
   
   % for LRDFIT equilibria, brsp and ecurrt are NaN.
   
     if length(gdata.brsp)==1 & isnan(gdata.brsp)
       fprintf('WARNING std_efit_units: brsp is NaN, coil currents cannot be processed\n')
       gfile_data.vc  = [];
       gfile_data.cc  = [];
       gfile_data.cc2 = [];
     elseif length(gdata.ecurrt)==1 & isnan(gdata.ecurrt)
       fprintf('WARNING std_efit_units: ecurrt is NaN, coil currents cannot be processed\n')
       gfile_data.vc  = [];
       gfile_data.cc  = [];
       gfile_data.cc2 = [];
     elseif (~isempty(gdata.brsp) & ~isempty(gdata.ecurrt) & length(gdata.brsp)>=(nfcoil+nves))
       gfile_data.vc = gdata.brsp(nfcoil+1:nfcoil+nves)*1.e-6;
       cc_terminal = [gdata.ecurrt; gdata.brsp(1:nfcoil)]*1.e-6;

       turns = zeros(2,1); %force to be column
       for k=1:nesum
 	  turns(k) = sum(gfile_data.ecturn);
       end
       for k=1:nfcoil
 	  idx = find(gfile_data.fcid==k);
 	  turns(k+nesum) = sum(gfile_data.turnfc(idx));
       end

       gfile_data.cc	 = cc_terminal.*turns;
       gfile_data.cc2	 = gfile_data.cc;

       gfile_data.gdef.cc =  'Fcoil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
       gfile_data.gdef.cc2 = 'Fcoil currents in MA-turns ';
       gfile_data.gdef.vc =  'vessel currents in MA';
     end


   case 'rtefit_v7a'   %first rtefit for NSTXU based on Sabbagh mhdin.dat dated 1/15/15 

     if isempty(nfcoil),   nfcoil=14; end  % default NSTXU values
     if isempty(nves),     nves=40; end     % default NSTXU values
     if isempty(nesum),    nesum=1; end    % default NSTXU values
     nec=8;	   %number of ohmic coil elements

     % 22 sub coils now so can't use with LRDFIT models
     % Order is from rtefit_NSTX_diagnostics.h
     % "PF1AU","PF1BU","PF1CU", "PF2U","PF3U","PF4U","PF5U","PF5L","PF4L","PF3L","PF2L","PF1CL","PF1BL","PF1AL",     
     gfile_data.fcid = [1.00000000E+00,   2.00000000E+00,   3.00000000E+00,   4.00000000E+00,... 
                        4.00000000E+00,   5.00000000E+00,   5.00000000E+00,   6.00000000E+00,...
                        6.00000000E+00,   7.00000000E+00,   7.00000000E+00,   8.00000000E+00,... 
                        8.00000000E+00,   9.00000000E+00,   9.00000000E+00,   1.00000000E+01,...
                        1.00000000E+01,   1.10000000E+01,   1.10000000E+01,   1.20000000E+01,... 
                        1.30000000E+01,   1.40000000E+01];
     gfile_data.turnfc = [6.40000000E+01,   3.20000000E+01,   2.00000000E+01,   1.40000000E+01,... 
                           1.40000000E+01,   1.50000000E+01,   1.50000000E+01,   9.00000000E+00,... 
                           8.00000000E+00,   1.20000000E+01,   1.20000000E+01,   1.20000000E+01,... 
                           1.20000000E+01,   9.00000000E+00,   8.00000000E+00,   1.50000000E+01,... 
                           1.50000000E+01,   1.40000000E+01,   1.40000000E+01,   2.00000000E+01,... 
                           3.20000000E+01,   6.40000000E+01];
     gfile_data.fcturn = [1, 1, 1,...
                          0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,...
			  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,...
			  1, 1, 1];			   
     gfile_data.vvid = [1, 1, 2, 2, 3, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7,...
 			8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 9, 9, 9, 9, 9, 9,...
 			9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10,...
 			11, 11, 11, 11, 11, 12, 12, 13, 13, 13, 14, 15, 16, 17, 18, 18, 18, 19, 19, 20,...
 			20, 20, 20, 20, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22,...
 			22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,...
 			24, 24, 24, 24, 25, 25, 25, 25, 25, 23, 23, 23, 23, 23, 23, 26, 26, 26, 27, 28,...
			29, 29, 30, 30, 31, 31, 32, 32, 33, 34, 35, 36, 37, 38, 39, 40]; 			
     gfile_data.vvfrac = [0.5000,  0.5000,  0.5000,  0.5000,  1.0000,  1.0000,  0.3333,  0.3333,  0.3333,  0.0625,...
 			  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,  0.2000,  0.2000,  0.2000,  0.2000,  0.2000,...
 			  0.2500,  0.2500,  0.2500,  0.2500,  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,...
 			  0.0625,  0.0625,  0.0625,  0.0625,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,...
 			  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,...
 			  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.5000,  0.5000,...
 			  0.2000,  0.2000,  0.2000,  0.2000,  0.2000,  0.5000,  0.5000,  0.3333,  0.3333,  0.3333,...
 			  1.0000,  1.0000,  1.0000,  1.0000,  0.3333,  0.3333,  0.3333,  0.5000,  0.5000,  0.2000,...
 			  0.2000,  0.2000,  0.2000,  0.2000,  0.5000,  0.5000,  0.0417,  0.0417,  0.0417,  0.0417,...
 			  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,...
 			  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,  0.0417,...
 			  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,...
 			  0.2500,  0.2500,  0.2500,  0.2500,  0.2000,  0.2000,  0.2000,  0.2000,  0.2000,  0.0625,...
 			  0.0625,  0.0625,  0.0625,  0.0625,  0.0625,  0.3333,  0.3333,  0.3333,  1.0000,  1.0000,...
 			  0.5000,  0.5000,  0.5000,  0.5000,  0.5000,  0.5000,  0.5000,  0.5000,  1.0000,  1.0000,...
 			  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000];
     gfile_data.ecid  = [1, 1, 1, 1, 1, 1, 1, 1];
     gfile_data.ecturn = [112.0,  110.0, 109.5, 108.5, 108.5, 109.5, 110.0,  112.0];

     turns = zeros(2,1); %force to be column
     for k=1:nesum
        turns(k) = sum(gfile_data.ecturn);
     end
     for k=1:nfcoil
        idx = find(gfile_data.fcid==k);
        turns(k+nesum) = sum(gfile_data.turnfc(idx));
     end
     
     %brsp has also the vessel currents
     cc_terminal = [gdata.ecurrt; gdata.brsp(1:nfcoil)]*1.e-6;
     gfile_data.cc     = cc_terminal.*turns;
     gfile_data.cc2    = gfile_data.cc;


     gfile_data.gdef.fcid    = 'Fcoil circuit id';
     gfile_data.gdef.fcturn = 'Multiplier used in generation of EFIT greens tables; if =1 brsp is in A-turn; fraction of total turns';
     gfile_data.gdef.turnfc = 'Multiplier of EFIT input currents, if =1 brsp is in terminal [A]';
     gfile_data.gdef.vvid   =  'VV identifier; Note: sab09 has 6 VV  groups but we maintain 56elm';
     gfile_data.gdef.vvfrac =  'VV fraction used in Green generation; sab09 WILL CHANGE LATER ';
     gfile_data.gdef.ecid   = 'Ohmic coil circuit id';
     gfile_data.gdef.ecturn = 'individual turns in ecoil(s)';
     gfile_data.gdef.cc     =  'Fcoil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
     gfile_data.gdef.cc2    = 'Fcoil currents in MA-turns ';
     gfile_data.gdef.vc     =  'vessel currents in MA';
     
   
   end %icase switch

     if(isfield(gdata,'gtime') | isfield(gdata,'time'))
 	  if(~isfield(gdata,'gtime') & isfield(gdata,'time'))
 	     gdata.gtime = gdata.time;
 	  end
 	  gfile_data.time = gdata.gtime;
 	  gfile_data.tms = gdata.gtime*1e+3;
      end

% ===========================================================================
  case {'NSTX'}
% ===========================================================================

  if(~exist('nfcoil','var'))
     nfcoil = [];
  end
  if(~exist('nesum','var'))
     nesum = [];
  end
  if(~exist('nves','var'))
     nves = [];
  end

  if old_efit==0 % NEW Apr 05 EFIT Config: 04202005Av1.0

    if isempty(nfcoil),   nfcoil=17; end  % default NSTX values
    if isempty(nves),     nves=35; end    % default NSTX values
    if isempty(nesum),    nesum=1; end    % default NSTX values

% Coil Note except for changes to pf1aul,
%      below I add a 2nd vector which adds from the Old to the New(Apr05)
    gfile_data.fcid= [[1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 11] ...
                      [12 13 14 15 16 17]];
    gfile_data.turnfc=[[20 14 14 15 15 9 8 12 12 12 12 9 8 15 15 14 14 20 32]...
                       [1 1 1 1 48 48]];

    gfile_data.fcturn = ...
     [[1 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1 1]...
      [1 1 1 1 1 1]];

% vv has some differences internal so new from "make"
    gfile_data.vvid = [ ...
        1 2 3 4 4 5 6 6 7 7 7 8 8 8 8 8 8 8 8 9 9 9 9 9 10 10 11 11 11 12 ...
        13 14 15 16 16 16 17 17 18 18 18 18 18 19 19 19 19 19 20 20 21 ...
        22 22 23 24 25 26 26 27 27 28 29 30 31 32 33 34 35];

    gfile_data.vvfrac = [ ...
   	1 1 1 0.5 0.5 1 ...
	0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 ...
	0.25 0.25 0.25 0.25 ...
	0.33333 0.33333 0.33333 0.33333 0.33333 0.33333 ...
	0.5 0.5 0.5 0.5 ...
	0.33333 0.33333 0.33333 0.33333 0.33333 0.33333 ...
	0.25 0.25 0.25 0.25 ...
	0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 ...
	1 0.5 0.5 1 1 1 0.5 0.5 0.5 0.5 1 1 1 1 1 1 1 1];

    gfile_data.vvfrac = [ ...
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
    gfile_data.ecid = ones(240,1);
    gfile_data.ecturn = [ ...
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
    gfile_data.gdef.ecturn = 'individual turns in ecoil(s)';

% To make the object have a consistent interpretation, this needs to be the
% fraction of current carried by each defined group????
%    gfile_data.ecturn = gfile_data.ecturn/sum(gfile_data.ecturn);
%    gfile_data.turnec = ?;  % WHAT IS THIS??????

% Changing cc and cc2 so they are amp-turns
% brsp is terminal current (A) NOT Amp-turns because NSTX uses groups and turns
% for each group. 

% for LRDFIT equilibria, brsp and ecurrt are NaN.

  if length(gdata.brsp)==1 & isnan(gdata.brsp)
    fprintf('WARNING std_efit_units: brsp is NaN, coil currents cannot be processed\n')
    gfile_data.vc  = [];
    gfile_data.cc  = [];
    gfile_data.cc2 = [];
  elseif length(gdata.ecurrt)==1 & isnan(gdata.ecurrt)
    fprintf('WARNING std_efit_units: ecurrt is NaN, coil currents cannot be processed\n')
    gfile_data.vc  = [];
    gfile_data.cc  = [];
    gfile_data.cc2 = [];
  elseif ~isempty(gdata.brsp) & ~isempty(gdata.ecurrt)
    gfile_data.vc = gdata.brsp(nfcoil+1:nfcoil+nves)*1.e-6;
    cc_terminal = [gdata.ecurrt; gdata.brsp(1:nfcoil)]*1.e-6;

    turns = zeros(2,1);	%force to be column
    for k=1:nesum
       turns(k) = sum(gfile_data.ecturn);
    end
    for k=1:nfcoil
       idx = find(gfile_data.fcid==k);
       turns(k+nesum) = sum(gfile_data.turnfc(idx));
    end

    gfile_data.cc     = cc_terminal.*turns;
    gfile_data.cc2    = gfile_data.cc;

    gfile_data.gdef.cc =  'Fcoil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
    gfile_data.gdef.cc2 = 'Fcoil currents in MA-turns ';
    gfile_data.gdef.vc =  'vessel currents in MA';
  end
  
% ----------------------------------------------------------------------------
 else % if old_efit==1 OLD Feb 2002 EFIT BELOW  Config_name= '02072002Av1.0';
% ----------------------------------------------------------------------------

    if isempty(nfcoil),   nfcoil=11; end    % default NSTX values
    if isempty(nves),     nves=30; end    % default NSTX values
    if isempty(nesum),    nesum=1; end       % default NSTX values

    gfile_data.fcid = [1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 9 9 10 11];
    gfile_data.turnfc = [48 14 14 15 15 9 8 12 12 12 12 9 8 15 15 14 14 48 32];
    gfile_data.fcturn = ...
     [1 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 1 1];

    gfile_data.vvid = [ ...
	1 2 3 4 4 5 6 6 6 6 6 6 6 7 7 7 7 8 8 8 9 9 9 10 10 11 11 12 12 12 ...
	13 13 13 14 14 14 14 15 15 15 15 15 15 15 16 17 17 18 19 20 21 21 ...
	22 22 23 24 25 26 27 28 29 30];
    gfile_data.vvfrac = [ ...
   	1 1 1 0.5 0.5 1 ...
	0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 ...
	0.25 0.25 0.25 0.25 ...
	0.33333 0.33333 0.33333 0.33333 0.33333 0.33333 ...
	0.5 0.5 0.5 0.5 ...
	0.33333 0.33333 0.33333 0.33333 0.33333 0.33333 ...
	0.25 0.25 0.25 0.25 ...
	0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 0.14286 ...
	1 0.5 0.5 1 1 1 0.5 0.5 0.5 0.5 1 1 1 1 1 1 1 1];

    gfile_data.ecid = ones(240,1);
    gfile_data.ecturn = [ ...
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
    gfile_data.gdef.ecturn = 'individual turns in ecoil(s)';

% To make the object have a consistent interpretation, this needs to be the
% fraction of current carried by each defined group????
%    gfile_data.ecturn = gfile_data.ecturn/sum(gfile_data.ecturn);
%    gfile_data.turnec = ?;  % WHAT IS THIS??????

% Changing cc and cc2 so they are amp-turns
% brsp is terminal current (A) NOT Amp-turns because NSTX uses groups and turns
% for each group. 

 if ~isempty(gdata.brsp) & ~isempty(gdata.ecurrt)
    gfile_data.vc = gdata.brsp(nfcoil+1:nfcoil+nves)*1.e-6;
    cc_terminal = [gdata.ecurrt; gdata.brsp(1:nfcoil)]*1.e-6;

    turns = zeros(2,1);	%force to be column
    for k=1:nesum
       turns(k) = sum(gfile_data.ecturn);
    end
    for k=1:nfcoil
       idx = find(gfile_data.fcid==k);
       turns(k+nesum) = sum(gfile_data.turnfc(idx));
    end

    gfile_data.cc     = cc_terminal.*turns;
    gfile_data.cc2    = gfile_data.cc;

    gfile_data.gdef.cc =  'Fcoil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
    gfile_data.gdef.cc2 = 'Fcoil currents in MA-turns ';
    gfile_data.gdef.vc =  'vessel currents in MA';
  end
% --------------------------------------------------
 end % if old_efit 
% --------------------------------------------------

  if(isfield(gdata,'gtime') | isfield(gdata,'time'))
     if(isfield(gfile_data,'gtime'))
        gfile_data.time = gdata.gtime;
        gfile_data.tms = gdata.gtime*1e+3;
     else
        gfile_data.time = gdata.time;
        gfile_data.tms = gdata.time*1e+3;
     end
  end

% ===========================================================================
  case {'ITER'}
% ===========================================================================

  if ~isempty(gdata.brsp)
    cc = [gdata.brsp*1e-6];
    gfile_data.cc     = cc;
    gfile_data.cc2    = cc;
    gfile_data.gdef.cc = 'PF coil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
    gfile_data.gdef.cc2 = 'PF coil currents in MA-turns ';
  end
    gfile_data.fcid = 1:12;

% use date of file to switch between turns in system: iter07 Vs 2010, See /u/leuer/efit/iter/
    date2010= datenum('01-Jan-2010'); % Runs prior to 2010 use old turns
    d= dir(filename);
    if datenum(d.date) < datenum('01-Jan-2010') % Runs prior to 2010 use old turns
      disp('% NOTE: read_gfile using OLD 2007 ITER TURNS: 249 106 185 169 217 425 548 548..')
      gfile_data.turnfc = [249 106 185 169 217 425 548 548 548 548 548 548];
    else
      disp('% NOTE: read_gfile using 2010 ITER TURNS: 249 115 186 170 217 459 553 553..')
      gfile_data.turnfc = [249 115 186 170 217 459 553 553 553 553 553 553];
    end   
    gfile_data.fcturn = ones(12,1);

    gfile_data.ecid = [];

%    disp('Using Build for ITER07')

% ===========================================================================
  case {'CTF'}
% ===========================================================================

%cc-vector with NO-ECOIL (MA-t):
 if ~isempty(gdata.brsp)
    cc = gdata.brsp*1.e-6;

    gfile_data.cc     = cc;
    gfile_data.cc2    = cc;
 
    gfile_data.gdef.cc =  '14 F coil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
    gfile_data.gdef.cc2 = '14 F coil currents in MA-turns ';
 end
    gfile_data.fcturn = ones(1,nfcoil);    % used to build efit greens tables
    gfile_data.turnfc = 50*ones(1,nfcoil); % Ipf*turnfc => efit greens table
    gfile_data.fcid = 1:14;
    gfile_data.ecid = [];
%    disp(['read_gfile_tok done reading CTF gfile: ',filename])

% ===========================================================================
  case {'FDF'}
% ===========================================================================

%cc-vector with NO-ECOIL (MA-t):
 if ~isempty(gdata.brsp)

% special reordering of brsp for symmetric EFIT runs:
    if size(gdata.brsp,1)==22
       sm= 0;
       mx= 0;
       for ii= 1:11 
          sm= (gdata.brsp(ii)-gdata.brsp(ii+11)).^2 + sm;
	  mx= max([mx,abs(gdata.brsp(ii)),abs(gdata.brsp(ii+1))]);
       end
       err= sqrt(sm/11)/mx; % diff is typically of order 6e-7 for symmetric
       if err <= 1e-5
          gfile_data.brsp= gfile_data.brsp([1:11 22:-1:12]); % reorder for sym efit
       end
     end

    cc = gdata.brsp*1.e-6;
    gfile_data.cc     = cc;
    gfile_data.cc2    = cc;
 
    gfile_data.gdef.cc =  'F coil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
    gfile_data.gdef.cc2 = 'F coil currents in MA-turns ';
 end
    gfile_data.fcturn = ones(1,nfcoil);    % used to build efit greens tables
  %Apparently turnfc must be ones, even though MHDIN.DAT claims = 50; 
  %  ==> Logic in rzrig or here may be broken...   DAH 5/31/08
% Below depends if MATLAB is build with fcnturn=1 or 50 JAL 6may2011
    gfile_data.turnfc = turnfc0*ones(1,nfcoil); % Ipf*turnfc => efit greens table
    gfile_data.fcid = 1:nfcoil;
    gfile_data.ecid = [];
%    disp(['read_gfile_tok done reading FDF gfile: ',filename])


% ===========================================================================
  case {'CFETR'}
% ===========================================================================

  if ~isempty(gdata.brsp)
    cc = [gdata.brsp*1e-6];
    gfile_data.cc     = cc;
    gfile_data.cc2    = cc;
    gfile_data.gdef.cc = 'PF coil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
    gfile_data.gdef.cc2 = 'PF coil currents in MA-turns ';
  end
    nfc=14;
    gfile_data.fcid = 1:nfc;
    gfile_data.turnfc = ones(nfc,1);
    gfile_data.fcturn = ones(nfc,1);
    gfile_data.ecid = [];
    
    disp('cfetr')

% ===========================================================================
  case {'PEGASUS'}
% ===========================================================================

 if ~isempty(gdata.brsp)
%   ncadd=4;    % fix cc so it has extra control coils
%   cc = [gdata.brsp*1e-6;zeros(ncadd,1)];
%   cc2= cc;

%   gfile_data.cc  = cc;
%   gfile_data.cc2 = cc2;
%   gfile_data.gdef.cc = 'PF coil currents in MA-turns ';
%   gfile_data.gdef.cc2 = 'PF coil currents in MA-turns ';
 end

% ???  Do we care that length(cc)~=length(turnfc) ???

%    gfile_data.fcturn = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
    gfile_data.turnfc = [2 3 5 5 5 5 3 2 5 5 280 4 4 1 1 1 1];
    gfile_data.fcid = [1 2 3 4 4 3 2 5 6 7 8 9 10 11 12 13 14];
    gfile_data.ecid = [];

    gfile_data.fcturn = [1 0.5 0.5 0.5 0.5 0.5 0.5 1 1 1 1 1 1 1 1 1 1]';

% ===========================================================================
  case {'SST'}
% ===========================================================================

    gfile_data.fcturn = ones(1,14);
    gfile_data.turnfc = [40 40 192 40 40 8 1 40 40 192 40 40 8 1];
    gfile_data.fcid = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
    gfile_data.ecid = [ ...                     % must be same as ecdata(5,:)
                        1 1 1 1 1 ];
    gfile_data.ecturn = [692 40 8 40 8]; % required if ecid is not empty
    gfile_data.gdef.ecturn = 'individual turns in ecoil(s)';

% Construct E-coil Number of turns to make output MA-turns:
    dum= min(2,nesum); % if nesum=1 then use ONLY 1 ECOIL current

  if ~isempty(gdata.brsp) & ~isempty(gdata.ecurrt)
% Construct E-coil Number of turns to make output MA-turns:
    dum= min(2,nesum); % if nesum=1 then use ONLY 1 ECOIL current
%  find reduced ecoil set data
      nee = length(unique(gfile_data.ecid));
      clear neturn
      for ii=1:nee
         idx = find(gfile_data.ecid==ii);
         neturn(ii,1) = sum(gfile_data.ecturn(idx));
      end

 %cc-vector with 2 e-coil segments (MA-t):
    cc2 = [gdata.ecurrt(1:dum).*neturn(1:dum);gdata.brsp]*1.e-6;

% CAUTION: we are not using the 6segment e-coil anymore so this is crippled
%cc-vector with all e-coil segs(2 or 6) (MA-t):
    cc= cc2;

    gfile_data.cc     = cc;
    gfile_data.cc2    = cc2;

    gfile_data.gdef.cc = 'E/F coil currents in MA-turns; convert to toksys cc0 using cc_efit_to_tok';
    gfile_data.gdef.cc2 = 'E/F coil currents in MA-turns ';
  end
% ===========================================================================
  case {'HL2M'} % Jim Leuer 2jul2015
% ===========================================================================
    if(~exist('nesum','var'))
       nesum = 1; 
    end

% include d3d build parameters (see build area => must be same)
    gfile_data.fcturn = [28 28 28 28 28 27 28 28 28 28 28 28 28 27 28 28];
    gfile_data.turnfc = ones(1,16);
    gfile_data.fcid =    1:16;
    gfile_data.ecid = [ 1 1];
    gfile_data.ecturn = [48 48]; % ? NOT SURE IF THIS IS 1 or 48?
    gfile_data.gdef.ecturn = 'individual turns in upper and lower ecoils';


  if ~isempty(gdata.brsp) & ~isempty(gdata.ecurrt)
% Construct E-coil Number of turns to make output MA-turns:
    dum= 1; % ONLY 1 ECOIL current
%  find reducted ecoil set data 
      nee = length(unique(gfile_data.ecid));
      clear neturn
      for ii=1:nee
         idx = find(gfile_data.ecid==ii);
         neturn(ii,1) = sum(gfile_data.ecturn(idx));
      end

 % cc-vector with 1 e-coil segment is in Amps, brsp is in Amps so conversion  (MA-t):
    cc2 = [gdata.ecurrt(1:dum).*neturn(1:dum);...
           gdata.brsp.*gfile_data.fcturn']*1.e-6; % this is now in MA-t

    cc= cc2;

    gfile_data.cc     = cc;
    gfile_data.cc2    = cc2;
 
    gfile_data.gdef.cc = 'E/F currents in MA-t; convert to toksys cc0 using std_efit_units';
    gfile_data.gdef.cc2 = 'E/F currents in MA-turns';
  end

  if(isfield(gdata,'gtime') | isfield(gdata,'time'))
     if(~isfield(gdata,'gtime') & isfield(gdata,'time'))
        gdata.gtime = gdata.time;
     end
     gfile_data.time=gdata.gtime*1e-3;
     gfile_data.tms= gdata.gtime;
  end
% Validated using test_build jal 7jul2014

% ===========================================================================
  otherwise
% ===========================================================================

  disp(['%ERROR: std_efit_units does not recognize tokamak:' tokamak])

end  % switch

  gfile_data.gdef.fcturn = 'Multiplier used in generation of EFIT greens tables; if =1 brsp is in A-turn';
  gfile_data.gdef.turnfc = 'Multiplier of EFIT input currents, if =1 brsp is in terminal [A]';
  gfile_data.gdef.fcid =   'Coil identifier used in construction of greens tables';
  gfile_data.gdef.ecid =   'Coil identifier used in construction of greens tables';

% ===========================================================================
% DERIVED DATA
% ===========================================================================

% Derive efit grid:

  rg= linspace(gdata.rgrid1, gdata.rgrid1+gdata.xdim, gdata.nw)';
  zg= linspace(gdata.zmid-gdata.zdim/2,gdata.zmid+gdata.zdim/2,gdata.nh)';
  dr= (rg(end)-rg(1))/(gdata.nw-1);
  dz= (zg(end)-zg(1))/(gdata.nh-1);
   
% Create useful fluxes (even for iecurr not = 2)
% Convert flux objects to REAL units:

  psizr  = -gdata.psirz'*2*pi*sign(gdata.cpasma);
  psimag = -gdata.ssimag*2*pi*sign(gdata.cpasma);
  psibry = -gdata.ssibry*2*pi*sign(gdata.cpasma);
  psibnd = -gdata.ssibry*2*pi*sign(gdata.cpasma);

% Convert pcurrt to nh x nw array and scale to MA/m^2:

 if ~isempty(gdata.pcurrt) & ~isfield(gfile_data,'jphi')
  jphi = gdata.pcurrt*1.e-6/(dz*dr);
  gfile_data.jphi   = jphi;
  gfile_data.gdef.jphi = 'current density on grid  MA/m^2';
 end

% Construct output object:

  gfile_data.rg= rg;
  gfile_data.zg= zg;
  gfile_data.dr= dr;
  gfile_data.dz= dz;
  gfile_data.psizr  = psizr;
  gfile_data.psimag = psimag;
  gfile_data.psibry = psibry;
  gfile_data.psibnd = psibnd;

  gfile_data.gdef.rg=  'radial coordinates of plasma grid';
  gfile_data.gdef.zg= 'vertical coordinates of plasma grid';
  gfile_data.gdef.dr= 'radial distance between plasma grid points';
  gfile_data.gdef.dz= 'vertical distance between plasma grid points';
  gfile_data.gdef.psizr = 'true total flux on grid in Wb';
  gfile_data.gdef.psimag = 'axis flux in true Wb';
  gfile_data.gdef.psibry = 'flux on plasma boundary';
  gfile_data.gdef.psibnd = 'boundary flux in true Wb  (psibry also defined same)';

