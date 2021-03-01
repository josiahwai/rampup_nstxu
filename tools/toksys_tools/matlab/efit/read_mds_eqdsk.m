function [eqdata,neq,eq,ier]= read_mds_eqdsk(shotnum,times,efit_source,tokamak,server)
 %
%  SYNTAX: [eqdata,neq,eq,ier] = read_mds_eqdsk(shotnum,times,efit_source,tokamak,server)
%
%  PURPOSE:  Read geqdsk and aeqdsk data from mdsplus, return data like
%		read_gfile_tok.m.
%	       (data needed from aeqdsk not well-defined yet)
%
%  INPUT: <default>
%	shotnum = shot number
%	times   = equilibrium time(s), one of:
%			if single float#: closest time available in mds (s)
%			if pair [t1,t2]: times between t1 and t2 (s)
%			if 'all': return all times in tree <'ALL'>
%	efit_source = source tree for EFIT (e.g. <'EFIT01'>, 'EFITRT1')
%	tokamak  = tokamak name.  <DIII-D>
%	server   = optional server name. Default is to use standard mds server for tokamak.  Use
%			'local' to get data from the GA test server (currently vidar).
%
%  OUTPUT:
%	eqdata = array of structures, each containing info from 1 eqdsk
%       (single time causes structure returned to be like read_gfile_tok; 
%		others return extended structures)
%	neq   = number of eqdsk's returned
%	eq    = entire tree containing efit equilibrium data
%	ier   = error code
% 
%  RESTRICTIONS: Only handles D3D, EAST, KSTAR, and NSTX(EFITRT only) right now.  Need
%	to add code for the other devices supported by read_gfile_tok.m.

%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	10/8/08 (generalizes read_mds_geqdsk.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Conventions:
% (1) The toksys convention is that equilibrium coil currents are always returned
%    in units of MA-turns (legacy from old codes), so for each device need to make sure
%    currents are returned in these units.
% (2) When an EFIT returns both coil and vessel currents in the brsp vector, this is
%    split into currents for coils (cc and cc2) and currents for vessel elements (vc) in the
%    output objects.

 if nargin<4
    tokamak='DIII-D';
 end
 if nargin<3
    efit_source='EFIT01';
 end
 if nargin < 2
    times = 'ALL';
 end

eqdata=[]; neq=[]; eq=[];

if(isstr(times))
   if(~strcmp(upper(times),'ALL'))
      wait('ERROR read_mds_eqdsk: invalid value for times')
      return;
   end
end

  jphi=[]; % needed to make sure this is a variable not a function

idx = findstr(':LOCAL',upper(tokamak));
if(~isempty(idx))
   tokamak = tokamak(1:idx-1);
   local_server=1;
else
   if(exist('server','var') & strcmp(server,'local'))
      local_server=1;
      clear server;
   else
      local_server=0;
   end
end

switch upper(tokamak)

% ===========================================================================
  case {'DIII-D','D3D','DIIID'}
% ===========================================================================

   if(~exist('server'))
      if(local_server)
         server = 'THOR';
      else
         server = 'DIII-D';
      end
   end

% ===========================================================================
  case {'east','EAST'}
% ===========================================================================

   if(~exist('server'))
      if(local_server)
         server = 'vidar';
      elseif shotnum < 27000
         server = 'vidar';
      else
         server = '202.127.204.12';
      end
   end

% ===========================================================================
  case {'kstar','KSTAR'}
% ===========================================================================

   if(~exist('server'))
      if(local_server)
         server = 'vidar';
      else
         server = getenv('KSTAR_PCS_MDS_SERVER');
      end
   end

% ===========================================================================
  case {'nstx','NSTX','NSTXU','nstxu'}
% ===========================================================================

   if(~exist('server'))
      if(local_server)
         server = 'vidar';
      else
         server = 'nstx';
	 %if strmatch(lower(tokamak),'nstxu'),server='nstxu';end
      end
   end

  otherwise

   fprintf('ERROR read_mds_eqdsk: unsupported device %s\n',tokamak)
   return;

end

[eqdata,neq,eq,ier] = read_mds_eq_func(shotnum,efit_source,server);
if(ier)
   return;
end
if(isempty(eqdata))
   ier=1;
   return;
end

%disp('obtained equilibria, doing unit conversions ...')
eq1 = rmfield(eqdata,'gdata');
eq1 = rmfield(eq1,'descriptions');
desc = eqdata.descriptions;
for k=1:length(eqdata.gdata)
    gdata = eqdata.gdata(k);
    eq1.gdata(k) = std_efit_units(gdata,upper(tokamak));
    eq1.time(k) = eq1.gdata(k).time;
    eq1.tms(k) = eq1.gdata(k).tms;
end
eqdata = eq1;
eqdata.descriptions = desc;
eqdata.descriptions.adata = 'A-eqdsk data';
eqdata.descriptions.time = 'times corresponding to gdata (s)';
eqdata.descriptions.tms = 'times corresponding to gdata (ms)';

% Select equilibria according to user request.

% NOTE: BELOW COULD BE PROBLEMATIC IF atime ~= gtime
if(isstr(times))	% must be 'all'
   return
elseif(length(times) > 1)
   idx = find(eqdata.time >= times(1) & eqdata.time <= times(2));
   eqdata.gdata = eqdata.gdata(idx);
   eqdata.adata = eqdata.adata(idx);
   eqdata.time = eqdata.time(idx);
   eqdata.tms = eqdata.tms(idx);
   neq = length(idx);
else
   [mm,mi] = min(abs(times-eqdata.time));
   if(mm > 100e-3)	% catch the obvious error of using wrong units on time
      wait('WARNING read_mds_eqdsk: requested time is > 10ms from nearest equilibrium.')
   end
   temp = eqdata; clear eqdata;
   eqdata = temp.gdata(mi);
   eqdata.adata = temp.adata(mi);
   eqdata.time =  temp.time(mi);
   eqdata.tms =   temp.tms(mi);
   eqdata = rmfield(eqdata,'gdef'); eqdata.gdef =  temp.gdata(mi).gdef;
   eqdata.descriptions =  temp.descriptions;
   neq = 1;
end

