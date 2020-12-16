function [equil,neq,eq] = read_mds_eq_func1(shotnum,tree,server,verbose)
 %
%  SYNTAX:  [equil,neq,eq] = read_mds_eq_func1(shotnum,tree,server)
%
%  PURPOSE:  Generic read of geqdsk and aeqdsk data, returning data in 
%	approximately the same form as read_gfile_func.m.
%	(data needed from aeqdsk not well-defined yet)
%
%  INPUT: <default>
%	shotnum = shot number
%	tree = mdsplus tree containing efit data
%	server = mdsplus server name
%       verbose= 1; print out progressive info. during MDS READ, 0=no print <[]>
%
%  OUTPUT:
%	equil = structure containing:
%		gdef = text description of variables
%		shotnum = shot number of equilibria	
%		time = array of times of equilibria
%		gdata = array of structures containing geqdsk data 
%			(same data produced by read_equil_func)
%		adata = array of structureus containing aeqdsk data
%	neq = number of equilibria found
%	eq  = equilibrium structure returned by get_mds_tree
%
%  RESTRICTIONS: The bzero and ecase entries of equil.gdata are empty.


% *******************************************************************************
%	!!! DO NOT MODIFY THIS CODE!!!!!!
%  INSTEAD, modify the calling code to not depend on time, tms units and replace
%  the call with call to read_mds_eq_func.m.
% *******************************************************************************

 
%  METHOD:
 
%  WRITTEN BY:  Mike Walker 	ON 	10/8/08 (generalizes read_mds_g_func)
%  JAL02feb2011 Corrected using read_gfile.m as TRUTH
%  jal22apr2011 time=s, tms=ms, gtime&atime=ms(d3d);s(nstx),
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @(#)read_mds_eq_func1.m	1.4 08/06/14  

wait('ERROR read_mds_eq_func1: THIS FUNCTION IS OBSOLETE - replace with read_mds_eq_func.m')
return;

   if nargin <= 3
     verbose= [];
   end

   equil=[];
   neq=0;
   eq=[];

   eq = get_mds_tree(shotnum,tree,server,-1,verbose);
   if(~isfield(eq,'results') && ~isfield(eq,'cpasma'))
      wait(['ERROR read_mds_eq_func: No equilibria stored in ' ...
	tree ' for shot ' int2str(shotnum)])
      eq = []; return;
   end
   if strmatch(tree,'EFIT_EAST') % Special case for non-standard tree
     eq = convert_efiteast_to_eqdsk(eq);
   end
   geq = eq.results.geqdsk;
   aeq = eq.results.aeqdsk;

% Some EFITs don't store measurements object.
   if(isfield(eq,'measurements'))
      meq = eq.measurements;
   end
   
   % Offline efits on NSTX store fitted currents in aeq
   if strcmp(upper(server),'NSTX') & ~exist('meq','var')
     if isfield(aeq,'ccbrsp') & isfield(aeq,'eccurt')
       meq.ccbrsp = NaN*ones(size(aeq.ccbrsp,1),length(geq.gtime));
       if min(size(aeq.eccurt)) == 1
         meq.cecurr = NaN*ones(1,length(geq.gtime));
       else
         meq.cecurr = NaN*ones(size(aeq.eccurt,1),length(geq.gtime));
       end
       for k=1:length(geq.gtime)
	 j = find(geq.gtime(k) == aeq.atime);
	 if length(j) == 1
           meq.ccbrsp(:,k) = aeq.ccbrsp(:,j);
	   if min(size(aeq.eccurt)) == 1
             meq.cecurr(:,k) = aeq.eccurt(j);
	   else
             meq.cecurr(:,k) = aeq.eccurt(:,j);
	   end
	 end
       end
     else
       disp('ERROR read_mds_eq_func: EFIT from NSTX not supported')
       disp('    because of missing data in the MDSplus structure.');
       wait
       return;
     end
     if isfield(aeq,'eccurt') % ?? what is this for
     end
   else

   end

% START LOAD OF EQUIL
   neq = length(geq.gtime);

   equil.gdef=gfile_def;
   equil.shotnum=shotnum;

% STANDARD: time => sec, tms => ms
% NSTX (EFIT01 EFITRT) have all (gtime, atime ...) in Seconds ==> time 
% Others like DIIID (EAST?NSTX?) have gtime,atime in ms!!  ==> tms
   if strcmp(upper(server),'NSTX')
     equil.time= geq.gtime;
     equil.tms=  geq.gtime*1e+3;
   elseif strcmp(upper(server),'KSTAR') | ~isempty(findstr(upper(tree),'KSTR')) |...
          strcmp(upper(server),'EAST') | ~isempty(findstr(upper(tree),'EAST'))
     equil.time= geq.gtime;
     equil.tms=  geq.gtime*1e+3;      
   
   else % DIIID
     equil.time=geq.gtime*1e-3;
     equil.tms= geq.gtime;
   end 
   equil.gdef.time= ['Time (sec); Note: gtime,atime = ms (D3D), sec (NSTX)'];
   equil.gdef.tms= ['Time in ms; '];

   for k=1:length(geq.gtime)
      equil.gdata(k).shotnum = shotnum;
      if strcmp(upper(server),'NSTX') | ...
         strcmp(upper(server),'KSTAR') | ~isempty(findstr(upper(tree),'KSTR')) | ...
         strcmp(upper(server),'EAST') | ~isempty(findstr(upper(tree),'EAST')) 
         equil.gdata(k).gtime = geq.gtime(k);
         equil.gdata(k).time =  geq.gtime(k);
         equil.gdata(k).tms=    geq.gtime(k)*1e+3;
      else
         equil.gdata(k).gtime = geq.gtime(k);
         equil.gdata(k).time =  geq.gtime(k)*1e-3;
         equil.gdata(k).tms=    geq.gtime(k);
      end
      equil.gdata(k).bzero = geq.bcentr(k); % vacuum toridal field %jal02022011
      equil.gdata(k).cpasma = geq.cpasma(k);
      equil.gdata(k).ecase = [];
      if(exist('meq')==1)
         equil.gdata(k).brsp = meq.ccbrsp(:,k);
         if(isfield(meq,'cecurr') && ~isstr(meq.cecurr))	% won't always exist
            if(min(size(meq.cecurr))==1)	% if only 1 ohmic coil
               equil.gdata(k).ecurrt = meq.cecurr(k);
            else
               equil.gdata(k).ecurrt = meq.cecurr(:,k);
            end
         else
            equil.gdata(k).ecurrt = [];
         end
      end
      equil.gdata(k).pprime = -sign(geq.cpasma(k))*geq.pprime(:,k); %jal02022011
      equil.gdata(k).ffprim = -sign(geq.cpasma(k))*geq.ffprim(:,k); %jal02022011
      equil.gdata(k).fpol = geq.fpol(:,k);
      equil.gdata(k).limitr = geq.limitr(k);
      equil.gdata(k).nbbbs = geq.nbbbs(k);

      psizr  = -geq.psirz(:,:,k)'*2*pi*sign(geq.cpasma(k)); % NOTE different
      psimag = -geq.ssimag(k)*2*pi*sign(geq.cpasma(k));     % than read_gfile.m
      psibry = -geq.ssibry(k)*2*pi*sign(geq.cpasma(k));     % But seems correct
      psibnd = -geq.ssibry(k)*2*pi*sign(geq.cpasma(k));     % jal02022011

      if(min(size(geq.r))==1)
         rg = geq.r; 
         zg = geq.z; 
      else	% KLUGE for NSTX
         rg = geq.r(:,1); 
         zg = geq.z(:,1); 
      end
      nw = length(rg); 
      nh = length(zg); 
      rgg = ones(nh,nw)*diag(rg);
      zgg = diag(zg)*ones(nh,nw);

      equil.gdata(k).nh=nh;
      equil.gdata(k).nw=nw;
      
      if ischar(geq.rbbbs) & isfield(geq,'rbdry')
        geq.rbbbs = geq.rbdry;
        geq.zbbbs = geq.zbdry;
        geq.nbbbs = geq.nbdry;
      end
      
      if(geq.cpasma(k)~=0 & geq.nbbbs(k)>1)
         [pcur,jphi,cpasma]= plasma_current(psizr,equil.gdata(k).pprime, ...
           equil.gdata(k).ffprim,psimag,psibnd,rg,zg, ...
           geq.rbbbs(:,k),geq.zbbbs(:,k),geq.nbbbs(k));
         equil.gdata(k).pcurrt =pcur;
      else
         equil.gdata(k).pcurrt = zeros(nh,nw);
      end

      equil.gdata(k).pres = geq.pres(:,k);
      equil.gdata(k).psirz  = geq.psirz(:,:,k);
      equil.gdata(k).qpsi   = geq.qpsi(:,k);
      equil.gdata(k).rbbbs  = geq.rbbbs(:,k);
      equil.gdata(k).rgrid1 = geq.rgrid1(k);
      equil.gdata(k).rmaxis = geq.rmaxis(k);
      if isfield(geq,'rzero') equil.gdata(k).rzero  = geq.rzero(k); end
      equil.gdata(k).ssibry = geq.ssibry(k);
      equil.gdata(k).ssimag = geq.ssimag(k);
      equil.gdata(k).xdim   = geq.xdim(k);
      equil.gdata(k).xlim   = geq.xlim;
      equil.gdata(k).ylim   = geq.ylim;
      equil.gdata(k).zbbbs  = geq.zbbbs(:,k);
      equil.gdata(k).zdim   = geq.zdim(k);
      equil.gdata(k).zmaxis = geq.zmaxis(k);
      equil.gdata(k).zmid   = geq.zmid(k);

if 0
      equil.gdata(k).cc2=[eq.measurements.cecurr(1:2,k);eq.measurements.ccbrsp(:,k)];
      equil.gdata(k).cc =[eq.measurements.cecurr(:,k);eq.measurements.ccbrsp(:,k)];

      equil.gdata(k).jphi = jphi;
      equil.gdata(k).psizr = psizr;
end

      idx = find(aeq.atime == geq.gtime(k));
      if(~isempty(idx))
         equil.adata(k).rxpt1 = aeq.rxpt1(idx);
         equil.adata(k).zxpt1 = aeq.zxpt1(idx);
         equil.adata(k).rxpt2 = aeq.rxpt2(idx);
         equil.adata(k).zxpt2 = aeq.zxpt2(idx);
      else		% no Xpts computed for this time
         equil.adata(k).rxpt1 = -9.99;
         equil.adata(k).zxpt1 = -9.99;
         equil.adata(k).rxpt2 = -9.99;
         equil.adata(k).zxpt2 = -9.99;
      end

   end
