function [equil,neq,eq,ier] = read_mds_eq_func(shotnum,tree,server,verbose)
 %
%  SYNTAX:  [equil,neq,eq,ier] = read_mds_eq_func(shotnum,tree,server)
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
%	ier = error code
%
%  RESTRICTIONS: The bzero and ecase entries of equil.gdata are empty.
 
%  METHOD:
 
%  WRITTEN BY:  Mike Walker 	ON 	10/8/08 (generalizes read_mds_g_func)
%  JAL02feb2011 Corrected using read_gfile.m as TRUTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)read_mds_eq_func.m	1.17 11/23/11

   if nargin <= 3
     verbose= [];
   end

   equil=[];
   neq=0;
   eq=[];

   [eq,ier] = get_mds_tree(shotnum,tree,server,-1,verbose);
   if(isfield(eq,'results'))
      geq = eq.results.geqdsk;
      if(ischar(geq.psirz) & strncmp(geq.psirz,'Java exception',14))
         wait(['ERROR read_mds_eq_func: No equilibria stored in ' ...
            tree ' on server ' server ' for shot ' int2str(shotnum)])
         ier=1; return;
      end
      aeq = eq.results.aeqdsk;
   elseif(isfield(eq,'psirz'))
      geq = eq;
      aeq = [];
   else
      wait(['ERROR read_mds_eq_func: No equilibria stored in ' ...
	tree ' on server ' server ' for shot ' int2str(shotnum)])
      eq = []; ier=1; return;
   end
   if(~isfield(geq,'gtime') & isfield(geq,'time'))
      geq.gtime = geq.time;
   end

% Some EFITs don't store measurements object.
   if(isfield(eq,'measurements'))
      meq = eq.measurements;
   end
   
% Some put coil currents in aeq.
   if exist('aeq','var') & ~exist('meq','var')
     if( isfield(aeq,'ccbrsp') & ~isstr(aeq.ccbrsp) & ~isempty(aeq.ccbrsp) & ...
	 isfield(aeq,'eccurt') & ~isstr(aeq.eccurt) & ~isempty(aeq.eccurt) )
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
     elseif( isfield(geq,'ccbrsp') & ~isstr(geq.ccbrsp) & ~isempty(geq.ccbrsp) ) 
       meq.ccbrsp = geq.ccbrsp;
     else
       wait('WARNING read_mds_eq_func: missing coil current data in EFIT.')
     end

   end

% START LOAD OF EQUIL
   neq = length(geq.gtime);

   equil.shotnum=shotnum;
   equil.server = server;
   equil.tree = tree;
   equil.time=geq.gtime;

   for k=1:length(geq.gtime)
      equil.gdata(k).shotnum = shotnum;
      equil.gdata(k).time = geq.gtime(k);
      equil.gdata(k).bzero = geq.bcentr(k); % vacuum toridal field %jal02022011
      equil.gdata(k).cpasma = geq.cpasma(k);
      equil.gdata(k).ecase = [];
      if(exist('meq')==1)
         equil.gdata(k).brsp = meq.ccbrsp(:,k);
         if(isfield(meq,'cecurr') & ~isstr(meq.cecurr))	% won't always exist
            if(numel(meq.cecurr)==length(geq.gtime))	% conclude only 1 ohmic coil
               %rtefit writes 1D vector of length ntimes if only one ohmic coil
               %offline efit always writes 2D vector with dimensions nOHcoilxntimes
               %handle both cases 
               tmp = reshape(meq.cecurr,length(meq.cecurr),1);
               equil.gdata(k).ecurrt = tmp(k);
            else
               equil.gdata(k).ecurrt = meq.cecurr(:,k);
            end
         else
            equil.gdata(k).ecurrt = [];
         end
      else
         equil.gdata(k).brsp = [];
         equil.gdata(k).ecurrt = [];
      end
      equil.gdata(k).pprime = -sign(geq.cpasma(k))*geq.pprime(:,k); %jal02022011
      equil.gdata(k).ffprim = -sign(geq.cpasma(k))*geq.ffprim(:,k); %jal02022011
      equil.gdata(k).fpol = geq.fpol(:,k);
      equil.gdata(k).limitr = geq.limitr(k);
      if(isfield(geq,'nbbbs'))
         equil.gdata(k).nbbbs = geq.nbbbs(k);
      elseif(isfield(geq,'nbdry'))
         equil.gdata(k).nbbbs = geq.nbdry(k);
      end

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
      
% If *bbbs objects don't exist or are error messages, use the *bdry objects:
      if(~isfield(geq,'rbbbs') | ischar(geq.rbbbs))
        if(isfield(geq,'rbdry'))
           geq.rbbbs = geq.rbdry;
           geq.zbbbs = geq.zbdry;
           geq.nbbbs = geq.nbdry;
        elseif(isfield(geq,'bdry'))
           geq.nbbbs = geq.nbdry;
           geq.rbbbs = squeeze(geq.bdry(1,:,:));
           geq.zbbbs = squeeze(geq.bdry(2,:,:));
        end
      end
      
      if(geq.cpasma(k)~=0 & geq.nbbbs(k)>1)
         [pcur,jphi,cpasma]= plasma_current(psizr,equil.gdata(k).pprime, ...
           equil.gdata(k).ffprim,psimag,psibnd,rg,zg, ...
           geq.rbbbs(:,k),geq.zbbbs(:,k),geq.nbbbs(k));
         equil.gdata(k).pcurrt =pcur;
      else
         equil.gdata(k).pcurrt = zeros(nh,nw);
      end

      equil.gdata(k).pres   = geq.pres(:,k);
      equil.gdata(k).psirz  = geq.psirz(:,:,k);
      equil.gdata(k).qpsi   = geq.qpsi(:,k);
      equil.gdata(k).rbbbs  = geq.rbbbs(:,k);
      equil.gdata(k).rgrid1 = geq.rgrid1(k);
      equil.gdata(k).rmaxis = geq.rmaxis(k);
      equil.gdata(k).rzero  = geq.rzero(k);
      equil.gdata(k).ssibry = geq.ssibry(k);
      equil.gdata(k).ssimag = geq.ssimag(k);
      equil.gdata(k).xdim   = geq.xdim(k);
      if(isfield(geq,'xlim'))
         equil.gdata(k).xlim   = geq.xlim;
         equil.gdata(k).ylim   = geq.ylim;
      elseif(isfield(geq,'lim'))
         equil.gdata(k).xlim   = geq.lim(1,:);
         equil.gdata(k).ylim   = geq.lim(2,:);
      end
      equil.gdata(k).zbbbs  = geq.zbbbs(:,k);
      equil.gdata(k).zdim   = geq.zdim(k);
      equil.gdata(k).zmaxis = geq.zmaxis(k);
      equil.gdata(k).zmid   = geq.zmid(k);
      equil.gdata(k).gdef   = gfile_def;

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
         equil.adata(k).betap = aeq.betap(idx);
         equil.adata(k).li    = aeq.li(idx);
      else		% no Xpts computed for this time
         equil.adata(k).rxpt1 = -9.99;
         equil.adata(k).zxpt1 = -9.99;
         equil.adata(k).rxpt2 = -9.99;
         equil.adata(k).zxpt2 = -9.99;
         equil.adata(k).betap = NaN;
         equil.adata(k).li    = NaN;
      end

   end
   equil.descriptions.shotnum = 'shot number from which data extracted';
   equil.descriptions.server = 'mds server from which data extracted';
   equil.descriptions.time= ['Time (units vary for different devices) '];
   equil.descriptions.adata = 'array of A-eqdsk structures';
   equil.descriptions.gdata = 'array of B-eqdsk structures';
