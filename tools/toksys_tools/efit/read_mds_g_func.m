function [gfile,neq,eq] = read_mds_g_func(shotnum,tree,server)
 %
%  SYNTAX:  [gfile,neq,eq] = read_mds_g_func(shotnum,tree,server)
%
%  PURPOSE:  Read geqdsk generic function, returning data in approximately the
%		same form as read_gfile_func.m.
%
%  INPUT:
%	shotnum = shot number
%	tree = mdsplus tree containing efit data
%	server = mdsplus server name
%
%  OUTPUT:
%	gfile = structure containing:
%		gdef = text description of variables
%		shotnum = shot number of equilibria	
%		time = array of times of equilibria
%		data = array of structures containing data 
%			(same data produced by read_gfile_func)
%	neq = number of equilibria found
%	eq  = equilibrium structure returned by eq_mds
%
%  RESTRICTIONS: The bzero and ecase entries of gfile.data are empty.
 
%  METHOD:
 
%  WRITTEN BY:  Mike Walker 	ON 	10/8/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)read_mds_g_func.m	1.2 07/10/09

   gfile=[];
   neq=0;
   eq=[];
   if(strcmp(upper(tree),'EFIT01') & strcmp(upper(server),'NSTX'))
      disp('ERROR read_mds_g_func: EFIT01 from NSTX not supported')
      disp('    because of missing data in the MDSplus structure.');
      wait
      return;
   end

   eq = eq_mds(shotnum,tree,server,-1);
   if(~isfield(eq,'results'))
      wait(['ERROR read_mds_g_func: No equilibria stored in ' ...
	tree ' for shot ' int2str(shotnum)])
      eq = []; return;
   end
   geq = eq.results.geqdsk;

% Some EFITs don't store measurements object.
   if(isfield(eq,'measurements'))
      meq = eq.measurements;
   end

   neq = length(geq.gtime);

   gfile.gdef=gfile_def;
   gfile.shotnum=shotnum;
   gfile.time=geq.gtime;
   for k=1:length(geq.gtime)
      gfile.data(k).shotnum = shotnum;
      gfile.data(k).time = geq.gtime(k);
      gfile.data(k).bzero = [];
      gfile.data(k).cpasma = geq.cpasma(k);
      gfile.data(k).ecase = [];
      if(exist('meq')==1)
         gfile.data(k).brsp = meq.ccbrsp(:,k);
         if(~isstr(meq.cecurr))	% won't always exist
            if(min(size(meq.cecurr))==1)	% if only 1 ohmic coil
               gfile.data(k).ecurrt = meq.cecurr(k);
            else
               gfile.data(k).ecurrt = meq.cecurr(:,k);
            end
         else
            gfile.data(k).ecurrt = [];
         end
      end
      gfile.data(k).pprime = -geq.pprime(:,k);
      gfile.data(k).ffprim = -geq.ffprim(:,k);
      gfile.data(k).fpol = geq.fpol(:,k);
      gfile.data(k).limitr = geq.limitr(k);
      gfile.data(k).nbbbs = geq.nbbbs(k);

      psizr  = -geq.psirz(:,:,k)'*2*pi*sign(geq.cpasma(k));
      psimag = -geq.ssimag(k)*2*pi*sign(geq.cpasma(k));
      psibry = -geq.ssibry(k)*2*pi*sign(geq.cpasma(k));
      psibnd = -geq.ssibry(k)*2*pi*sign(geq.cpasma(k));

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

      [pcur,jphi,cpasma]= plasma_current(psizr,gfile.data(k).pprime, ...
           gfile.data(k).ffprim,psimag,psibnd,rg,zg, ...
           geq.rbbbs(:,k),geq.zbbbs(:,k),geq.nbbbs(k));

      gfile.data(k).nh=nh;
      gfile.data(k).nw=nw;
      gfile.data(k).pcurrt =pcur;

      gfile.data(k).pres   = geq.pres(:,k);
      gfile.data(k).psirz  = geq.psirz(:,:,k);
      gfile.data(k).qpsi   = geq.qpsi(:,k);
      gfile.data(k).rbbbs  = geq.rbbbs(:,k);
      gfile.data(k).rgrid1 = geq.rgrid1(k);
      gfile.data(k).rmaxis = geq.rmaxis(k);
      gfile.data(k).rzero  = geq.rzero(k);
      gfile.data(k).ssibry = geq.ssibry(k);
      gfile.data(k).ssimag = geq.ssimag(k);
      gfile.data(k).xdim   = geq.xdim(k);
      gfile.data(k).xlim   = geq.xlim;
      gfile.data(k).ylim   = geq.ylim;
      gfile.data(k).zbbbs  = geq.zbbbs(:,k);
      gfile.data(k).zdim   = geq.zdim(k);
      gfile.data(k).zmaxis = geq.zmaxis(k);
      gfile.data(k).zmid   = geq.zmid(k);

if 0
      gfile.data(k).cc2=[eq.measurements.cecurr(1:2,k);eq.measurements.ccbrsp(:,k)];
      gfile.data(k).cc =[eq.measurements.cecurr(:,k);eq.measurements.ccbrsp(:,k)];

      gfile.data(k).jphi = jphi;
      gfile.data(k).psizr = psizr;
end

   end
