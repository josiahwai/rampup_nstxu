function [model_data,cc0,vc0,dbg_objs]=rzrig( ...
	equil_data,tokamak,vac_objs,idoplots,idoncal,idxvv,iwait,iverbose)
 %
%  USAGE:  function [model_data,cc0,vc0,dbg_objs]=rzrig( ...
%		equil_data,tokamak,vac_objs,idoplots,idoncal,idxvv,iwait)
%
%  PURPOSE: Script to calculate vertical and radial response using rigid
%           current-conserving model of plasma for GENERAL systems.
%	    Takes basic data from EFIT g-file or external data
%	    input (cc, jphi; eg from Corsica) and needs environment objects of
%	    same names as D3D environment. cc and jphi are ASSUMED to be in
%	    MA when read from gfile or defined externally.
%
%  INPUTS:
%    equil_data = structure containing equilibrium information
%    tokamak = device to construct model objects for (e.g. 'NSTX','KSTAR',etc)
%    vac_objs = structure containing: 
%		mcc, mvv, mcv, mpc, mpv, resc, resv, zg, rg, ecdata
%     		(Note that imks and iterminal in this structure define the
%		 units of the data objects to be produced. imks=1 gives MKS
%		 units, otherwise units are MA,uH,uOhms.  iterminal=1 gives
%		 terminal mode, 0 gives lumped.)
%    idoplots= (optional) flag to select plotting: 1=plot, 0=don't(default)
%    idoncal = (optional) flag to select calc of decay index: 1=calc, 
%			0=don't(default)
%    idxvv = (optional) indices of VV elements to use in force calculation 
%			(default=1:nvv)
%    iwait = (optional) 1(default)=wait with error messages, 0=no wait
%    iverbose = (optional) Flag to display messages, default is 0
%
%  OUTPUTS:
%    model_data = quantities describing vertical/radial force balance&couplings
%                       - data objects
%                       - units = units of data objects
%                       - desc  = descriptions of data objects
%    cc0 = equilibrium coil currents (size matching objects in vac_objs, with
%	   units as specified by imks and iterminal in vac_objs)
%    vc0 = equilibrium vessel currents (size matching objects in vac_objs, with
%          units as specified by imks and iterminal in vac_objs)
%
%  Also displays max positive eigenvalue for various conductor modifications.
%
%  RESTRICTIONS:
%     Plasma current density must be zero within 1 grid of upper & lower walls
%       (for 2-sided calculations of gradients).
%     Circuits defaulted to be lumped (one-turn) elements (see ccnturn def 
%	below). If magnetic mapping objects (mcc, mpc, etc...) are in 
%	terminal mode (iterminal==1), ccnturn must correspond to turns in them.
%     cc coming from read_gfile (or wherever) assumed to be lumped mode.
%
%  METHOD:  
%     Described in Walker/Humphreys, Valid Coordinate Systems for Linearized
%     Plasma Shape Response Models in Tokamaks, FS&T,Nov.06; GA report GA-A25042

%  WRITTEN BY:  Dave Humphreys  ON	8/10/03
%
%  MODIFICATION HISTORY:
%    DAH  8/10/03  Assembling machine-independent version of rrigd3d from
%		   zrig.m and rrigd3d.m.
%    DAH  8/21/03  Making rzresp.m out of present rzrig.m, which is to be
%			called from calc_rzresp to connect with d3d_sim
%			model-building tools.
%    DAH  1/27/04  Modified dFrdr to have correct hoop force expression 
%			(see rzresp.m in ~humphrys/matlab/plresp).
%    ASW  4/30/04  Going back to old hoop force expression, adding drdli
%                       nesum is 1 by default unless given.
%    DAH  7/29/04  Adding imks flag to do calculations for MKS inputs, 
%			produce MKS in all output quantities.  
%    ASW  5/13/09  Adding MSE signals and loop voltage loops
%	
%  NOTES:
%   7/30/04 Must still put in terminal mode option, clean up usage of cc, turns
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)rzrig.m	1.22 09/09/11

include_cross_terms = 0;

   imks = vac_objs.imks;
   iterminal = vac_objs.iterminal;

   mcc = vac_objs.mcc;
   mvv = vac_objs.mvv;
   mcv = vac_objs.mcv;
   mpc = vac_objs.mpc;
   mpv = vac_objs.mpv;
   mpl = vac_objs.mpl;
   gpb = vac_objs.gpb;
   resc = vac_objs.resc;
   resv = vac_objs.resv;
   zg = vac_objs.zg;
   rg = vac_objs.rg;
   if(isfield(vac_objs,'ecdata'))
      ecdata = vac_objs.ecdata;
   else
      ecdata = [];
   end

% Get dimensions:
   sss=size(mvv); nvv=sss(1);
   sss=size(mcc); ncc=sss(1);
   nss=ncc+nvv;  %total # of conductors PF + VV
   nz=length(zg); nr=length(rg);  %grid dimensions
   ngg=nz*nr;

% Inputs and defaults:
   if nargin<6, idxvv=1:nvv; end, nvvx=length(idxvv);
   if exist('idoplots')~=1, idoplots=0; end
   if exist('idoncal')~=1, idoncal=0; end
   if exist('rzrespeq')~=1, rzrespeq=1; end %flag to enable reads of equilib
   if exist('iwait')~=1, iwait=1; end %default to waits enabled
   if exist('iverbose')~=1, iverbose=0; end

% Prelims:
   mu0 = 4*pi*1e-7;
   twopi = 2*pi;
   dz=abs(zg(2)-zg(1));  
   dr=abs(rg(2)-rg(1));
   rgg = ones(nz,1)*rg';
   zgg = zg*ones(1,nr);
   rggv = rgg(:);
   zggv = zgg(:);
   twopir = twopi*rggv;

%Derived values:
   if imks, iscale=1e6; else iscale=1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equilibrium reading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if rzrespeq  %only if eq read enabled, extract EFIT data from equil_data:
      equil_I = cc_efit_to_tok(vac_objs,equil_data);
      cc0 = equil_I.cc0t;
      vc0 = equil_I.vc0t;
   end   %end if rzrespeq

  if(abs(equil_data.dz-dz)>1e-8)
    if iwait==1
     wait(['WARNING rzrig: tokamak vacuum objects dz= ' num2str(dz) ...
		', not equal to EFIT dz=' num2str(equil_data.dz)])
    elseif iwait==0
     disp(['WARNING rzrig: tokamak vacuum objects dz= ' num2str(dz) ...
		', not equal to EFIT dz=' num2str(equil_data.dz)])
    end
  end
  if(abs(equil_data.dr-dr)>1e-8)
    if iwait==1
     wait(['WARNING rzrig: tokamak vacuum objects dr= ' num2str(dr) ...
		', not equal to EFIT dr=' num2str(equil_data.dr)])
    elseif iwait==0
     disp(['WARNING rzrig: tokamak vacuum objects dr= ' num2str(dr) ...
		', not equal to EFIT dr=' num2str(equil_data.dr)])
    end
  end

%Jphi-related variables:
  jphi0 = iscale*equil_data.jphi;	%jphi0 now imks
  cphi = jphi0*equil_data.dz*equil_data.dr;  %cphi now imks
  ip0 = sum(cphi(:));			%ip0 now imks
  if(ip0==0)
     wait('ERROR rzrig: plasma current jphi=0 - plasma response calculation is invalid')
     return;
  end
  idxpl = find(cphi(:)~=0);

% Extract "shape" quantities from cphi:
   a0 = (max(rggv(idxpl))-min(rggv(idxpl)))/2;
   b0 = (max(zggv(idxpl))-min(zggv(idxpl)))/2;
   if(a0==b0)   % handle "point" distributions of current
      kap0=1;
   else
      kap0 =  b0/a0;
   end
   z0 = sum(jphi0(:).*zggv)./sum(jphi0(:));  %Zcentroid
   r0 = sum(jphi0(:).*rggv)./sum(jphi0(:));  %centroid major radius

% Plot current density and equilbrium field:
  if idoplots
   figure(1),clf,hold on
   lev=linspace(min(min(jphi0)),max(max(jphi0)),30)';
   contour(rg,zg,jphi0,lev,'-')
   hold on
   title('Jphi, Psi(vac) contours')
   if(isempty(vc0))
      psit = reshape(mpc*cc0,nz,nr);
   else
      psit = reshape(mpc*cc0+mpv*vc0,nz,nr);
   end
   contour(rg,zg,psit,20,'--')
   axis('image')
  end

%----------------------------Vertical response--------------------------
% Calculate growth rate using zrigd3d algorithm for Matlab:

 %ZSHIFT algorithm: DFSDZ and DFZDIS
if iverbose
   disp(['Derived from cphi, Ip = ',num2str(sum(sum(cphi)))])
   if exist('equil_data')==1 %cpasma assumed in Amps
     if imks   %compare Amps
      disp([' From EFIT cpasma, Ip = ',num2str(equil_data.cpasma)])
     else	    %compare MA
      disp([' From EFIT cpasma, Ip = ',num2str(equil_data.cpasma*1e-6)])
     end
   end
end
   dcdz = ([zeros(1,nr);cphi(1:nz-1,:)]-[cphi(2:nz,:);zeros(1,nr)])./(2*dz);
   dfvdz = mpv'*dcdz(:);
   dfcdz = mpc'*dcdz(:);
   dfsdz = [dfcdz;dfvdz];
   mss = [[mcc mcv];[mcv' mvv]];

   if exist('mpl')==1 & ~isempty(mpl)
     dfldz = mpl'*dcdz(:);    %needed for building diagnostic response to dZ
   else
     dfldz = 0;
   end
   if exist('gpb')==1 & ~isempty(gpb)
     dbpdz = gpb'*dcdz(:);
   else
     dbpdz = 0;
   end

% Differential force Fz (units = N/A for MKS, else N/MA):
% *** NOTE that dFzdis is NOT EQUAL to dfsdz except in mks units.***

   dFzdis = dfsdz'*1e6/iscale; 	%dFz is the diff. force Fz

 %Decay Index Force term (summation from Ief):
   mpc_zz=zeros(ngg,ncc);
   for ii=1:ncc
     tmp=reshape(mpc(:,ii),nz,nr);
     tmp = (  [zeros(1,nr);tmp(1:nz-1,:)] - 2*tmp ...
                       + [tmp(2:nz,:);zeros(1,nr)])./(dz^2);
     mpc_zz(:,ii) = tmp(:);
   end
   if(~isempty(vc0))
      mpv_zz=zeros(ngg,nvv);
      for ii=1:nvv
        tmp=reshape(mpv(:,ii),nz,nr);
        tmp = (  [zeros(1,nr);tmp(1:nz-1,:)] - 2*tmp ...
                       + [tmp(2:nz,:);zeros(1,nr)])./(dz^2);
        mpv_zz(:,ii) = tmp(:);
      end
   end

 %Decay Index: 
     mpc_r=zeros(ngg,ncc);
     for ii=1:ncc
       tmp=reshape(mpc(:,ii),nz,nr);
       tmp = ( -[zeros(nz,1) tmp(:,1:nr-1)] ...
                         + [tmp(:,2:nr) zeros(nz,1)] )./(2*dr);
       mpc_r(:,ii) = tmp(:);
     end
     bz = reshape((mpc_r*cc0)./twopir,nz,nr);

     if(~isempty(vc0))
        mpv_r=zeros(ngg,nvv);
        for ii=1:nvv
          tmp=reshape(mpv(:,ii),nz,nr);
          tmp = ( -[zeros(nz,1) tmp(:,1:nr-1)] ...
                         + [tmp(:,2:nr) zeros(nz,1)] )./(2*dr);
          mpv_r(:,ii) = tmp(:);
        end
        bz = bz+reshape((mpv_r*vc0)./twopir,nz,nr);
     end

 %  if idoncal   %commenting out 8/25/07: always calculate now...
     bz0 = interp2(rg,zg,bz,r0,z0);
     bzav = sum(bz(:).*cphi(:))/ip0;
     dbrdz = -(mpc_zz*cc0)./twopir;
     if(~isempty(vc0))
        dbrdz = dbrdz - (mpv_zz*vc0)./twopir;
     end
     ndig = -rggv.*dbrdz./bz0;      %decay index over grid
     ndecay0 = interp2(rg,zg,reshape(ndig,nz,nr),r0,z0);
     ndecay_av = sum(ndig.*cphi(:))/ip0;
     ndecay = ndecay_av;   %n/ncrit== ndecay/ncrit;
     ncrit = trace( inv(mss)*(dfsdz*dfsdz') )/(2*pi*bz0*ip0);
     if iverbose
     disp(['ndecay = ',num2str(ndecay)])
     disp(['ncrit = ',num2str(ncrit)])
     disp(['n/ncrit = ',num2str(ndecay/ncrit)])
     end
 %  end

   if(isempty(vc0))
      dFzdz = sum(cphi(:).*(mpc_zz*cc0))*1e6/iscale;   %N/m
   else
      dFzdz = sum(cphi(:).*(mpc_zz*cc0+mpv_zz*vc0))*1e6/iscale;   %N/m
   end
   br = reshape(vac_objs.gbr2c*cc0 + vac_objs.gbr2v*vc0,nz,nr);

 %Vertical Response:
   dzdis = -dFzdis./dFzdz;

%----------------------------Radial response--------------------------
 %RSHIFT algorithm: DFSDR and DFRDIS
   tmp = cphi;
% derivative of plasma currents on grid w.r.t. Rc
   dcdr = ( [zeros(nz,1) tmp(:,1:nr-1)] ...
                    - [tmp(:,2:nr) zeros(nz,1)] )./(2*dr);  %was -,+ <2/24/98

% Be careful to get the sign right for dFzdr only if significant size.
   dFzdr = -2*pi*sum(sum(rgg.*dcdr.*br))*1e6/iscale;

   dfvdr = mpv'*dcdr(:);
   dfcdr = mpc'*dcdr(:);
   dfsdr = [dfcdr;dfvdr];
   if exist('mpl')==1 & ~isempty(mpl)
     dfldr = mpl'*dcdr(:);    %needed for building diagnostic response to dR
   else 
     dfldr = 0;
   end
   if exist('gpb')==1 & ~isempty(gpb)
     dbpdr = gpb'*dcdr(:);
   else
     dbpdr = 0;
   end
   if isfield(vac_objs,'gmsebrp') & ~isempty(vac_objs.gmsebrp)
     dBrMSEdr = vac_objs.gmsebrp*dcdr(:);
     dBzMSEdr = vac_objs.gmsebzp*dcdr(:);
     dBrMSEdz = vac_objs.gmsebrp*dcdz(:);
     dBzMSEdz = vac_objs.gmsebzp*dcdz(:);
   else
     dBrMSEdr = [];
     dBzMSEdr = [];
     dBrMSEdz = [];
     dBzMSEdz = [];
   end
   if isfield(vac_objs,'mph') & ~isempty(vac_objs.mph)
     dflvdr = vac_objs.mph'*dcdr(:);
     dflvdz = vac_objs.mph'*dcdz(:);
   else
     dflvdr = [];
     dflvdz = [];
   end

% Here, dfsdr = [mpc_r'; mpv_r']*cphi(:) - see valid coords paper

   dFrdis = dfsdr'*1e6/iscale;    %dfr is the diff. force Fr, so units are N/MA
   dbzdr = gradient(bz,rg,zg);
   dFrdr_1 = twopi*sum(sum(cphi.*bz))*1e6/iscale;
   dFrdr_2 = twopi*sum(sum(rgg.*cphi.*dbzdr))*1e6/iscale;
   dFrdr_3 = mu0*ip0^2/(2*r0) * (1e6/iscale)^2;
   dFrdr = twopi*sum(sum(cphi.*bz))*1e6/iscale ... 
	+ twopi*sum(sum(rgg.*cphi.*dbzdr))*1e6/iscale  ...
        + mu0*ip0^2/(2*r0) * (1e6/iscale)^2;

% Be careful to get the sign right for dFrdz only if significant size.
   dFrdz = twopi*sum(sum(rgg.*dcdz.*bz))*1e6/iscale;

% start TEST
% from paper (noting that dcdr is derivative w.r.t Rc):
%   dFrdr = twopi*sum(sum(rgg.*dcdr.*bz))*1e6/iscale ...
%	+ mu0*ip0^2/(2*r0)*(1e6/iscale)^2;

% from calc_rzresp:
%   dFrdr = twopi*sum(sum(cphi.*bz))*1e6/iscale ... 
%	+ twopi*sum(sum(rgg.*cphi.*dbzdr))*1e6/iscale  ...
%	+ twopi*sum(sum(rgg.*cphi.*bz))*1e6/r0/iscale;
% end TEST

   drdis = -dFrdis/dFrdr;

%betap=0.69; li = 0.91;
%Gamma = log(8*r0/(a0*sqrt(kap0)))+(2*kap0/(kap0^2+1))*betap+(li/2)-1.5;
%Frhoop = mu0*ip0^2/2 * Gamma;
%Frvac = twopi*sum(sum(rgg.*cphi.*bz));

%---------------------------------------------------------------------
%Ip response:
  %dFrdip = (sum(twopir.*bz(:)) - 2*sum(twopir.*bz(:)))*1e6; %Fr/Ip+dFhp/dIp   
  %Frhoop = (mu0*ip0^2*1e6/2)*(log(8*r0/(a0*sqrt(kap0)))+betap00+li00/2-1.5);
   Frbz = twopi*sum(sum(rgg.*cphi.*bz))*1e6/iscale;   %rad force due to Bz [N]
   dFrdip = -Frbz/ip0;   %dFr/dIp = (Frbz+2Frhoop)/Ip0 = -Frbz/Ip0 
   dFzdip = 0;
   drdip = -dFrdip/dFrdr;
   dzdip = 0;

if 0	% for testing only right now
   dR_Iseq = drdis * [equil_I.cc0t; equil_I.vc0t];
   drdip = - dR_Iseq/ip0;
   dZ_Iseq = dzdis * [equil_I.cc0t; equil_I.vc0t];
   dzdip = - dZ_Iseq/ip0;
   iverbose=1;
end

   dfsdip = [mpc'*cphi(:)/ip0; mpv'*cphi(:)/ip0];

%---------------------------------------------------------------------
%Betap response:
   dFrdbetap = (2*kap0/(1+kap0^2))*(mu0*ip0^2)/2 * (1e6/iscale)^2;
   dFzdbetap = 0;
   drdbetap = -dFrdbetap/dFrdr;
   dzdbetap = 0;

%---------------------------------------------------------------------
%li response:
%   drdli = ***motion from force change + flattening motion in pl frame***
%   dzdli = ***flattening only***
    dFrdli=(mu0*ip0^2)/4*(1e6/iscale)^2; % works if li only changes hoop force!
    dFzdli= 0;
    drdli = -dFrdli/dFrdr;
    dzdli = 0;

% compare calculation with cross-terms:

if(include_cross_terms)
   drdis0 = drdis; dzdis0 = dzdis;
   drdip0 = drdip; dzdip0 = dzdip;
   drdbetap0 = drdbetap; dzdbetap0 = dzdbetap;
   drdli0 = drdli; dzdli0 = dzdli;
   temp = -inv([dFrdr dFrdz; dFzdr dFzdz])*[dFrdis; dFzdis];
   temp = -inv([dFrdr dFrdz; dFzdr dFzdz]) ...
	 *[dFrdis dFrdip dFrdbetap dFrdli;
	   dFzdis dFzdip dFzdbetap dFzdli];
   n=size(temp,2);
   drdis = temp(1,1:n-3);
   dzdis = temp(2,1:n-3);
   drdip = temp(1,n-2);
   dzdip = temp(2,n-2);
   drdbetap = temp(1,n-1);
   dzdbetap = temp(2,n-1);
   drdli = temp(1,n);
   dzdli = temp(2,n);

   figure(100),clf
   plot(drdis0)
   hold on
   plot(drdis,'r--')
   plot(dzdis0)         
   plot(dzdis,'m--')
   hold off
   title('drdis(r) and dzdis(m) with (dashed) and without (solid) cross-terms')
   xlabel('conductor number')
   ylabel('m/A')

   fprintf('\t\tdrd*0\t\tdrd*\t\tdzd*0\t\tdzd*\n')
   fprintf('ip\t%e\t%e\t%e\t%e\t\n',drdip0,drdip,dzdip0,dzdip)
   fprintf('betap\t%e\t%e\t%e\t%e\t\n',drdbetap0,drdbetap,dzdbetap0,dzdbetap)
   fprintf('li\t%e\t%e\t%e\t%e\t\n',drdli0,drdli,dzdli0,dzdli)

   wait
% NOTE: need to add corresponding cross-terms to Ip circuit eqn
end

% Verify force balance:
   test1 = drdis(1:length(cc0))*cc0;
   if(~isempty(vc0))
      test1 = test1 + drdis(length(cc0)+[1:length(vc0)])*vc0;
   end
   test = (test1 + drdip*ip0)/test1;
   if(abs(test)>1e-12) & (iverbose | iwait)
     fprintf('WARNING rzrig: drdi equation force balance is violated (%f)\n',...
				test)
     if iwait==1, wait, elseif iwait==0, end
   end
   test1 = dzdis(1:length(cc0))*cc0;
   if(~isempty(vc0))
      test1 = test1 + dzdis(length(cc0)+[1:length(vc0)])*vc0;
   end
   test = (test1 + dzdip*ip0)/test1;
   if(abs(test)>1e-12) & (iverbose | iwait)
     fprintf('WARNING rzrig: dzdi equation force balance is violated (%f)\n',...
				test)
     if iwait==1, wait, elseif iwait==0, end
   end


%---------------------------------------------------------------------
 %Eigenproblems for system w/o Ip dynamic equation:
   rvv = diag(resv);
   rss = diag([resc;resv]);
   mss = [[mcc mcv];[mcv' mvv]];

   xmatz = dfsdz*dzdis;
   xmatr = dfsdr*drdis;

   if(include_cross_terms)
      xmatz0 = dfsdz*dzdis0;
      xmatr0 = dfsdr*drdis0;
      xmatz = dfsdz*dzdis;
      xmatr = dfsdr*drdis;
   end

   lstar = mss + xmatz + xmatr;
   lstari = inv(lstar);

  %Stability margin (Portone form; cf NF 45 (2005) 926):
   tmp = -inv(mss)*lstar;
   d = eigsort(tmp);
   m_s = d(1);
   amat = -lstari*rss;
   [vecs,vals]=eigsort(amat);
   gamma = real(vals(1));
   if iverbose
     disp(['Stability margin = ',num2str(m_s)])
     disp(['High passive root = ',num2str(real(vals(1))),' + i', ...
                                	num2str(imag(vals(1)))])
     disp(' ')
     disp(['Gamma_z of unmodified system =',num2str(gamma), ...
				  ' , Tau_z=',num2str(1/gamma)])
     disp(' ')
   end

% objects needed to add Ip state:
% cc is always in MA-turns, mpc is in Wb/A when imks=iterminal=1
% ?????????????????????????????????????????????????????????????
% cc0 gives length 20 for d3d, while cc gives length 24.  What do 
% these correspond to for EAST and KSTAR?
% ?????????????????????????????????????????????????????????????
   psivac = mpc*cc0;   %vacuum flux 
   if(~isempty(vc0))
      psivac = psivac + mpv*vc0;   %vacuum flux 
   end
   dfpdr = sum(dcdr(:).*psivac)/ip0;   %gradient in (vac)flux linked by plasma
   dfpdz = sum(dcdz(:).*psivac)/ip0;
   dfpdrC = -dfpdr; % sign correction to get dfpdrC, which is what's needed
   dfpdzC = -dfpdz; % sign correction to get dfpdrC, which is what's needed

desc = struct( ...
'drdis', 'partial derivative of centroid Rc w.r.t. "stabilizing" elts', ...
'dzdis', 'partial derivative of centroid Zc w.r.t. "stabilizing" elts', ...
'dFrdis', 'change in radial force w.r.t. conductor currents', ...
'dFzdis', 'change in vertical force w.r.t. conductor currents', ...
'drdip', 'partial derivative of centroid Rc w.r.t. Ip', ...
'dzdip', 'partial derivative of centroid Zc w.r.t. Ip', ...
'dFrdip', 'partial derivative of radial force w.r.t. Ip', ...
'dfsdip', 'change in flux at conductors w.r.t. Ip', ...
'drdbetap', 'partial derivative of centroid Rc w.r.t. betap', ...
'dzdbetap', 'partial derivative of centroid Zc w.r.t. betap', ...
'dFrdbetap', 'partial derivative of radial force w.r.t. betap', ...
'drdli', 'partial derivative of centroid Rc w.r.t. li', ...
'dFrdli', 'partial derivative of radial force w.r.t. li', ...
'dFrdr', 'partial derivative of radial force w.r.t. Rc', ...
'dFzdr', 'partial derivative of vertical force w.r.t. Rc', ...
'dfsdr', 'partial derivative of flux at stabilizing elts w.r.t. Rc', ...
'dfcdr', 'partial derivative of flux at coils w.r.t. Rc', ...
'dfvdr', 'partial derivative of flux at vessel elts w.r.t. Rc', ...
'dcdr', 'partial derivative of plasma current (on grid) w.r.t. Rc', ...
'dfldr', 'change in flux at flux loops due to change in Rc', ...
'dflvdr', 'change in *flux* at loop voltage loops due to change in Rc', ...
'dbpdr', 'change in field at B probes due to change in Rc', ...
'dBrMSEdr', 'change in radial field at MSE points due to change in Rc', ...
'dBzMSEdr', 'change in vertical field at MSE points due to change in Rc', ...
'dFrdz', 'partial derivative of radial force w.r.t. Zc', ...
'dFzdz', 'partial derivative of vertical force w.r.t. Zc', ...
'dfsdz', 'partial derivative of flux at stabilizing elts w.r.t. Zc', ...
'dfcdz', 'partial derivative of flux at coils w.r.t. Zc', ...
'dfvdz', 'partial derivative of flux at vessel elts w.r.t. Zc', ...
'dfpdrC', 'partial derivative of flux at plasma w.r.t. Zc', ...
'dfpdzC', 'partial derivative of flux at plasma w.r.t. Zc', ...
'dcdz', 'partial derivative of plasma current (on grid) w.r.t. Zc', ...
'dfldz', 'change in flux at flux loops due to change in Zc', ...
'dflvdz', 'change in *flux* at loop voltage loops due to change in Zc', ...
'dbpdz', 'change in field at B probes due to change in Zc', ...
'dBrMSEdz', 'change in radial field at MSE points due to change in Zc', ...
'dBzMSEdz', 'change in vertical field at MSE points due to change in Zc', ...
'dbrdz', 'change in radial field due to change in Zc', ...
'bz0', 'vertical field at curr centroid ref pt r0, z0', ...
'bzav', 'cphi weighted av of vacuum vertical field', ...
'cphi', 'current in individual grid elements', ...
'ip0','total equilibrium current', ...
'r0','current weighted centroid radial position (m)', ...
'z0','current weighted centroid vertical position (m)', ...
'a0',' ', ...
'b0',' ', ...
'kap0','elongation of ??', ...
'mss','vacuum mutual inductance matrix', ...
'xmatz','contribution to mutual inductance due to plasma z motion', ...
'xmatr', 'contribution to mutual inductance due to plasma r motion', ...
'ndecay', 'decay index (avg over plasma): compare with ncrit', ...
'ncrit', 'critical decay index: ndecay/ncrit=1 corr to ideal instab', ...
'm_s', 'Portone stability margin: #1 eigenvalue of -inv(mss)*lstar' );

if(imks)
   dcdz_units = 'Amp/m';
   dFzdis_units = 'N/Amp';
   dFrdis_units = 'N/Amp';
   dzdis_units = 'm/Amp';
   ip0_units = 'Amp';
   cphi_units = 'Amp';
   dzdis_units = 'm/Amp';
   drdis_units = 'm/Amp';
else
   dcdz_units = 'MAmp/m';
   dFzdis_units = 'N/MAmp';
   dFrdis_units = 'N/MAmp';
   dzdis_units = 'm/MAmp';
   ip0_units = 'MAmp';
   cphi_units = 'MAmp';
   dzdis_units = 'm/MAmp';
   drdis_units = 'm/MAmp';
end
if(iterminal)
else
end

units = struct( ...
'drdis',drdis_units, ...
'dzdis',dzdis_units, ...
'dFrdis',dFrdis_units, ...
'dFzdis',dFzdis_units, ...
'drdip','??', ...
'dzdip','??', ...
'dFrdip','N/A', ...
'dfsdip','??', ...
'drdbetap','m/unit(betap)', ...
'dzdbetap','m/unit(betap)', ...
'dFrdbetap','N/unit(betap)', ...
'drdli','m/unit(li)', ...
'dFrdli','N/unit(li)', ...
'dFrdr','N/m', ...
'dFzdr','N/m', ...
'dfsdr','??', ...
'dfcdr','??', ...
'dfvdr','N/m', ...
'dcdr','??', ...
'dfldr','??', ...
'dflvdr','??', ...
'dbpdr','??', ...
'dBrMSEdr','T/m', ...
'dBzMSEdr','T/m', ...
'dFrdz','N/m', ...
'dFzdz','N/m', ...
'dfsdz','??', ...
'dfcdz','??', ...
'dfvdz','Wb/m', ...
'dfpdrC','Wb/m', ...
'dfpdzC','Wb/m', ...
'dcdz',dcdz_units, ...
'dfldz','??', ...
'dflvdz','??', ...
'dbpdz','??', ...
'dBrMSEdz','T/m', ...
'dBzMSEdz','T/m', ...
'dbrdz','T/m', ...
'bz0','T', ...
'bzav','T', ...
'cphi',cphi_units, ...
'ip0',ip0_units, ...
'ndecay', 'dim-less', ...
'ncrit', 'dim-less', ...
'm_s', 'dim-less');

%save rzrig_data

model_data = struct( ...
'drdis',drdis, ...		%	OK (comparison with imks, iterminal)
'dzdis',dzdis, ...              %       OK
'dFrdis',dFrdis, ...              %       OK
'dFzdis',dFzdis, ...              %       OK
'drdip',drdip, ...              %       OK
'dzdip',dzdip, ...              %       OK
'dFrdip',dFrdip, ...              %       OK
'dfsdip',dfsdip, ...              %       OK
'drdbetap',drdbetap, ...              %       OK
'dzdbetap',dzdbetap, ...              %       OK
'dFrdbetap',dFrdbetap, ...              %       OK
'drdli',drdli, ...              %       OK
'dzdli',dzdli, ...              %       not checked
'dFrdli',dFrdli, ...              %       OK
'dFrdr',dFrdr, ...              %       OK
'dFzdr',dFzdr, ...              %       not checked
'dfsdr',dfsdr, ...              %       OK
'dfcdr',dfcdr, ...              %       OK
'dfvdr',dfvdr, ...              %       OK
'dcdr',dcdr, ...              %       OK
'dfldr',dfldr, ...              %       not checked
'dflvdr',dflvdr, ...              %       not checked
'dbpdr',dbpdr, ...              %       not checked
'dBrMSEdr',dBrMSEdr, ...              %       not checked
'dBzMSEdr',dBzMSEdr, ...              %       not checked
'dFrdz',dFrdz, ...              %       not checked
'dFzdz',dFzdz, ...              %       OK
'dfsdz',dfsdz, ...              %       OK
'dfcdz',dfcdz, ...              %       OK
'dfvdz',dfvdz, ...              %       OK
'dfpdrC',dfpdrC, ...              %       OK
'dfpdzC',dfpdzC, ...              %       OK
'dcdz',dcdz, ...              %       OK
'dfldz',dfldz, ...              %       not checked
'dflvdz',dflvdz, ...              %       not checked
'dbpdz',dbpdz, ...              %       not checked
'dBrMSEdz',dBrMSEdz, ...              %       not checked
'dBzMSEdz',dBzMSEdz, ...              %       not checked
'dbrdz',dbrdz, ...
'bz0',bz0, ...
'bzav',bzav, ...
'cphi',cphi, ...              %       OK
'ip0',ip0, ...              %       OK
'r0',r0, ...              %       OK
'z0',z0, ...              %       OK
'a0',a0, ...              %       OK
'b0',b0, ...              %       OK
'kap0',kap0, ...              %       OK
'mss',mss, ...
'xmatz',xmatz, ...
'xmatr', xmatr, ...
'ndecay', ndecay, ...
'ncrit', ncrit, ...
'm_s', m_s, ...
'units',units, ...
'descriptions',desc);

if(nargout > 3)
   dbg_objs = struct( ...
      'vals',vals, ...
      'vecs',vecs, ...
      'mpc_r',mpc_r, ...
      'dFrdr_1',dFrdr_1, ...
      'dFrdr_2',dFrdr_2, ...
      'dFrdr_3',dFrdr_3, ...
      'dbzdr',dbzdr);
   if(exist('mpv_r'))
      dbg_objs.mpv_r = mpv_r;
   end
end

%save rzrig_data_new
