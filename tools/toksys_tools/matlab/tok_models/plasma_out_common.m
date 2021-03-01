 %
%  SYNTAX:  plasma_out_common 
%
%  PURPOSE:  common code for all "plasma_output" routines
%
%  INPUT:
%	equil_data = 
%	rzrig_data = 
%	irzresp = set to:
%		  0 = use no rigid response
%		  1 = use only rigid r response
%		  2 = use both rigid r and z response
%	num_Ecoils
%	rgg
%	zg
%	rg
%	mvv
%	jphi0
%	fcnturn
%	cc0 = coil currents from EFIT (MA-turns)
%	include_Ip 
%	EB_out (optional, default=1 if num_Ecoils=5)
%       plasma_obj_file = name of plasma_objects.mat file (if include_Ip=1)
%
%  OUTPUT: (units depend on input objects - need to correct comments)
%   	cvnturn 
%   	drdi 	(m/MA)
%   	dzdi 	(m/MA)
%   	CRZmat 	(m/A)
%   	dcrzdi 	(A/A)
%   	dcrzdip (A/A)
%   	dcrzdbetap (A/?)
%   	dcrzdli  (A/?)
%	ir
%	iz
%	idxpf 
%	cphi0 = EFIT equilibrium plasma current distribution (MA)
%	cc0 = cc0 converted to MA from MA-turns
%   	mu0 = 0.4*pi;
%   	twopi = 2*pi;
%   	twopir = twopi*rgg(:)
%   	dz = (m)
%   	dr = (m)
%   	nvv = size(mvv,1)
%       ncurrents
%       irscl
%       izscl
%	ccnturn
 
%  RESTRICTIONS:
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	9/10/98
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)plasma_out_common.m	1.3 09/01/09

   switch irzresp

      case 0
         irscl=0; izscl=0;
      case 1
         irscl=1; izscl=0;
      case 2
         irscl=0; izscl=1;
      case 3
         irscl=1; izscl=1;
      otherwise
         wait('ERROR plasma_out_common: invalid irzresp value\n');
         return;
   end

   ncurrents=vacuum_objs.nc+vacuum_objs.nv;
   nvv = vacuum_objs.nv;

   mu0 = 0.4*pi;
   twopi = 2*pi;
   twopir = twopi*vacuum_objs.rgg(:);
   dz = vacuum_objs.zg(2)-vacuum_objs.zg(1);
   dr = vacuum_objs.rg(2)-vacuum_objs.rg(1);

   drdi = rzrig_data.drdis;		% m/MA-turn
   dzdi = rzrig_data.dzdis;		% m/MA-turn

   cvnturn = [vacuum_objs.ecnturn; vacuum_objs.fcnturn; ones(nvv,1)];
   nxx = length(cvnturn);
   if length(drdi)~=nxx 	%handle old objects (before RDP struct add 8/99)
      disp(['WARNING: OLD d(*)di being used... adding zeros to be ', ...
   						'compatible with nvv!!!'])
      drdi = [drdi zeros(1,nxx-length(drdi))];
      dzdi = [dzdi zeros(1,nxx-length(dzdi))];
      wait
   end

% EFIT equilibrium plasma current distribution (MA)
   cphi0 = equil_data.jphi*dz*dr; 

   equil_I = cc_efit_to_tok(vacuum_objs,equil_data);
   cc0 = equil_I.cc0e;

% Create RZ Response:
   drdi = irscl*drdi;                                         % m/A
   dzdi = izscl*dzdi;                                         % m/A
   CRZmat = [drdi; dzdi];                                     % m/A

   dcdr = rzrig_data.dcdr;
   dcdz = rzrig_data.dcdz;
   drdbetap = rzrig_data.drdbetap;
   dzdbetap = rzrig_data.dzdbetap;
   drdli = rzrig_data.drdli;
% Create Plasma Current motion response for selected RZ response:
   dcrzdi = dcdr(:)*drdi + dcdz(:)*dzdi; 		% (MA/m)*(m/MA) = A/A 

dcrzdbetap = (dcdr(:)*drdbetap+dcdz(:)*dzdbetap); % A/m*m/beta? = A/beta?
dcrzdli    = dcdr(:)*drdli; 			% A/m * m/li? = A/li?

dcrzdip = irscl*dcdr(:)*rzrig_data.drdip;		% A/A

