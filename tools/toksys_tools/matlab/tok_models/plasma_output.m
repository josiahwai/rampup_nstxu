function output_objs = plasma_output(build_inputs, equil_data, ...
		vacuum_objs,isoflux_defns, rzrig_data)
 %
%  SYNTAX:  plasma_output
%
%  PURPOSE: Calculate DIII-D plasma response output equation objects for
%	rigid plasma.  Outputs are control point fluxes and field at X-point.
%	Run plasma_dynamics.m first to create the file plasma_objects.mat.
%
%  INPUTS:
%   build_inputs = structure containing:
%             tokamak: name of the device (e.g. 'd3d', 'nstx', 'east', etc)
%            ichooseq: 
%      output_irzresp: 
%                       0 = use no rigid response
%                       1 = use only rigid r response
%                       2 = use only rigid z response
%                       3 = use both rigid r and z response
%        vacuum_model: 
%             iplcirc: set=1 to include Ip variation model (optional, default 0)
%              Te_res: 
%              li_res: 
%                  Rp: 
%            Zeff_res: 
%             ccnturn: 
% efit_gfile OR efit_source
%                fcid: 
%                vvid: 
%              vvfrac: 
%        calc_resp: 
%       make_resp_dir: 
%         efit_source: 
%             shotnum: 
%                time: 
%          num_Ecoils: for D3D, set to 2 or 6; all other devices, 0 (default)
%   vacuum_objs   =
%   isoflux_defns = structure containing isonms, ref_num, bgrdnms
%   rzrig_data    = 
%   scale_cc_resp = (optional)
%   scale_vv_resp = (optional)

%  ARE THESE USED??
%   From plasma_objects.mat (only needed if iplcirc=1 and Rp=0):
%	Mpd = mutuals from plasma filament to conductors (H)
%	Lp = plasma self-inductance(H) 
%	Rp = plasma resistance(Ohms)
%	Xpp = (H)
%	Xpd = (H)
%   plasma_obj_file = name of plasma_objects.mat file (if needed, see above)
%   plasma_obj_log = Text file containing parameter settings which created
%	plasma_objects file.  Must have same prefix as plasma_objects file
%	and suffix = ".log" (file name will be constructed if needed)
%
%     efit_gfile = EFIT g0-file name
%	OR
%	efit_source
%
%  OUTPUTS:
%    output_objs = structure (saved in  out_objs_<shot>_<time>.mat) containing:
%     nptiso = array of number of points in each isoflux segment (isonms)
%     nptbgd = array of number of points in each bgrid segment (bgrdnms)
%              (define here totiso=sum(nptiso), totbgd=sum(nptbgd),
%		   ncond = necoils + nfcoils + nvv)
%     nptbgd2 = array of number of points in each bgrid segment (bgrdnms2)
%     CXmatI = Ctrl. pt. and gridpt flux response due to plasma 
%	    	motion from I_cond (=d[Psi_isflux, Psi_grid]pl/dIcond), 
%	    	Wb/A  (size = (totiso+totbgd) x ncond)
%     CXmatBR, CXmatBZ = gridpt. Br, Bz response due to plasma 
%	    	motion from I_cond (=d[Br_grid, Bz_grid,...]pl/dIcond), 
%	    	T/A (size = totbgd x ncond)
%     CmatI  = Ctrl. pt. and gridpt flux, response due to plasma+cond from 
%		I_cond (=d[Psi_segments,Psi_grid(s)]tot/dIcond), Wb/A 
%		(size = (totiso+totbgd) x ncond). For indices representing grid,
%		first index represents indexing through the grid by rows first.
%     CmatBR, CmatBZ = Grid pt. Br, Bz response due to plasma+cond from 
%		I_cond (=d[Br, Bz,...]tot/dIcond), T/A (size = totbgd x ncond)
%		First index represents indexing through the grid by rows first.
%     CRZmat = Plasma R, Z response to I_cond (m/A)
%     Psi0   = equilibrium values of flux (Wb) on all control pts
%		(size = totiso+totbgd x 1)
%     BR0, BZ0 = equilibrium values of field (T) on control grid
%		(size = totbgd x 1)
%     CYmatI = Ctrl. pt. and gridpt flux response due to plasma 
%	    	current change (=d[Psi_isflux,Psi_grid]/dIpl)
%		Wb/A (size = totiso+totbgd x 1)
%     CYmatBR, CYmatBZ = grid pt. Br, Bz response due to plasma 
%	    	current change (=d[Brgrid, Bzgrid]pl/dIpl)
%		T/A (size = totbgd x ncond)
%     CmatIbetap = ctrl pt and gridpt flux response 
%		due to changes in betap (=d[Psi_iso, Psi_grid]pl/dbetap), 
%		Wb/unit-beta (size = totiso+totbgd x 1)
%     CmatBRbetap, CmatBZbetap = gridpt Br, Bz response 
%		due to changes in betap (=d[Bgrid]pl/dbetap), T/unit-beta
%		(size = totbgd x 1)
%     CmatIli = ctrl pt and gridpt flux response due to 
%		changes in li (=d[Psi_iso,Psi_grid  ]pl/dli), Wb/unit-li 
%		(size = totiso+totbgd x 1)
%     CmatBRli, CmatBZli = gridpt Br, Bz response due to 
%		changes in li (=d[Bgrid ]pl/dli), T/unit-li
%		(size = totbgd x 1)
%     rbgd = r coordinates of grid pts
%     zbgd = z coordinates of grid pts
%     rbgd2 = r coordinates of grid pts, grid 2
%     zbgd2 = z coordinates of grid pts, grid 2
%  Log file out_objs_<shot>_<time>.log
%
%  RESTRICTIONS:
%   (1) Works only for D3D and NSTX (so far).
%   (2) Ip model only works for single filament (idxfil)
%   (3) Data file data_file must exist and contain Green function data for
%     control point response to current in conductors and on grid.
%   (4) NOTE that sizes of output objects are defined by the data in the 
%	"sizedata" array - may not match size in tok_data_struct OR size in 
%	tok_system.
%
%  METHOD:
%  (1) Green's functions from coils, vessel, and plasma grid points to
%	the isoflux control points are stored in files and read in by
%	read_response.m - "help read_response" for Green function conventions.
%  (2) plasma current variation objects only have 1 d.o.f. but a distributed
%      current; a current wtd average is used.

%  WRITTEN BY:  Mike Walker 	ON	10/24/97 from Dave's make_rig_plresp.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)plasma_output.m	1.4 09/02/09

% NOTE: Indices of fcoils, ecoils, and vessel are FIXED in the PCS files. 
% They do not change when we change number of coils being modeled.

  struct_to_ws(build_inputs);
  struct_to_ws(isoflux_defns);

  shotstr = int2str(build_inputs.shotnum);
  timestr = int2str(build_inputs.time);
  out_file_name = ['out_objs_' shotstr '_' timestr]

clear GRbgd_pf

  niso = size(isonms,1);
  if(~isempty(bgrdnms))
     nbgrd = 1;
  end
  if(size(bgrdnms,1)==2)
     nbgrd2 = 1;
  end

  out_log_file = [out_file_name '.log']
  if(iplcirc & Rp==0)
     if(~exist('plasma_obj_file'))
        disp(['plasma_output ERROR: '...
		' plasma_obj_file must be provided when iplcirc=1 and Rp=0'])
        return
     end
     plasma_obj_log = [plasma_obj_file '.log']
     err = unix(['cat ' plasma_obj_log ' >> ' out_log_file]);
     if(err~=0) wait; end;
  else
     fid_log = fopen(out_log_file,'a+');
     fprintf(fid_log,'\n')
     fprintf(fid_log,'No plasma object file data provided.\n')
     fprintf(fid_log,'\n')
     fclose(fid_log)
  end
  
  irzresp = output_irzresp;

  fid_log = fopen(out_log_file,'a+');
  fprintf(fid_log,'Inputs to plasma_output,  %s \n',date);
%  fprintf(fid_log,'  num_Ecoils = %d\n',num_Ecoils);
  fprintf(fid_log,'  output_irzresp = %d\n',output_irzresp);
  fprintf(fid_log,'  fcnturn = ');
  fprintf(fid_log,' %d',vacuum_objs.fcnturn); fprintf(fid_log,'\n');
  fprintf(fid_log,'  nvv	= %d\n',vacuum_objs.nv);
  if(iplcirc & Rp==0)
     fprintf(fid_log,'  plasma_obj_file = %s\n', plasma_obj_file);
     fprintf(fid_log,'  plasma_obj_log = %s\n', plasma_obj_log);
  end
  if(exist('efit_gfile'))
     fprintf(fid_log,'  efit_gfile = %s\n',efit_gfile);
  else
%     fprintf(fid_log,'  efit_source = %s\n',efit_source);
  end
  fprintf(fid_log,'  iplcirc = %d\n',iplcirc);

  fprintf(fid_log,'  isonms = \n');
  for k=1:size(isonms,1) fprintf(fid_log,'      %s \n',isonms(k,:)); end;
  fprintf(fid_log,'\n');

  fprintf(fid_log,'  bgrdnms = \n');
  for k=1:size(bgrdnms,1) fprintf(fid_log,'      %s \n',bgrdnms(k,:)); end;
  fprintf(fid_log,'\n');

  fprintf(fid_log,'  niso = %d\n',niso);
  fprintf(fid_log,'  nbgrd = %d\n',nbgrd);
  fprintf(fid_log,'  nbgrd2 = %d\n',nbgrd2);
  fprintf(fid_log,'  out_file_name =  %s\n', out_file_name);
  fprintf(fid_log,'  out_log_file  =  %s\n', out_log_file);
  fclose(fid_log)

  plasma_out_common		% common code for all plasma_out routines

% Do initial plot (on which to plot all segments):
   figure(1),clf
   plot_tok_geo(vacuum_objs,[],equil_data);
   hold on

% need to create some room for very large data objects below

clear mpc mpc_r mpc_z mpc_zz mpl mpv mpv_r gpb gbv gbc rgg zgg tmp mlc mlv mvv
clear jphi jphi0 itmp ecdata cphi br br_r bz bz_z

% Need to do a little work to get dcrzdi to match efit currents, since
% number of coils or vessel elts in efit might not match our objects.

% Need to get sizedata first:
if(calc_resp)
   idx = findstr('segment',isonms(1,:));
   seg_grid_name = isonms(1,1:idx(1)-2);
   [sizedata,control_pts,mutuals] = ...
          calc_isoflux_response(seg_grid_name,make_resp_dir,vacuum_objs);
else
   if(ref_num(1)==0)
      filename = 'unused'
   else
      filename = [remove_space(isonms(1,:))]
   end
   [sizedata,control_pts,mutuals] = read_response(filename,gridsize);
end

% replace "<" with "~=" here when d3d vacuum objects uses 6 E coils:
if(sizedata.nfcoils+sizedata.necoils<vacuum_objs.nc & exist('fcid')~=1)
   wait('ERROR plasma_output: incompatible fcoil vacuum and equil objects')
   return;
end
if(sizedata.nvessel~=vacuum_objs.nv & exist('vvid')~=1)
   wait('ERROR plasma_output: incompatible vessel vacuum and equil objects')
   return;
end

ncoils = sizedata.nfcoils + sizedata.necoils;
if(num_Ecoils==2)	% special case for D3D
   ncoils = sizedata.nfcoils + 2;
end
if(exist('fcid')==1)
   cproj = zeros(vacuum_objs.nc,ncoils+sizedata.nvessel);
   for k=1:sizedata.necoils
      cproj(k,k)=1;
   end
   for k=1:vacuum_objs.nc-sizedata.necoils
      cproj(k+sizedata.necoils,fcid(k)+sizedata.necoils)=1;
   end
else
   cproj = eye(vacuum_objs.nc,ncoils+sizedata.nvessel); % default
end

if(exist('vvid')==1)
   if(exist('vvfrac')~=1)
       vvfrac = ones(size(vvid));
   end
   vproj = zeros(vacuum_objs.nv,sizedata.nvessel);
   for k=1:vacuum_objs.nv
      vproj(k,vvid(k))=vvfrac(k);
   end
   vproj = [zeros(vacuum_objs.nv,ncoils) vproj];
else
   vproj = [zeros(vacuum_objs.nv,ncoils) eye(vacuum_objs.nv,sizedata.nvessel)];
end
proj = [cproj; vproj];
dcrzdi = dcrzdi*proj;

% Do all the isoflux segment CXmat parts:
 CXmatI=[]; CmatI=[]; Psi0=[]; riso=[]; ziso=[]; nptiso=[]; Psi0p=[]; Psi0c=[];
 CYmatI=[]; CmatIbetap=[]; CmatIli=[];

 for ii=1:niso
   if(calc_resp)
      if(ref_num(ii)==0)
         sizedata.ncntrlpts = 0;
      else
         idx = findstr('segment',isonms(ii,:));
         seg_grid_name = isonms(ii,1:idx(1)-2);
         [sizedata,control_pts,mutuals] = ...
          calc_isoflux_response(seg_grid_name,make_resp_dir,vacuum_objs);
      end
   else
      if(ref_num(ii)==0)
         filename = 'unused'
      else
         filename = [remove_space(isonms(ii,:))]
      end
      [sizedata,control_pts,mutuals] = read_response(filename,gridsize);
   end

   if(sizedata.ncntrlpts > 0)
      idxf = 1:sizedata.nfcoils; 
      idxe = sizedata.nfcoils+[1:sizedata.necoils];
      idxv = sizedata.nfcoils+sizedata.necoils+[1:sizedata.nvessel];
      idxp = sizedata.nfcoils+sizedata.necoils+sizedata.nvessel ...
		+ [1:sizedata.ngrid];

% This is to handle cases where non-standard units are used in computing the
% mutuals and Greens functions:
      if(exist('scale_cc_resp')==1)
         mutuals([idxe idxf],:) = diag(scale_cc_resp)*mutuals([idxe idxf],:);
      end
      if(exist('scale_vv_resp')==1)
         mutuals(idxv,:) = diag(scale_vv_resp)*mutuals(idxv,:);
      end

      rcntrlpt = control_pts(1,:)';
      zcntrlpt = control_pts(2,:)'; 
      Miso_c = [mutuals(idxe,:)' mutuals(idxf,:)']*twopi;	% Wb/A-turn
      if num_Ecoils==2
         Miso_c(:,1) = sum(Miso_c(:,[1 3 5])')' / ccnturn(1);
         Miso_c(:,2) = sum(Miso_c(:,[2 4 6])')' / ccnturn(2);
         Miso_c = Miso_c(:,[1 2 7:24]);
      end
      Miso_c = Miso_c.*(ones(size(Miso_c,1),1)*ccnturn');		% Wb/A
      Miso_v = mutuals(idxv,:)'*twopi;					% Wb/A
      Miso_p = mutuals(idxp,:)'*twopi;  				% Wb/A
      Miso_pf= (cphi0(:)'*mutuals(idxp,:))'/sum(cphi0(:))*twopi;	% Wb/A

      CXmatI1 = Miso_p*dcrzdi;                            		% Wb/A
      CYmatI1   = Miso_p*dcrzdip;
      CmatI1 = CXmatI1 + [Miso_c Miso_v];			% Wb/A
      CXmatI = [CXmatI; CXmatI1];					% Wb/A

      CYmatI   = [CYmatI; CYmatI1+Miso_pf];
      CmatI = [CmatI; CmatI1];                                  	% Wb/A
      CmatIbetap = [CmatIbetap; Miso_p*dcrzdbetap];
      CmatIli = [CmatIli; Miso_p*dcrzdli];

% Note that cc0 is already converted to MA by plasma_out_common
      Psi0 = [Psi0; Miso_p*cphi0(:)*1e6+Miso_c*cc0];		% Wb 
      Psi0p = [Psi0p; Miso_p*cphi0(:)*1e6];
      Psi0c = [Psi0c; Miso_c*cc0];
      riso=[riso; rcntrlpt];
      ziso=[ziso; zcntrlpt];
      nptiso = [nptiso; length(rcntrlpt)]

      plot(rcntrlpt,zcntrlpt,'m.'); pause(1);
   else
      nptiso = [nptiso; 0]
   end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Xpoint grid 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do all the bgrid rows CXmat parts:
 CXmatBZ=[]; CmatBZ=[]; BZ0=[];
 CXmatBR=[]; CmatBR=[]; BR0=[];
 CYmatBR=[]; CYmatBZ=[]; CmatBZbetap=[]; CmatBRbetap=[];
 CmatBZli=[]; CmatBRli=[];
 rbgd=[]; zbgd=[]; nptbgd=[]; Mbgd_c0=[];
 for ii=1:nbgrd
   if(calc_resp)
      idx = findstr('grid',bgrdnms(ii,:));
      seg_grid_name = bgrdnms(ii,1:idx(1)-2);
      [sizedata,control_pts,mutuals,brgreens,bzgreens] = ...
          calc_isoflux_response(seg_grid_name,make_resp_dir,vacuum_objs);
   else
      filename = [remove_space(bgrdnms(ii,:))]
      [sizedata,control_pts,mutuals,brgreens,bzgreens] = ...
			read_response(filename,gridsize);
   end

   idxf = 1:sizedata.nfcoils; 
   idxe = sizedata.nfcoils+[1:sizedata.necoils];
   idxv = sizedata.nfcoils+sizedata.necoils+[1:sizedata.nvessel];
   idxp = sizedata.nfcoils+sizedata.necoils+sizedata.nvessel ...
		+ [1:sizedata.ngrid];

% This is to handle cases where non-standard units are used in computing the
% mutuals and Greens functions:
   if(exist('scale_cc_resp')==1)
      mutuals([idxe idxf],:) = diag(scale_cc_resp)*mutuals([idxe idxf],:);
      brgreens([idxe idxf],:) = diag(scale_cc_resp)*brgreens([idxe idxf],:);
      bzgreens([idxe idxf],:) = diag(scale_cc_resp)*bzgreens([idxe idxf],:);
   end
   if(exist('scale_vv_resp')==1)
      mutuals(idxv,:) = diag(scale_vv_resp)*mutuals(idxv,:);
      brgreens(idxv,:) = diag(scale_vv_resp)*brgreens(idxv,:);
      bzgreens(idxv,:) = diag(scale_vv_resp)*bzgreens(idxv,:);
   end

   rcntrlpt = control_pts(1,:);
   zcntrlpt = control_pts(2,:);

   GRbgd_c = [brgreens(idxe,:)' brgreens(idxf,:)'];		% T/A-turn
   if num_Ecoils==2
      GRbgd_c(:,1) = sum(GRbgd_c(:,[1 3 5])')' / ccnturn(1);
      GRbgd_c(:,2) = sum(GRbgd_c(:,[2 4 6])')' / ccnturn(2);
      GRbgd_c = GRbgd_c(:,[1 2 7:24]);
   end
   GRbgd_c = GRbgd_c.*(ones(size(GRbgd_c,1),1)*ccnturn');       % T/A
   GRbgd_v = brgreens(idxv,:)';					% T/A
   GRbgd_p = brgreens(idxp,:)';					% T/A
%->
   GRbgd_pf = (cphi0(:)'*brgreens(idxp,:))'/sum(cphi0(:));	% T/A
   clear brgreens
   CXmatBR1 = GRbgd_p*dcrzdi;                                   % T/A
   CXmatBR = [CXmatBR; CXmatBR1];                           	% T/A
   CYmatBR1 = GRbgd_p*dcrzdip;
   CYmatBR = [CYmatBR; CYmatBR1+GRbgd_pf];
   CmatBR1 = CXmatBR1 + [GRbgd_c GRbgd_v];           % T/A
   CmatBR = [CmatBR; CmatBR1];                           	% T/A
   CmatBRbetap = [CmatBRbetap; GRbgd_p*dcrzdbetap];
   CmatBRli = [CmatBRli; GRbgd_p*dcrzdli];
%<-
% Note that cc0 is already converted to MA by plasma_out_common
   BR0 = [BR0; GRbgd_p*cphi0(:)*1e6 + GRbgd_c*cc0];% T

   GZbgd_c = [bzgreens(idxe,:)' bzgreens(idxf,:)'];             %T/A-turn
   if num_Ecoils==2
      GZbgd_c(:,1) = sum(GZbgd_c(:,[1 3 5])')' / ccnturn(1);
      GZbgd_c(:,2) = sum(GZbgd_c(:,[2 4 6])')' / ccnturn(2);
      GZbgd_c = GZbgd_c(:,[1 2 7:24]);
   end
   GZbgd_c = GZbgd_c.*(ones(size(GZbgd_c,1),1)*ccnturn');       %T/A
   GZbgd_v = bzgreens(idxv,:)';
   GZbgd_p = bzgreens(idxp,:)';
%->
   GZbgd_pf = (cphi0(:)'*bzgreens(idxp,:))'/sum(cphi0(:));	% T/A
   clear bzgreens
   CXmatBZ1 = GZbgd_p*dcrzdi;                                   % T/A
   CXmatBZ = [CXmatBZ; CXmatBZ1];				% T/A
   CYmatBZ1 = GZbgd_p*dcrzdip;
   CYmatBZ = [CYmatBZ; CYmatBZ1+GZbgd_pf];
   CmatBZ1 = CXmatBZ1 + [GZbgd_c GZbgd_v];
   CmatBZ = [CmatBZ; CmatBZ1];					% T/A
   CmatBZbetap = [CmatBZbetap; GZbgd_p*dcrzdbetap];
   CmatBZli = [CmatBZli; GZbgd_p*dcrzdli];
%<-
% Note that cc0 is already converted to MA by plasma_out_common
   BZ0 = [BZ0; GZbgd_p*cphi0(:)*1e6 + GZbgd_c*cc0];% T

   Mbgd_c = [mutuals(idxe,:)' mutuals(idxf,:)']*twopi;		%Wb/A-turn

   if num_Ecoils==2
% For D3D, the values for Mbgd_c that are stored in PCS assume the E-coil
% is a 1-turn object (and therefore assumes #turns multiplies the data).
% In plasma_out_common.m, we generate E-coil currents that are TERMINAL
% current values, and therefore Mbgd_c needs to be scaled to account for that.
% ARGUMENT IS WRONG - WHAT IS CORRECT ARGUMENT???
      Mbgd_c(:,1) = sum(Mbgd_c(:,[1 3 5])')' / ccnturn(1);
      Mbgd_c(:,2) = sum(Mbgd_c(:,[2 4 6])')' / ccnturn(2);
      Mbgd_c = Mbgd_c(:,[1 2 7:24]);
   end
% QUESTION IS: HOW DO NSTX OBJECTS HANDLE ECOIL?

   Mbgd_c = Mbgd_c.*(ones(size(Mbgd_c,1),1)*ccnturn');		% Wb/A
   Mbgd_c0 = [Mbgd_c0; Mbgd_c];
   Mbgd_v = mutuals(idxv,:)'*twopi;				%Wb/A
   Mbgd_p = mutuals(idxp,:)'*twopi;				%Wb/A
%->
   Mbgd_pf= (cphi0(:)'*mutuals(idxp,:))'/sum(cphi0(:))*twopi;	% Wb/A
   clear mutuals
   CXmatI1 = Mbgd_p*dcrzdi;                                     % Wb/A
   CXmatI = [CXmatI; CXmatI1];					% Wb/A
   CYmatI1 = Mbgd_p*dcrzdip;
   CYmatI  = [CYmatI; CYmatI1+Mbgd_pf];
   CmatI1 = CXmatI1 + [Mbgd_c Mbgd_v];
   CmatI = [CmatI; CmatI1];					% Wb/A
   CmatIbetap = [CmatIbetap; Mbgd_p*dcrzdbetap];
   CmatIli = [CmatIli; Mbgd_p*dcrzdli];
%<-
% Note that cc0 is already converted to MA by plasma_out_common
   Psi0 = [Psi0; Mbgd_p*cphi0(:)*1e6 + Mbgd_c*cc0];% Wb 
   Psi0p = [Psi0p; Mbgd_p*cphi0(:)*1e6];
   Psi0c = [Psi0c; Mbgd_c*cc0];
   rbgd= rcntrlpt;
   zbgd= zcntrlpt;
   nptbgd = length(rcntrlpt);

   plot(rcntrlpt,zcntrlpt,'c.'); pause(1);
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Xpoint grid 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do all the bgrid rows CXmat parts:
 CXmatBZ2=[]; CmatBZ2=[]; BZ20=[];
 CXmatBR2=[]; CmatBR2=[]; BR20=[];
 CYmatBR2=[]; CYmatBZ2=[]; CmatBZ2betap=[]; CmatBR2betap=[];
 CmatBZ2li=[]; CmatBR2li=[];
 rbgd2=[]; zbgd2=[]; nptbgd2=[]; Mbgd2_c0=[];
 for ii=1:nbgrd2
   if(calc_resp)
      idx = findstr('grid',bgrdnms(2,:));
      seg_grid_name = bgrdnms(2,1:idx(1)-2);
      [sizedata,control_pts,mutuals,brgreens,bzgreens] = ...
          calc_isoflux_response(seg_grid_name,make_resp_dir,vacuum_objs);
   else
      filename = [remove_space(bgrdnms(2,:))]
      [sizedata,control_pts,mutuals,brgreens,bzgreens] = ...
			read_response(filename,gridsize);
   end

   idxf = 1:sizedata.nfcoils; 
   idxe = sizedata.nfcoils+[1:sizedata.necoils];
   idxv = sizedata.nfcoils+sizedata.necoils+[1:sizedata.nvessel];
   idxp = sizedata.nfcoils+sizedata.necoils+sizedata.nvessel ...
		+ [1:sizedata.ngrid];

% This is to handle cases where non-standard units are used in computing the
% mutuals and Greens functions:
   if(exist('scale_cc_resp')==1)
      mutuals([idxe idxf],:) = diag(scale_cc_resp)*mutuals([idxe idxf],:);
      brgreens([idxe idxf],:) = diag(scale_cc_resp)*brgreens([idxe idxf],:);
      bzgreens([idxe idxf],:) = diag(scale_cc_resp)*bzgreens([idxe idxf],:);
   end
   if(exist('scale_vv_resp')==1)
      mutuals(idxv,:) = diag(scale_vv_resp)*mutuals(idxv,:);
      brgreens(idxv,:) = diag(scale_vv_resp)*brgreens(idxv,:);
      bzgreens(idxv,:) = diag(scale_vv_resp)*bzgreens(idxv,:);
   end

   rcntrlpt = control_pts(1,:);
   zcntrlpt = control_pts(2,:);

   GRbgd2_c = [brgreens(idxe,:)' brgreens(idxf,:)'];		% T/A-turn
   if num_Ecoils==2
      GRbgd2_c(:,1) = sum(GRbgd2_c(:,[1 3 5])')' / ccnturn(1);
      GRbgd2_c(:,2) = sum(GRbgd2_c(:,[2 4 6])')' / ccnturn(2);
      GRbgd2_c = GRbgd2_c(:,[1 2 7:24]);
   end
   GRbgd2_c = GRbgd2_c.*(ones(size(GRbgd2_c,1),1)*ccnturn');       % T/A
   GRbgd2_v = brgreens(idxv,:)';					% T/A
   GRbgd2_p = brgreens(idxp,:)';					% T/A
%->
   GRbgd2_pf = (cphi0(:)'*brgreens(idxp,:))'/sum(cphi0(:));	% T/A
   clear brgreens
   CXmatBR21 = GRbgd2_p*dcrzdi;                                   % T/A
   CXmatBR2 = [CXmatBR2; CXmatBR21];                           	% T/A
   CYmatBR21 = GRbgd2_p*dcrzdip;
   CYmatBR2 = [CYmatBR2; CYmatBR21+GRbgd2_pf];
   CmatBR21 = CXmatBR21 + [GRbgd2_c GRbgd2_v];           % T/A
   CmatBR2 = [CmatBR2; CmatBR21];                           	% T/A
   CmatBR2betap = [CmatBR2betap; GRbgd2_p*dcrzdbetap];
   CmatBR2li = [CmatBR2li; GRbgd2_p*dcrzdli];
%<-
% Note that cc0 is already converted to MA by plasma_out_common
   BR20 = [BR20; GRbgd2_p*cphi0(:)*1e6 + GRbgd2_c*cc0];% T

   GZbgd2_c = [bzgreens(idxe,:)' bzgreens(idxf,:)'];             %T/A-turn
   if num_Ecoils==2
      GZbgd2_c(:,1) = sum(GZbgd2_c(:,[1 3 5])')' / ccnturn(1);
      GZbgd2_c(:,2) = sum(GZbgd2_c(:,[2 4 6])')' / ccnturn(2);
      GZbgd2_c = GZbgd2_c(:,[1 2 7:24]);
   end
   GZbgd2_c = GZbgd2_c.*(ones(size(GZbgd2_c,1),1)*ccnturn');       %T/A
   GZbgd2_v = bzgreens(idxv,:)';
   GZbgd2_p = bzgreens(idxp,:)';
%->
   GZbgd2_pf = (cphi0(:)'*bzgreens(idxp,:))'/sum(cphi0(:));	% T/A
   clear bzgreens
   CXmatBZ21 = GZbgd2_p*dcrzdi;                                   % T/A
   CXmatBZ2 = [CXmatBZ2; CXmatBZ21];				% T/A
   CYmatBZ21 = GZbgd2_p*dcrzdip;
   CYmatBZ2 = [CYmatBZ2; CYmatBZ21+GZbgd2_pf];
   CmatBZ21 = CXmatBZ21 + [GZbgd2_c GZbgd2_v];
   CmatBZ2 = [CmatBZ2; CmatBZ21];					% T/A
   CmatBZ2betap = [CmatBZ2betap; GZbgd2_p*dcrzdbetap];
   CmatBZ2li = [CmatBZ2li; GZbgd2_p*dcrzdli];
%<-
% Note that cc0 is already converted to MA by plasma_out_common
   BZ20 = [BZ20; GZbgd2_p*cphi0(:)*1e6 + GZbgd2_c*cc0];% T

   Mbgd2_c = [mutuals(idxe,:)' mutuals(idxf,:)']*twopi;		%Wb/A-turn
   if num_Ecoils==2
% ?????????????????????????????????????????????????????????????????????
% D3D-specific handling of turns ?????
% ?????????????????????????????????????????????????????????????????????
      Mbgd2_c(:,1) = sum(Mbgd2_c(:,[1 3 5])')' / ccnturn(1);
      Mbgd2_c(:,2) = sum(Mbgd2_c(:,[2 4 6])')' / ccnturn(2);
      Mbgd2_c = Mbgd2_c(:,[1 2 7:24]);
   end
   Mbgd2_c = Mbgd2_c.*(ones(size(Mbgd2_c,1),1)*ccnturn');		% Wb/A
   Mbgd2_c0 = [Mbgd2_c0; Mbgd2_c];
   Mbgd2_v = mutuals(idxv,:)'*twopi;				%Wb/A
   Mbgd2_p = mutuals(idxp,:)'*twopi;				%Wb/A
%->
   Mbgd2_pf= (cphi0(:)'*mutuals(idxp,:))'/sum(cphi0(:))*twopi;	% Wb/A
   clear mutuals
   CXmatI1 = Mbgd2_p*dcrzdi;                                     % Wb/A
   CXmatI = [CXmatI; CXmatI1];					% Wb/A
   CYmatI1 = Mbgd2_p*dcrzdip;
   CYmatI  = [CYmatI; CYmatI1+Mbgd2_pf];
   CmatI1 = CXmatI1 + [Mbgd2_c Mbgd2_v];
   CmatI = [CmatI; CmatI1];					% Wb/A
   CmatIbetap = [CmatIbetap; Mbgd2_p*dcrzdbetap];
   CmatIli = [CmatIli; Mbgd2_p*dcrzdli];
%<-
% Note that cc0 is already converted to MA by plasma_out_common
   Psi0 = [Psi0; Mbgd2_p*cphi0(:)*1e6 + Mbgd2_c*cc0];% Wb 
   Psi0p = [Psi0p; Mbgd2_p*cphi0(:)*1e6];
   Psi0c = [Psi0c; Mbgd2_c*cc0];
   rbgd2=[rbgd2; rcntrlpt];
   zbgd2=[zbgd2; zcntrlpt];
   nptbgd2 = [nptbgd2; length(rcntrlpt)]

   plot(rcntrlpt,zcntrlpt,'c.'); pause(1);
 end

   
if iplcirc
%  If 0 plasma resistance, Ip relation becomes algebraic and outputs come from
%  coil and vessel currents "through" plasma.
   if Rp==0
      CYmatI  = CYmatI * (-(Mpd+Xpd)/(Lp+Xpp));
      CYmatBR = CYmatBR* (-(Mpd+Xpd)/(Lp+Xpp));
      CYmatBZ = CYmatBZ* (-(Mpd+Xpd)/(Lp+Xpp));
   end
end

if(nptbgd2~=0) 
   ngrids=2;
else
   ngrids=1;
end
output_objs = struct( ...
'CRZmat',CRZmat, ...
'nsegments',niso, ...
'ngrids',ngrids, ...
'nptiso',nptiso, ...
'riso',riso, ...
'ziso',ziso, ...
'cc0',cc0, ...
'ccnturn',ccnturn, ...
'Psi0',Psi0 , ...  
'CmatI',CmatI, ...
'CXmatI',  CXmatI, ...
'CYmatI',CYmatI, ...
'CmatIbetap',CmatIbetap, ...
'CmatIli',CmatIli, ...
'nptbgd',nptbgd, ...
'rbgd',rbgd, ...
'zbgd',zbgd, ...
'BR0',BR0, ...
'BZ0',BZ0, ...
'CmatBR',CmatBR, ...
'CmatBZ',CmatBZ, ...
'CXmatBR',CXmatBR, ...
'CXmatBZ',CXmatBZ, ...
'CYmatBR',CYmatBR, ...
'CYmatBZ',CYmatBZ, ...
'CmatBRbetap',CmatBRbetap, ...
'CmatBZbetap',CmatBZbetap, ...
'CmatBRli',CmatBRli, ...
'CmatBZli',CmatBZli, ...
'nptbgd2',nptbgd2, ...
'rbgd2',rbgd2, ...
'zbgd2',zbgd2, ...
'BR20',BR20, ...
'BZ20',BZ20, ...
'CmatBR2',CmatBR2, ...
'CmatBZ2',CmatBZ2, ...
'CXmatBR2',CXmatBR2, ...
'CXmatBZ2',CXmatBZ2, ...
'CYmatBR2',CYmatBR2, ...
'CYmatBZ2',CYmatBZ2, ...
'CmatBR2betap',CmatBR2betap, ...
'CmatBZ2betap',CmatBZ2betap, ...
'CmatBR2li',CmatBR2li, ...
'CmatBZ2li',CmatBZ2li, ...
'dcrzdi', dcrzdi, ...
'dcrzdip', dcrzdip, ...
'dcrzdbetap', dcrzdbetap, ...
'dcrzdli', dcrzdli, ...
'iplcirc',iplcirc);

segment_green_fns = struct( ...
'Miso_c',  Miso_c, ...
'Miso_v', Miso_v, ...
'Miso_p', Miso_p, ...
'Miso_pf', Miso_pf);

grid_green_fns = struct( ...
'GRbgd_c', GRbgd_c, ... 
'GRbgd_v', GRbgd_v, ...
'GRbgd_p', GRbgd_p, ... 
'GRbgd_pf',GRbgd_pf, ...
'GZbgd_c', GZbgd_c, ... 
'GZbgd_v', GZbgd_v, ... 
'GZbgd_p', GZbgd_p, ...
'GZbgd_pf', GZbgd_pf, ...
'Mbgd_c', Mbgd_c, ... 
'Mbgd_v', Mbgd_v, ... 
'Mbgd_p', Mbgd_p, ...
'Mbgd_pf', Mbgd_pf, ...
'GRbgd2_c', GRbgd2_c, ... 
'GRbgd2_v', GRbgd2_v, ...
'GRbgd2_p',GRbgd2_p, ... 
'GRbgd2_pf',GRbgd2_pf, ...
'GZbgd2_c', GZbgd2_c, ... 
'GZbgd2_v', GZbgd2_v, ... 
'GZbgd2_p', GZbgd2_p, ...
'GZbgd2_pf', GZbgd2_pf, ...
'Mbgd2_c', Mbgd2_c, ... 
'Mbgd2_v', Mbgd2_v, ... 
'Mbgd2_p', Mbgd2_p, ...
'Mbgd2_pf', Mbgd2_pf);

output_objs.segment_green_fns = segment_green_fns;
output_objs.grid_green_fns = grid_green_fns;

output_objs.description = struct( ...
'CRZmat','Plasma R, Z response to conductor currents (m/A) ', ...
'nsegments','number of control segments ', ...
'ngrids','number of X point grids ', ...
'nptiso','array of number of points in each isoflux segment (isonms) ', ...
'riso','radial coordinate of isoflux control points, all segments', ...
'ziso','vertical coordinate of isoflux control points, all segments', ...
'Psi0','equilibrium values of flux (Wb) on all control pts (totiso+totbgd x 1) ', ...  
'CmatI','Ctrl. pt. and gridpt flux, total response due to plasma+cond from I_cond (=d[Psi_segments,Psi_grid(s)]tot/dIcond), Wb/A ((totiso+totbgd) x ncond). For indices representing grid, first index represents indexing through the grid by rows first. ', ...
'CXmatI','Ctrl. pt. and gridpt flux response due only to plasma motion from I_cond (=d[Psi_isflux, Psi_grid]pl/dIcond), Wb/A  ((totiso+totbgd) x ncond) ', ...
'CYmatI','Ctrl. pt. and gridpt flux response due to plasma current change (=d[Psi_isflux,Psi_grid]/dIpl) Wb/A (totiso+totbgd x 1) ', ...
'CmatIbetap','ctrl pt and gridpt flux response due to changes in betap (=d[Psi_iso, Psi_grid]pl/dbetap), Wb/unit-beta (totiso+totbgd x 1) ', ...
'CmatIli','ctrl pt and gridpt flux response due to changes in li (=d[Psi_iso,Psi_grid  ]pl/dli), Wb/unit-li (totiso+totbgd x 1) ', ...
'nptbgd','array of number of points in each bgrid segment (bgrdnms) ', ...
'rbgd','r coordinates of grid pts ', ...
'zbgd','z coordinates of grid pts ', ...
'BR0','equilibrium values of field (T) on control grid (totbgd x 1) ', ...
'BZ0','equilibrium values of field (T) on control grid (totbgd x 1) ', ...
'CmatBR','Grid pt. Br total response due to plasma+cond from I_cond (=d[Br, Bz,...]tot/dIcond), T/A (totbgd x ncond) First index represents indexing through the grid by rows first. ', ...
'CmatBZ','Grid pt. Bz total response due to plasma+cond from I_cond (=d[Br, Bz,...]tot/dIcond), T/A (totbgd x ncond) First index represents indexing through the grid by rows first. ', ...
'CXmatBR','gridpt. Br response due only to plasma motion from I_cond (=d[Br_grid, Bz_grid,...]pl/dIcond), T/A (totbgd x ncond) ', ...
'CXmatBZ','gridpt. Bz response due only to plasma motion from I_cond (=d[Br_grid, Bz_grid,...]pl/dIcond), T/A (totbgd x ncond) ', ...
'CYmatBR','grid pt. Br response due to plasma current change (=d[Brgrid, Bzgrid]pl/dIpl) T/A (totbgd x ncond) ', ...
'CYmatBZ','grid pt. Bz response due to plasma current change (=d[Brgrid, Bzgrid]pl/dIpl) T/A (totbgd x ncond) ', ...
'CmatBRbetap','gridpt Br response due to changes in betap (=d[Bgrid]pl/dbetap), T/unit-beta (totbgd x 1) ', ...
'CmatBZbetap','gridpt Bz response due to changes in betap (=d[Bgrid]pl/dbetap), T/unit-beta (totbgd x 1) ', ...
'CmatBRli','gridpt Br response due to changes in li (=d[Bgrid ]pl/dli), T/unit-li (totbgd x 1) ', ...
'CmatBZli','gridpt Bz response due to changes in li (=d[Bgrid ]pl/dli), T/unit-li (totbgd x 1) ', ...
'nptbgd2','array of number of points in each bgrid segment (bgrdnms2) ', ...
'rbgd2','r coordinates of grid pts, grid 2 ', ...
'zbgd2','z coordinates of grid pts, grid 2 ', ...
'BR20',' ', ...
'BZ20',' ', ...
'CmatBR2',' ', ...
'CmatBZ2',' ', ...
'CXmatBR2',' ', ...
'CXmatBZ2',' ', ...
'CYmatBR2',' ', ...
'CYmatBZ2',' ', ...
'CmatBR2betap',' ', ...
'CmatBZ2betap',' ', ...
'CmatBR2li',' ', ...
'CmatBZ2li',' ', ...
'dcrzdi',' ', ...
'dcrzdip',' ', ...
'dcrzdbetap',' ', ...
'dcrzdli',' ', ...
'GRbgd_p','Greens function from plasma current on equilibrium grid to Br on Xpoint grid 1 ', ... 
'GRbgd_c', 'Greens function from coil currents to Br on Xpoint grid 1 ', ... 
'GRbgd_v',  'Greens function from vessel currents to Br on Xpoint grid 1 ', ...
'GRbgd_pf', 'Greens function from ?? to Br on Xpoint grid 1 ', ...
'GZbgd_p','Greens function from plasma current on equilibrium grid to Bz on Xpoint grid 1 ', ...
'GZbgd_c', 'Greens function from coil currentsto Bz on Xpoint grid 1 ', ... 
'GZbgd_v', 'Greens function from vessel currents to Bz on Xpoint grid 1 ', ... 
'GZbgd_pf', 'Greens function from ?? to Bz on Xpoint grid 1 ', ...
'iplcirc','flag, if=1 the model includes plasma circuit');

if(nbgrd2~=0)
   output_objs.description.GRbgd2_p = ...
 	'Greens function from plasma current on equilibrium grid to Br on Xpoint grid 1 ';
   output_objs.description.GRbgd2_c = ...
 	'Greens function from coil currents to Br on Xpoint grid 1 ';
   output_objs.description.GRbgd2_v = ...
 	'Greens function from vessel currents to Br on Xpoint grid 1 ';
   output_objs.description.GRbgd2_pf= ...
 	'Greens function from ?? to Br on Xpoint grid 1 ';
   output_objs.description.GZbgd2_p = ...
 	'Greens function from plasma current on equilibrium grid to Bz on Xpoint grid 1 ';
   output_objs.description.GZbgd2_c = ...
 	'Greens function from coil currentsto Bz on Xpoint grid 1 ';
   output_objs.description.GZbgd2_v = ...
 	'Greens function from vessel currents to Bz on Xpoint grid 1 '; 
   output_objs.description.GZbgd2_pf= ...
 	'Greens function from ?? to Bz on Xpoint grid 1 ';
end

output_objs.build_inputs = build_inputs;

eval(['save ' out_file_name ' output_objs'])

clear Miso_p GRbgd_p GZbgd_p
