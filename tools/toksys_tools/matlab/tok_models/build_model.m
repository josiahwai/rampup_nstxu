 %
%  SYNTAX:  build_model
%
%  PURPOSE:  Generate data objects which can be used to construct complete
%	system model, based on linearized plasma response.
%
%  INPUT:
%	shot    = shot number
%	tmodel  = time to derive model from (in ms)
%	tokamak = name of device to build model for
%	tok_data_struct = vacuum data objects for this device
%	output_irzresp = irzresp for output equation object generation
%	eqdir   = directory where efit equilibrium files are located
%	efit_source = tree in mdsplus to get equilibrium from
%		(only one of either eqdir or efit_source should be specified)
%	outflag = flag to determine what output files to produce
%		0 = none
%		1 = run plasma_output, save results
%		2 = run plasma_output2, save results
%		3 = run both plasma_output, plasma_output2 - save (default)
%       Rp = plasma resistance (only needed if outflag>0)
%
%  OUTPUT: files = 
%	drzdi_<shot>_<tmodel>.mat
%	out_objs_<shot>_<tmodel>.mat
%	out2_objs_<shot>_<tmodel>.mat
  
%  RESTRICTIONS:  Moving files with escapes to the operating system does not
%	always work with matlab 5.2 - do this with external commands to be sure.
% 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker   ON      ??/98
%  UPDATED BY:  Mike Walker 	ON 	7/13/09	to generalize to all tokamaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)build_model.m	1.4 09/02/09

if(~exist('eqdir') & ~exist('efit_source'))
  wait('ERROR build_model: either eqdir or efit_source must be specified')
end

% from D3D version:
%	idoe6   = set to 0 for shots between 11/97 and ??, otherwise 1
%load /m/GAtools/tokamaks/d3d/make/after_ADP/d3d_obj_mks_struct_6565.mat

if ~exist('makeout')
   makeout = 3;
end

izstab = 0 
shotnum=shot;

filetype = 2;

if(~exist('Te_res')), Te_res = 4000; end;   %eV, for plasma r
if(~exist('li_res')), li_res = 0.5; end;     %self inductance
if(~exist('Zeff_res')), Zeff_res = 1.5; end;    %Zeff for res

if exist('eqdir')
   if(tmodel<1000)
      strtemp = [int2str(shot) '.00' int2str(tmodel)];
   else
      strtemp = [int2str(shot) '.0' int2str(tmodel)];
   end
   if(shotnum < 100000)
      efit_gfile = [eqdir '/g0' strtemp];
      efit_afile = [eqdir '/a0' strtemp];
   else
      efit_gfile = [eqdir '/g' strtemp];
      efit_afile = [eqdir '/a' strtemp];
   end

   filename = efit_gfile
   equil_data = read_gfile_tok(filename,tokamak);

   filename = efit_afile
   read_afile
else
   equil_data = read_mds_eqdsk(shotnum,tmodel,efit_source,tokamak);
end

% Derive equilibrium parameters from EFIT data:
ip0 = equil_data.cpasma*1e-6;   %nominal Ip (MA)
if 0
z0 = zmagx*0.01;     %nominal vertical position (m)
r0 = rmagx*0.01;     %nominal major radius  (m)
a0 = aout*0.01;      %equilibrium minor radius (m)
kap0 = eout;         %equilibrium X-point elongation
betap = betap;       %equilibrium Betap
li = ali;            %equilibrium li
end

num_Ecoils=length(tok_data_struct.ecnturn);

mu0 = 0.4*pi;
twopi = 2*pi;
twopir = twopi*tok_data_struct.rgg(:);

dz=abs(tok_data_struct.zg(2)-tok_data_struct.zg(1));
dr=abs(tok_data_struct.rg(2)-tok_data_struct.rg(1));

[rzrig_data,cc0,vc0,dbg_objs]=rzrig(equil_data,tokamak,tok_data_struct,1);

%norm(drdi-drdix)
%norm(dzdi-dzdis)
%wait('debug 1')

if outflag > 0
   irzresp = output_irzresp;
   ichooseq = 1;
   make_out_objs
end

%norm(drdi-drdix)
%norm(dzdi-dzdis)
%wait('debug 2')

plotfile = [int2str(shot) '_equil_' int2str(tmodel) '.psc'];
figure(1)

write_summary_file
