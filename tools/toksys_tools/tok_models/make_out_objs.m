 %
%  SYNTAX:  make_out_objs
%
%  PURPOSE:  make output object files using both plasma_output.m and 
%	plasma_output2.m
%
%  INPUT:
%    tokamak
%    ichooseq
%    shotnum = shot number to generate model from
%    tmodel = time of equilibrium to generate model from (ms)
%    outflag = flag to determine which set(s) of output objects to generate
%               0 = none
%               1 = run plasma_output, save results
%               2 = run plasma_output2, save results
%               3 = run both plasma_output, plasma_output2 - save (default)
%    output_irzresp = set to one of:
%             0 = use no rigid response
%             1 = use only rigid r response
%             2 = use both rigid r and z response
%    vacuum_model', 0, ...
%    iplcirc', 1, ...
%    num_Ecoils = number of E-coils (2 or 5)
%    eqdir   = directory where nominal EFIT equilibrium files are located
%    Te_res   =
%    li_res   = 
%    Rp       = plasma resistance
%    Zeff_res = 
%    isonms = names of isoflux segments
%    bgdnms = names of grids used to compute magnetic field near X points
%    scale_cc_resp = scaling vector for coil currents in response objects
%    scale_vv_resp = scaling vector for vessel currents in response objects
%
%  OUTPUT files:
%	out_eqn_objects.mat = objects from plasma_output.m
%	out2_eqn_objects.mat = objects from plasma_output2.m
%
%  RESTRICTIONS:  Assumes that drzdi files is of form 
%			drzdi_<shotnum>_<tmodel>.mat
 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	9/16/98
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  @(#)make_out_objs.m	1.4 09/02/09

include_Ip=1;

%drzdi_file = ['drzdi_' int2str(shotnum) '_' int2str(tmodel) '.mat']

build_inputs = struct( ...
   'tokamak', tokamak, ...
   'ichooseq', ichooseq, ...
   'output_irzresp',output_irzresp, ...  
   'vacuum_model', 0, ...
   'iplcirc', 1, ...
   'Te_res', Te_res, ...
   'li_res', li_res, ...
   'Rp', Rp, ...	
   'Zeff_res', Zeff_res);

if(ichooseq>0 & exist('eqdir')==1)
   if tmodel>999
      tstring = ['.0' int2str(tmodel)];
   else
      tstring = ['.00' int2str(tmodel)];
   end
   if shotnum < 100000
      efit_gfile = [eqdir '/g0' int2str(shotnum) tstring]
   else
      efit_gfile = [eqdir '/g' int2str(shotnum) tstring]
   end
   build_inputs.efit_gfile = efit_gfile;
elseif(ichooseq>0 & exist('efit_source'))
   build_inputs.efit_source = efit_source;
   build_inputs.shotnum = shotnum;
else
   wait('ERROR make_out_objs: either eqdir or efit_source must be specified')
   return;
end

% Need to customize based on how coils are assumed connected for purpose
% of equilibrium reconstruction.

switch tokamak

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Can we use the scale_cc_resp input to get rid of the D3D-specific code
% here and in plasma_output.m????
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   case {'D3D','DIII-D','DIIID','d3d','diii-d','diiid'}

      build_inputs.time = round(equil_data.time);	% time is in ms for D3D

% PCS isoflux response assumes 6 E coils and 25 vessel elements.
      if num_Ecoils>=5
         build_inputs.ccnturn = [ones(6,1); tok_data_struct.fcnturn];
      else
         build_inputs.ccnturn = [ones(2,1); tok_data_struct.fcnturn];
      end
      build_inputs.ccnturn = [tok_data_struct.ecnturn; tok_data_struct.fcnturn];

      build_inputs.calc_resp= 0;
%      build_inputs.ecid = [1 2 1 2 1 2]; % rtefit response stores 6 Ecoil segs
      build_inputs.num_Ecoils = num_Ecoils;

% RTEFIT responses are stored assuming 6 E-coils: 
%		ECOILA ECOILB E567UP E567DN E89DN E89UP 
%      build_inputs.scale_cc_resp = ones(24,1);
%      build_inputs.scale_cc_resp(1:6) = 1./[48 48 6 6 7 7];

   case {'NSTX','nstx'}

      build_inputs.time = round(equil_data.time*1000);% time is seconds for NSTX
      fcid = tok_data_struct.def_connect.fcid;
      nfc = length(unique(fcid));
      build_inputs.ccnturn = zeros(nfc,1);

      if(isfield(tok_data_struct,'ecdata') & ~isempty(tok_data_struct.ecdata))
         nec = length(unique(tok_data_struct.ecdata(5,:)));
         for k=1:nec
            build_inputs.ccnturn(k) = sum(tok_data_struct.ecdata(5,:)'==k);
         end
      else
         nec = 0;
      end
      for k=1:nfc
         idx = find(fcid==k);
         build_inputs.ccnturn(k+nec) = sum(tok_data_struct.ccnturn(idx+nec));
      end

      build_inputs.fcid = tok_data_struct.def_connect.fcid;
      build_inputs.vvid = tok_data_struct.def_connect.vvid;
      build_inputs.vvfrac = tok_data_struct.def_connect.vvfrac;
      build_inputs.calc_resp= 1;
      build_inputs.make_resp_dir = '/pcshome/walker/nstx_pcs/';
      build_inputs.num_Ecoils = 0;
      build_inputs.scale_cc_resp = 1./build_inputs.ccnturn;
 
   otherwise
      wait(['ERROR make_out_objs: tokamak = ' tokamak ' not supported'])
      return;

end

build_inputs.gridsize = 100*tok_data_struct.nr + tok_data_struct.nz;

isoflux_defns = struct( ...
'isonms', isonms, ...
'ref_num',ref_num, ...
'bgrdnms', bgrdnms);

if(exist('bgrdnms2'))
   isoflux_defns.bgrdnms2 = bgrdnms2;
end

if (outflag==1 | outflag ==3)
   output_objs = plasma_output(build_inputs, equil_data, ...
		tok_data_struct, isoflux_defns, rzrig_data);	
end

if (outflag==2 | outflag ==3)
   ichooseq = 1;
   output2_objs = plasma_output2(equil_data,rzrig_data, ...
				tok_data_struct,build_inputs);
end


