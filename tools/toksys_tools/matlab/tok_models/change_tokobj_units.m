function tok_data_struct = change_tokobj_units(tok_data_struct,imks,iterminal)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE: mod_data_struct=change_tokobj_units(tok_data_struct,imks,iterminal)
%
%  PURPOSE: Create a copy of the tokamak vacuum data objects data structure
%		with specified units and turns convention.
%
%  INPUTS:
%    tok_data_struct = data structure, contents defined by creation in
%		 	make_tok_objects.m
%    imks = if 1, create data objects in MKS units, 0 means currents in MA
%    iterminal = if 1, all multi-turn coils are interpreted to be in
%		"terminal mode", e.g. currents are as measured at input or
%		output terminals, 0 means "lumped mode", e.g. current is
%		total cross-sectional current in multi-turn coil.
%
%  OUTPUTS:
%    mod_data_struct = same data structure, but with modified data objects

%  RESTRICTIONS:
% 
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker	ON	5/4/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(nargin<3)
   wait('ERROR change_tokobj_units: must enter 3 input arguments')
end

if(imks==0 & tok_data_struct.imks==1)
   units_multiplier = 1e6;
elseif(imks==1 & tok_data_struct.imks==0)
   units_multiplier = 1e-6;
else
   units_multiplier = 1;
end

if(iterminal==0 & tok_data_struct.iterminal==1)
   turnmat = inv(diag(tok_data_struct.ccnturn));
   ccnturn_mult = 1./(tok_data_struct.ccnturn.^2);
elseif(iterminal==1 & tok_data_struct.iterminal==0)
   turnmat = diag(tok_data_struct.ccnturn);
   ccnturn_mult = tok_data_struct.ccnturn.^2;
else
   turnmat = 1;
   ccnturn_mult = 1;
end

tok_data_struct.resv = tok_data_struct.resv * units_multiplier;
tok_data_struct.mvv  = tok_data_struct.mvv * units_multiplier;
tok_data_struct.mpv  = tok_data_struct.mpv * units_multiplier;
tok_data_struct.mlv  = tok_data_struct.mlv * units_multiplier;
tok_data_struct.gbv  = tok_data_struct.gbv * units_multiplier;
tok_data_struct.mpl  = tok_data_struct.mpl * units_multiplier;
tok_data_struct.gpb  = tok_data_struct.gpb * units_multiplier;
tok_data_struct.mpp  = tok_data_struct.mpp * units_multiplier;
tok_data_struct.gbr2p = tok_data_struct.gbr2p * units_multiplier;
tok_data_struct.gbz2p = tok_data_struct.gbz2p * units_multiplier;

tok_data_struct.resc = tok_data_struct.resc * units_multiplier; 
tok_data_struct.mcv  = tok_data_struct.mcv * units_multiplier;
tok_data_struct.mcc  = tok_data_struct.mcc * units_multiplier;
tok_data_struct.mpc  = tok_data_struct.mpc * units_multiplier; 
tok_data_struct.mlc  = tok_data_struct.mlc * units_multiplier; 
tok_data_struct.gbc  = tok_data_struct.gbc * units_multiplier;
tok_data_struct.gbr2c = tok_data_struct.gbr2c * units_multiplier; 
tok_data_struct.gbz2c = tok_data_struct.gbz2c * units_multiplier;
tok_data_struct.gbr2v = tok_data_struct.gbr2v * units_multiplier;
tok_data_struct.gbz2v = tok_data_struct.gbz2v * units_multiplier;

tok_data_struct.resc = tok_data_struct.resc.*ccnturn_mult;
tok_data_struct.mcv = turnmat*tok_data_struct.mcv;
tok_data_struct.mcc = turnmat*tok_data_struct.mcc*turnmat;
tok_data_struct.mpc = tok_data_struct.mpc*turnmat;
tok_data_struct.mlc = tok_data_struct.mlc*turnmat;
tok_data_struct.gbc = tok_data_struct.gbc*turnmat;
tok_data_struct.gbr2c = tok_data_struct.gbr2c*turnmat;
tok_data_struct.gbz2c = tok_data_struct.gbz2c*turnmat;

tok_data_struct.imks = imks;
tok_data_struct.iterminal = iterminal;
