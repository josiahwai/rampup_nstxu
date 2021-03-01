function units = make_units_struct(imks,iterminal)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SYNTAX:  units = make_units_struct(imks,iterminal)
%
%  PURPOSE:  Make units data structure for standard format tokamak data struct.
%
%  INPUT:
%    imks = if 1, input objects in datafiles are mks, else 0 (default=mA,uH)
%    iterminal = if 1, input objects are terminal units, else 0 (default)
%
%  OUTPUT:
%    units = units of data objects 
  
%  RESTRICTIONS: 
%
%  METHOD:  
%
%  WRITTEN BY:  Mike Walker 	ON 	12/16/04
%
%  MODIFICATIONS:
%     2008-01-21  NWE   Added units for LV mutuals (mph,mhv,mhc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if imks 
   gbc_units = 'T/A';
   gbv_units = 'T/A';
   gpb_units = 'T/A';
   gbr2p_units = 'T/A';
   gbz2p_units = 'T/A';
   mcc_units = 'Wb/A';
   mcv_units = 'Wb/A';
   mvv_units = 'Wb/A';
   mpc_units = 'Wb/A';
   mpv_units = 'Wb/A';
   mlc_units = 'Wb/A';
   mlv_units = 'Wb/A';
   mpl_units = 'Wb/A';
   mhc_units = 'Wb/A';
   mhv_units = 'Wb/A';
   mph_units = 'Wb/A';
   mpp_units = 'Wb/A';
   rc_units =  'Ohms';
   rv_units =  'Ohms';
else
   gbc_units = 'T/MA';
   gbv_units = 'T/MA';
   gpb_units = 'T/MA';
   gbr2p_units = 'T/MA';
   gbz2p_units = 'T/MA';
   mcc_units = 'Wb/MA';
   mcv_units = 'Wb/MA';
   mvv_units = 'Wb/MA';
   mpc_units = 'Wb/MA';
   mpv_units = 'Wb/MA';
   mpp_units = 'Wb/MA';
   mlc_units = 'Wb/MA';
   mlv_units = 'Wb/MA';
   mpl_units = 'Wb/MA';
   mhc_units = 'Wb/MA';
   mhv_units = 'Wb/MA';
   mph_units = 'Wb/MA';
   rc_units =  'micro-Ohms';
   rv_units =  'micro-Ohms';
end

if iterminal
   gbc_units = [gbc_units ' (terminal)'];
   mcc_units = [mcc_units ' (terminal)'];
   mpc_units = [mpc_units ' (terminal)'];
   mlc_units = [mlc_units ' (terminal)'];
   mhc_units = [mhc_units ' (terminal)'];
   rc_units =  [rc_units  ' (terminal)'];
else
   gbc_units = [gbc_units '-turn']; 
   mcc_units = [mcc_units '-turn']; 
   mpc_units = [mpc_units '-turn']; 
   mlc_units = [mlc_units '-turn']; 
   mhc_units = [mhc_units '-turn']; 
   rc_units =  [rc_units  ' (lumped)'];
end

units = struct( ...
'gbc', gbc_units, ...
'gbv', gbv_units, ...
'gpb', gpb_units, ...
'gbr2p',gbr2p_units, ...
'gbz2p',gbz2p_units, ...
'mcc', mcc_units, ...
'mcv', mcv_units, ...
'mvv', mvv_units, ...
'mpc', mpc_units, ...
'mpv', mpv_units, ...
'mlc', mlc_units, ...
'mlv', mlv_units, ...
'mpl', mpl_units, ...
'mhc', mlc_units, ...
'mhv', mlv_units, ...
'mph', mpl_units, ...
'mpp', mpp_units, ...
'resc', rc_units, ...
'resv', rv_units);
