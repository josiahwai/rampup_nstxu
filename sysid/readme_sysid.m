
% DESCRIPTION OF CONTENTS:

%%
% To fit the plasma resistance, use fit_plasma_resistance.m - this fits the
% plasma resistance Rp(t) and inductance Lp(t) in order to match the ip
% dynamics. For 2015-16 campaign Rp(t) and Lp(t) vary from shot to shot 
% by ~factor of 2. For an accurate fit, would need to fit directly to that
% shot. 

%%
% To fit the vacuum system, the workflow is:

% 1. run extFit_SysId.m - fits external power supply resistances 
% and inductances based on dedicated power supply test shots. Results are 
% saved in ext_fit.mat
 
% 2. run coil_vessel_sysid.m - uses the prev fit and modifies the coil AND vessel
% inductances and resistances, using non-dedicated shots. Parameters are
% adjusted so that simulated coil and vessel currents match their efit01
% values. This can take >1hr to run. Results are saved in
% coil_vessel_fit.mat

%%
% To verify the combined model, run verify_sysid_fits.m. This simulates with
% the fitted parameters and compares it to actual efit data. 

% Josiah Wai 11/2021














