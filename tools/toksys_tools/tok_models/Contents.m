%
%  Tokamak System Modeling Tools (Humphreys)
% 
%  build_model         - 
%  build_tokamak_system - Function to build state space for control design and
%                        analysis. Calculates vertical growth rates along the
%                        way. Note that build_tok_sys is intended to be called
%                        by a machine-specific script like build_east_sys.m
%  calc_ssop_from_I    - Compute steady-state operating point (ynom,Inom,Vnom)
%  cc_efit_to_tok      - convert EFIT conductor currents to toksys
%  change_tokobj_units - modify units of vacuum data objects in tok_data_struct
%  gbr2p_x_vec         - multiply compressed object gbr2p by a vector
%  get_signals         - Fetch signals in "output_signals" field of system model struct
%  lim2vv              - Construct a vacuum vessel model from limiter surface specification
%  makelpcur           -  make elliptical cross-section, flat to parabolic plasma 
%                         current distr. defined on standard EPGenv grid.
%  make_out_objs - 
%  make_tok_objects    - Function to calculate all Green functions and other
%                        data objects needed for build_tokamak_system.m.
%                        Note that make_tok_objects is intended to be called by 
%                        a machine-specific script like make_east_objects.m
%  mpp_x_vec           - multiply compressed object mpp or gbz2p by a vector
%  plasma_output       - 
%  plasma_output2      - 
%  plasma_out_common   - common code for both plasma_output, plasma_output2
%  proj_turn           - Make proj. vector/matrix for coil connections and turn specs.
%  rzrig               - creates plasma rigid model response objects
%  scale_equil_response- scale plasma response model with equilibrium Ip
%  load_tok_objects    - Load vacuum data objects for selected device configuration
%  plot_tok_geo_config - call plot_tok_geo for multiple vacuum data objects
%
%  gspert   - return linear plasma response
%  gseq     - return an equilibrium, given the conductor currents and pprime, ffprim
%  gsevolve - like gseq but also returns time derivatives
