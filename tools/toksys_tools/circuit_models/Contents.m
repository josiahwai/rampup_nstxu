% Conductor circuit modeling
%  
% User functions:
%  BJacobians
%  BJacobians_source
%  bld_filter          - Build a 4-pole Bessel or Chebyshev filter model
%  cccirc_to_netlist   - convert cccirc coil circuit description to netlist description
%  coil_mutuals        - Calculate mut inductance between pairs of coils/elements.
%  coilR               - Calculate resistance of toroidal current carrying element
%  compare_fresp
%  compare_multimodel
%  compare_new
%  compare_sim
%  dynamic_equil.m
%  ecoil_turns.dat
%  elt_to_corners     - Convert "element form" ([Z;R;dZ;dR;AC;AC2]) to set of element corners 
%  get_plasma_greens
%  inverse_filter     - Inverse filter data filtered by one of PCS antialiasing filters. 
%  Jacsind
%  load_netlist       - load a netlist from a netlist file
%  make_units_struct
%  M_B_vs_ecoils
%  M_B_vs_fcoils
%  M_B_vs_vessel
%  MJacobians
%  model_from_netlist - generate ABCD model from netlist input 
%  msecalc            - Calculate mutuals from SL (saddle loops) to E coils.
%  msfcalc            - Calculate mutuals from SL (saddle loops) to F coils.
%  mutind_sl          - Calculate mut inductance between a coil and saddle loop. 
%  plot_mpp
%  read_SL
%  val_full_model
%
% Helper functions:
%  echelon            - Define branches in tree and co-tree in netlist modeling
%  makcturnzr         - Discretization of coils/vessel for calc of Green fn or mut induct
%  make_griddata      - Construct the plasma filaments grid
%  make_incidence
%  partn_cond         - partition conductors into "turns" mut induct and Green fns calcs
%  test_mpp_x_vec
%  update_constraint_eqns
