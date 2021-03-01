% Magnetics Tools 
% 
% bgreens_distrib B-field between d3d parallelogram type coils and Bprobes
% bgreens_fine    Like bgreens_distrib but uses "FINE" for lowest level (fast)
% bgreens_conway  Like bgreens_distrib but uses "Conway" for lowest level
% calc_ohmic_res  Coil resistance needed to produce uniform plasma flux plateau 
% calc_bgreens    Br, Bz green functions calculation
% calc_mpp        calculate mutuals between plasma grid points
% cfil_3d_b       B-field at [X,Y,Z] from filament circle coil
% cir_filament    Axisymmetric magnetics from a point to a point
% crect2bprobe    B-probe signal from uniform current dens. rectangle coil
% crect_3d        B-field at [X,Y,Z] from rectangular X-section circle coil
% fine            Flux and B field from axisymmetric rectangular X sec coil
% fil_a           Vector potential Az from filament in Z direction, with 1 Amp
% fil_3d_b        B-field at [X,Y,Z] from straight filament
% filr_3d_b       B-field at [X,Y,Z] from straight filament with finite radius
% green_near      Br, Bz green function for point close to rectang. coil
% induc_rect      Inductance of a rectangle of sides axb and wire radius r
% mag_fil_pt      3D magnetic fields at x,y,z pt from 3-D filament
% mindbf          Mutual Inductance and Br, Bz from loop to point far from loop
% mindbf_gen      Mutual Inductance and Br,m loop, handles <close> points
% mut_fil_fil     Mutual Inductance from a filament set A to filament set B
% mutind          Mutual Inductance of two coaxial circular current loops
% mutind_near     low level function, mutuals calculation, current loops
% mutind_fine     like mutind_near but uses FINE for low level magnetics (FAST)
% mutind_conway   like mutind_near but uses "Conway" for low level magnetics
% mut_coil_set    Mutual inductance rect. axi. coil system: FINE & 3rd ord quad
% mut_ladder      computes mutual inductance of a ladder structure
% mutind_distrib  Mutuals between d3d parallelogram type coils
% ohmic_dist      Calculate Ohmic distribution in coils based on plasma shape
% rectl           self inductance for rectangular current source
% rect_2d         B-field at [X,Y] from straight rectangular X-sect. conductor
% rect_3d_b       B-field at [X,Y,Z] from box shaped straight filament
% self_fil        Self Inductance of filament set (ONLY APPROX)
% sind		  Self Inductance of a axisymmetic elongated coil (Uniform J)
