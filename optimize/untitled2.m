% build the output cmat matrix, dy = cmat * dx
% Note that this eqn is valid only in the linearized frame

% function cmat_data = build_cmat_data_vacuum(eq, circ, tok_data_struct, targ_geo)

ccc
load('matlab.mat')

circ = nstxu2016_circ(tok_data_struct);
struct_to_ws(tok_data_struct);

% coil currents response (identity)
cmat_data.x = eye(circ.nx);
      
% vacuum flux-response == mutual inductances
dcphidip = eq.pcurrt(:) / sum(eq.pcurrt(:));
mp_ip = mpp * dcphidip;

dpsizrdis = [mpc mpv mp_ip];


% control points and boundary
for j = 1:circ.ns
  dpsizr = reshape(dpsizrdis(:,j), nz, nr);  
  [dpsicpdis(:,j), dpsicpdis_r(:,j), dpsicpdis_z(:,j)] = bicubicHermite(rg, zg, dpsizr, targets.rcp, targets.zcp);  
  [dpsibrydis(j), dpsibrydis_r(j), dpsibrydis_z(j)]  = bicubicHermite(rg, zg, dpsizr, targets.rbdef, targets.zbdef);      
end

cmat_data.dpsicpdix = (dpsicpdis - dpsibrydis) * circ.Pxx;
cmat_data.dpsibrydix = dpsibrydis * circ.Pxx;
cmat_data.drx1dix = drxdis * circ.Pxx;
cmat_data.dzx1dix = dzxdis * circ.Pxx;



























































