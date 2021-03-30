% build the output cmat matrix, dy = cmat * dx
% Note that this eqn is valid only in the linearized frame

function cmat_data = build_cmat_data_vacuum(eq, circ, tok_data_struct, targ_geo)

nz = tok_data_struct.nz;
nr = tok_data_struct.nr;
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;

% coil currents response (identity)
cmat_data.x = eye(circ.nx);

% x-point calcs
[rx, zx] = isoflux_xpFinder(eq.psizr, 0.65, -1.1, eq.rg, eq.zg);
[~, ~, ~, psi_rr, psi_zz, psi_rz] = bicubicHermite(eq.rg, eq.zg, eq.psizr, rx, zx);
Jinv = inv([psi_rr psi_rz; psi_zz psi_rz]);
      
% vacuum flux-response == mutual inductances
dcphidip = eq.pcurrt(:) / sum(eq.pcurrt(:));
mp_ip = tok_data_struct.mpp * dcphidip;

dpsizrdis = [tok_data_struct.mpc tok_data_struct.mpv mp_ip];

% control points and boundary
for j = 1:circ.ns
  dpsizr = reshape(dpsizrdis(:,j), nz, nr);
  
  dpsicpdis(:,j) = bicubicHermite(rg, zg, dpsizr, targ_geo.cp.r, targ_geo.cp.z);
  
  [dpsibrydis(j), dpsirdis, dpsizdis]  = bicubicHermite(rg, zg, dpsizr, rx, zx);
    
  drzxdis = -Jinv * [dpsirdis; dpsizdis]; 
  drxdis(j) = drzxdis(1);
  dzxdis(j) = drzxdis(2);
end

cmat_data.dpsicpdix = (dpsicpdis - dpsibrydis) * circ.Pxx;
cmat_data.dpsibrydix = dpsibrydis * circ.Pxx;
cmat_data.drx1dix = drxdis * circ.Pxx;
cmat_data.dzx1dix = dzxdis * circ.Pxx;



























































