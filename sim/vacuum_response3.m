% build the output cmat matrix, dy = cmat * dx
% Note that this eqn is valid only in the linearized frame

function response = vacuum_response3(target, tok_data_struct)

circ = nstxu2016_circ(tok_data_struct);
struct_to_ws(tok_data_struct);
mpc = tok_data_struct.mpc;

% coil currents response to themselves (identity)
disdis = eye(circ.ncx + circ.nvx);
      
% vacuum flux-response == mutual inductances
dpsizrdis = [mpc mpv] * circ.Pxx(1:end-1,1:end-1);


% flux response at control points and boundary-defining point
for j = 1:(circ.ncx+circ.nvx)
  dpsizr = reshape(dpsizrdis(:,j), nz, nr);  
  dpsicpdis(:,j) = bicubicHermite(rg, zg, dpsizr, target.rcp, target.zcp);  
  [dpsibrydis(j), dpsibrydis_r(j), dpsibrydis_z(j)] = bicubicHermite(rg, zg, dpsizr, target.rbdef, target.zbdef); 
  [dpsixdis(j), dpsixdis_r(j), dpsixdis_z(j)] = bicubicHermite(rg, zg, dpsizr, target.rx, target.zx); 
  dpsi_ingapdis(j) = bicubicHermite(rg, zg, dpsizr, target.r_ingap, target.z_ingap); 
end

response = variables2struct(disdis, dpsicpdis, dpsibrydis, dpsibrydis_r, dpsibrydis_z, dpsixdis, dpsixdis_r, dpsixdis_z, dpsi_ingapdis); 

end





















































