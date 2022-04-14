% build the output cmat matrix, dy = cmat * dx
% Note that this eqn is valid only in the linearized frame

function response = vacuum_response(pla, targets, tok_data_struct)

circ = nstxu2016_circ(tok_data_struct);
struct_to_ws(tok_data_struct);
mpc = tok_data_struct.mpc;

% coil currents response to themselves (identity)
disdis = eye(circ.nx);
      
% vacuum flux-response == mutual inductances
% dcphidip = pla.pcurrt(:) / sum(pla.pcurrt(:));
% mp_ip = mpp * dcphidip;
% mp_ip = mpp * dcphidip * 0;
mp_ip = zeros(nz*nr, 1);

dpsizrdis = [mpc mpv mp_ip] * circ.Pxx;


% flux response at control points and boundary-defining point
for j = 1:circ.nx
  dpsizr = reshape(dpsizrdis(:,j), nz, nr);  
  dpsicpdis(:,j) = bicubicHermite(rg, zg, dpsizr, targets.rcp, targets.zcp);  
  [dpsibrydis(j), dpsibrydis_r(j), dpsibrydis_z(j)]  = bicubicHermite(rg, zg, dpsizr, targets.rbdef, targets.zbdef);      
end

response = variables2struct(disdis, dpsicpdis, dpsibrydis, dpsibrydis_r, dpsibrydis_z); 

end





















































