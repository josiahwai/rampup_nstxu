% build the output cmat matrix, dy = cmat * dx
% Note that this eqn is valid only in the linearized frame

function cmat_data = build_cmat_data(eq, circ, tok_data_struct, targ_geo)

nz = tok_data_struct.nz;
nr = tok_data_struct.nr;
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;

% coil currents response (identity)
cmat_data.x = eye(circ.nx);

r = gspert(eq,tok_data_struct);

dpsizrdis = [r.dpsizrdis r.dpsizrdip];
dpsibrydis = [r.dpsibrydis r.dpsibrydip];
  
% change in flux error at control points
for j = 1:circ.ns
  dpsizr = reshape(dpsizrdis(:,j), nz, nr);
  dpsicpdis(:,j) = bicubicHermite(rg, zg, dpsizr, targ_geo.cp.r, targ_geo.cp.z);
end
cmat_data.dpsicpdix = (dpsicpdis - dpsibrydis) * circ.Pxx;
cmat_data.dpsibrydix = dpsibrydis * circ.Pxx;

% motion of primary x-point
drxdix = [r.drxdis r.drxdip] * circ.Pxx;
dzxdix = [r.dzxdis r.dzxdip] * circ.Pxx;
cmat_data.drx1dix = drxdix(1,:);
cmat_data.dzx1dix = dzxdix(1,:);



fn = fieldnames(cmat_data);
for i = 1:length(fn)
  cmat_data.(fn{i}) = double(cmat_data.(fn{i}));
end

























































