% Given the plasma current distribution and info about target boundary,
% estimate the coil currents/applied flux.

function app = update_psi_app(pla, target, tok_data_struct, opts)

struct_to_ws(tok_data_struct);
mpc = tok_data_struct.mpc;
circ = nstxu2016_circ(tok_data_struct);
psizr_pla = pla.psizr_pla;

r = vacuum_response(pla, target, tok_data_struct);

if opts.use_vessel_currents
  i = [circ.iicx circ.iivx];
  n = circ.ncx + circ.nvx;
else
  i = circ.iicx;
  n = circ.ncx;
end

% form the equations A*currents = b
% b = flux error at control points, flux gradient at x-point
% A = response of errors to the currents

A = [eye(n);  
  r.dpsicpdis(:,i) - r.dpsibrydis(:,i);
  r.dpsibrydis_r(i);
  r.dpsibrydis_z(i)];

[psibry_pla, psibry_pla_r, psibry_pla_z] = bicubicHermite(...
  rg,zg,psizr_pla,target.rbdef, target.zbdef);

psi_cp_err = bicubicHermite(rg,zg,psizr_pla,target.rcp, target.zcp)';

b = -[zeros(n,1);
  psi_cp_err - psibry_pla;
  psibry_pla_r;
  psibry_pla_z];

% weights
wt.ic = ones(1,circ.ncx) * 1e-6;
wt.ic(circ.iicx_remove) = 1e8;
wt.iv = [];
if opts.use_vessel_currents
  wt.iv = ones(1,circ.nvx) * 1e-6;
end
wt.cp = ones(1,length(target.rcp));
if target.islimited
  wt.x = [0 0];
else
  wt.x = [1 1];
end
weights = diag([wt.ic wt.iv wt.cp wt.x]);

% solve weighted least squares problem to obtain currents
currents = pinv(weights*A)*(weights*b);

% measure applied fluxes
icx = currents(circ.iicx);
mpcx = mpc*circ.Pcc;
psizr_app = mpcx*icx;

if opts.use_vessel_currents
  ivx = currents(circ.iivx);
  mpvx = mpv*circ.Pvv;
  psizr_app = psizr_app + mpvx*ivx;
else
  ivx = zeros(circ.nvx,1);
end

psizr_app = reshape(psizr_app, nr, nz);

psizr = psizr_app + psizr_pla;

app = variables2struct(psizr, psizr_app, psizr_pla, icx, ivx);













