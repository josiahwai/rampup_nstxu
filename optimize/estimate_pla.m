% INPUTS:
% bry - struct with fields 'r' and 'z', position of target boundary points
% ip - target plasma current
% wmhd - target plasma thermal energy
% tok_data_struct - struct with machine geometry
% eq0 - optional struct with previous iteration of fluxes, not used if
%       opts.cold_start = 1
% opts - optional struct with fields: plotit - flag to make plots,
%        cold_start - flag, if (1) then estimate flux on grid, if (0) then use
%        grid flux from eq0
%
%
% OUTPUT
% pla - struct containing estimates of the plasma current and flux distribution
%       an estimate

function pla = estimate_pla(target, tok_data_struct, eq0, opts)
% ccc
% load('matlab.mat')
% opts.plotit = 1;
% opts.cold_start = 0;

if ~exist('opts', 'var'), opts = struct; end
if ~isfield(opts, 'plotit'), opts.plotit = 0; end
if ~isfield(opts, 'cold_start'), opts.cold_start = 0; end
if ~isfield(opts, 'debug_use_actual'), opts.debug_use_actual = 0; end



circ = nstxu2016_circ(tok_data_struct);
struct_to_ws(tok_data_struct);
mpc = tok_data_struct.mpc;
dr = mean(diff(rg));
dz = mean(diff(zg));
dA = dr * dz;


%%
if opts.cold_start

  % ================================================================
  % Cold start: estimate plasma current distribution without knowing
  % previous iteration of equilibrium.
  % ================================================================
  % On the first iteration, need to estimate the plasma flux without knowing
  % the previous iteration info (boundary, external flux, etc). Use the
  % target boundary and target Ip to find a plasma current distribution:
  % j = jhat * (1 - x^a)^b, where x is the distance ratio from centroid to
  % boundary (0 at centroid, 1 at boundary), a and b are constants, and jhat
  % is a constant scaled to match the target Ip.

  i = isnan(target.rcp);
  target.rcp(i) = [];
  target.zcp(i) = [];

  % sort and interpolate target boundary vertices, for well-defined polygon
  rz = solveTSP([target.rcp(:) target.zcp(:)], 0); % sort via traveling salesman
  rz = [rz; rz(1,:)];
  bry = fnplt(cscvn(rz'));  % spline interpolation
  rbbbs = bry(1,:);
  zbbbs = bry(2,:);
  P = polyshape(rbbbs, zbbbs);

  [rc, zc] = centroid(P);

  % at each grid point, find normalized distance between centroid and
  % boundary
  rgridc = rgg + dr/2;
  zgridc = zgg + dz/2;
  dist2cent = sqrt((rgridc(:)-rc).^2 + (zgridc(:)-zc).^2);
  [xy,dist2bry] = distance2curve([rbbbs(:) zbbbs(:)], [rgridc(:) zgridc(:)]);
  x = dist2cent ./ (dist2cent + dist2bry);

  in = inpolygon(rgridc(:), zgridc(:), P.Vertices(:,1), P.Vertices(:,2));
  x(~in) = nan;
  x = reshape(x, nr, nz);

  ip = target.ip;

  a = 1.87;
  b = 1.5;
  j0 = (1-x.^a).^b;
  jhat = ip / (nansum(j0(:)) * dr * dz);
  jphi = j0 * jhat;
  jphi(isnan(jphi)) = 0;
  pcurrt = jphi * dr * dz;
  jphi = jphi / 1e6; % A/m^2 -> MA/m^2

  psizr_pla = mpp * pcurrt(:);
  psizr_pla = reshape(psizr_pla, nr, nz);

  % Estimate the applied flux & total flux
  % weighted least squares fit to find coil currents that would create this
  % boundary. A*Ic = b.

  dum.pcurrt = pcurrt;
  r = vacuum_response(dum, target, tok_data_struct);

  A = [eye(circ.ncx_keep);
    r.dpsicpdis(:,circ.iicx_keep)-r.dpsibrydis(:,circ.iicx_keep);
    r.dpsibrydis_r(circ.iicx_keep);
    r.dpsibrydis_z(circ.iicx_keep)];

  [psibry_pla, psibry_pla_r, psibry_pla_z] = bicubicHermite(...
    rg,zg,psizr_pla,target.rbdef, target.zbdef);

  psi_cp_err = bicubicHermite(rg,zg,psizr_pla,target.rcp, target.zcp)';

  b = -[zeros(circ.ncx_keep,1);
    psi_cp_err - psibry_pla;
    psibry_pla_r;
    psibry_pla_z];

  weights = diag([ones(circ.ncx_keep,1)*1e-6; ...
    ones(length(target.rcp), 1); ...
    ones(2,1)]);

  icx = pinv(weights*A)*(weights*b);
  mpcx = mpc*circ.Pcc;
  mpcx = mpcx(:,circ.iicx_keep);

  psizr_app = reshape(mpcx * icx, nr, nz);

  psizr = psizr_app + psizr_pla;

  if opts.plotit    
    figure
    hold on
    contour(rg,zg,psizr_app+psizr_pla,50)
    scatter(target.rcp, target.zcp)
    scatter(target.rbdef, target.zbdef)
  end
  
  rbbbs = target.rcp;
  zbbbs = target.zcp;
  pla = variables2struct(psizr_pla, pcurrt, jphi, rbbbs, zbbbs);

  eq_opts.plotit = 0;
  eq = eq_analysis(psizr, pla, tok_data_struct, eq_opts);
  
else % warm start
  eq = eq0;  
end

%%
if ~opts.debug_use_actual

  % Warm start: now that we have an estimate for plasma flux on 
  % grid, can use P' & FF' profiles to update the equilibrium (plasma current
  % distribution. We will use standard P' and FF' profiles and scale these
  % (2 free parameters) in order to match the target Wmhd and target Ip.

  mu0 = pi*4e-7;
  psin = linspace(0,1,nr);

  % find boundary
  %   rz = solveTSP([target.rcp(:) target.zcp(:)], 0); % sort via traveling salesman
  %   rz = [rz; rz(1,:)];
  %   bry = fnplt(cscvn(rz'));  % spline interpolation
  %   rbbbs = bry(1,:);
  %   zbbbs = bry(2,:);
  rbbbs = eq.rbbbs;
  zbbbs = eq.zbbbs;

  in = inpolygon(rgg, zgg, eq.rbbbs, eq.zbbbs);
  psin_grid = (eq.psizr - eq.psimag) / (eq.psibry - eq.psimag);
  psin_grid(~in) = nan;

  % profiles to be scaled
  % notation: zero suffix is the unscaled quantity
  [pprime0, ffprim0] = load_standard_efit_profiles;


  % scale P' profile to match Wmhd
  % ------------------------------
  pres0 = cumtrapz(psin, pprime0) * (eq.psibry - eq.psimag) / (2*pi);  % pressure profile to be scaled
  pres0 = pres0-pres0(end);
  pres0_grid = interp1(psin, pres0, psin_grid(:));

  tmp_wmhd = nansum(3 * pi * rgg(:) .* pres0_grid * dA);

  pprime_coef = target.wmhd / tmp_wmhd;

  pres = pprime_coef * pres0;
  pprime = pprime_coef * pprime0;

  pprime_grid = interp1(psin, pprime, psin_grid(:));
  pprime_grid = reshape(pprime_grid, nr, nz);


  % scale FF' profile to match Ip
  % ------------------------------

  % current from pprime term
  jphi_pressure = rgg .* pprime_grid;
  pcurrt_pressure = jphi_pressure * dA;
  ip_pressure = nansum(pcurrt_pressure(:));

  % current from ffprim term, use target Ip to scale the FF' profile
  ffprim0_grid = interp1(psin, ffprim0, psin_grid(:));
  ffprim0_grid = reshape(ffprim0_grid, nr, nz);

  jphi0_ffprim = ffprim0_grid ./ (rgg * mu0);
  pcurrt0_ffprim = jphi0_ffprim * dA;
  ip0_ffprim = nansum(pcurrt0_ffprim(:));

  ip_ffprim = target.ip - ip_pressure;

  ffprim_coeff = ip_ffprim / ip0_ffprim;

  ffprim = ffprim0 * ffprim_coeff;
  ffprim_grid = ffprim0_grid * ffprim_coeff;

  % current from scaled profiles
  jphi = rgg .* pprime_grid + ffprim_grid ./ (rgg * mu0);
  jphi(isnan(jphi)) = 0;
  pcurrt = jphi * dr * dz;
  jphi = jphi / 1e6;   % A/m^2 -> MA/m^2

  % write to output
  psizr_pla = mpp * pcurrt(:);
  psizr_pla = reshape(psizr_pla, nr, nz);
  area = polyarea(rbbbs, zbbbs);


else   % DEBUGGING USE ONLY (opts.debug_use_actual = 1)
  pcurrt = opts.ref_eq.pcurrt;
  psizr_pla = reshape(mpp * pcurrt(:), nr, nz);
  rbbbs = opts.ref_eq.rbbbs;
  zbbbs = opts.ref_eq.zbbbs;
  i = rbbbs ~= 0 & zbbbs ~= 0;
  rbbbs = rbbbs(i);
  zbbbs = zbbbs(i);
  area = polyarea(rbbbs, zbbbs);
  jphi = opts.ref_eq.jphi;
  pprime = opts.ref_eq.pprime;
  ffprim = opts.ref_eq.ffprim;
  pres = opts.ref_eq.pres;
end

pla = variables2struct(psizr_pla, pcurrt, jphi, area, rbbbs, zbbbs);


%%
% make some plots
if opts.plotit
  
  % plot profiles
  figure
  hold on
  plot(pprime)
  try plot(opts.ref_eq.pprime), catch; end
  legend('Estimate', 'EFIT', 'fontsize', 18)
  
  figure
  hold on
  plot(pres)
  try plot(opts.ref_eq.pres), catch; end
  legend('Estimate', 'EFIT', 'fontsize', 18)

  figure
  hold on
  plot(ffprim)
  try plot(opts.ref_eq.ffprim), catch; end
  legend('Estimate', 'EFIT', 'fontsize', 18)
  
  figure
  subplot(121)
  contourf(jphi)
  colorbar
  title('Estimate')
  
  try
  subplot(122)
  contourf(opts.ref_eq.jphi)
  colorbar
  title('EFIT')
  catch 
  end
  
  try    
    psi_app = reshape(mpc * opts.ref_eq.ic + mpv*opts.ref_eq.iv, nr, nz);
    psizr = psi_app + psizr_pla;
    
    figure
    hold on
    
%     [~,cs] = contour(rg,zg,psizr,20,'--r');
%     contour(rg,zg,opts.ref_eq.psizr,cs.LevelList, 'b')
    
    psibry = bicubicHermite(rg, zg, opts.ref_eq.psizr, target.rbdef, target.zbdef);
    contour(rg, zg, opts.ref_eq.psizr, [psibry psibry], 'b')
    
    psibry = bicubicHermite(rg, zg, psizr, target.rbdef, target.zbdef);
    contour(rg, zg, psizr, [psibry psibry], '--r')
    
    axis equal
    set(gcf, 'Position', [992 183 431 622])
    scatter(target.rcp, target.zcp, 'k', 'filled')
    
    
  catch
  end
  
end













