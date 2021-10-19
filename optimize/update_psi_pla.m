% Given the previous iteration of equilibrium, update the plasma
% current distribution. If cold-starting (no equilibrium given), then
% performs a very rough estimate. 

function pla = update_psi_pla(eq, target, tok_data_struct, opts)


struct_to_ws(tok_data_struct);
dr = mean(diff(rg));
dz = mean(diff(zg));
dA = dr*dz;

if opts.cold_start || ~exist('eq', 'var')  
  
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

  % sort and interpolate target boundary vertices, to get a very rough
  % estimate of the plasma boundary
  i = isnan(target.rcp);
  target.rcp(i) = [];
  target.zcp(i) = [];
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
  
  pla = variables2struct(psizr_pla, pcurrt, rbbbs, zbbbs);
  
else % warm start
  
  % =====================================================================
  % Find plasma current distribution according to Grad-Shafranov equation
  % =====================================================================
  % P' and FF' profiles are found by scaling 2 standard profiles
  % characteristic of EFIT01, in order to match Ip, Wmhd targets
  
  mu0 = pi*4e-7;
  psin = linspace(0,1,nr);
  
  in = inpolygon(rgg, zgg, eq.rbbbs, eq.zbbbs);
  psin_grid = (eq.psizr - eq.psimag) / (eq.psibry - eq.psimag);
  psin_grid(~in) = nan;

  % load the profiles to be scaled
  % notation: zero suffix is the unscaled quantity    
  np = 1;
  nf = 1; 
  plotit = 0;
  if isfield(opts, 'time')
    trange = opts.time + [-0.04 0.04];
  else
    trange = [0 inf];
  end
  p = profile_test2(np, nf, plotit, trange);
  pprime0 = p.pprime;
  
  ffprim0 = p.ffprim_mean;
  % ffprim0 = p.ffprim;  
  % [pprime0, ffprim0] = load_standard_efit_profiles;

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
  rbbbs = eq.rbbbs;
  zbbbs = eq.zbbbs;
  
  pla = variables2struct(psizr_pla, pcurrt, rbbbs, zbbbs, jphi, ffprim, pprime, pres);

end












































