% INPUTS:
% bry - struct with fields 'r' and 'z', position of target boundary points
% ip - target plasma current
% wmhd - target plasma thermal energy 
% tok_data_struct - struct with machine geometry
% eq - optional struct with 
% opts - optional struct with fields: plotit - flag to make plots


% OUTPUT
% pla - struct containing estimates of the plasma current and flux distribution
%       an estimate 

function pla = estimate_pla(target, tok_data_struct, eq, opts)
% ccc
% load('matlab.mat')
% clear eq

if ~exist('eq','var'), eq = []; end
if ~exist('opts', 'var'), opts = struct; end
if ~isfield(opts, 'plotit'), opts.plotit = 0; end
if ~isfield(opts, 'first_iteration'), opts.first_iteration = 0; end

struct_to_ws(tok_data_struct);
mpc = tok_data_struct.mpc; 
mpp = mpp_full;
dr = mean(diff(rg));
dz = mean(diff(zg));
dA = dr * dz;

%%

% ==========================================================
% Estimate current distribution without knowing applied flux
% ==========================================================
% On the first iteration, need to estimate the plasma flux without knowing
% the previous iteration info (boundary, external flux, etc). Use the
% target boundary and target Ip to find a plasma current distribution: 
% j = jhat * (1 - x^a)^b, where x is the distance ratio from centroid to
% boundary (0 at centroid, 1 at boundary), a and b are constants, and jhat
% is a constant scaled to match the target Ip. 

if isempty(eq) || opts.first_iteration % no previous iteration info
  
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

  psizr_pla = mpp * pcurrt(:);
  psizr_pla = reshape(psizr_pla, nr, nz);

else

  %%
  % ===================================
  % Iterate plasma current distribution 
  % ===================================
  % Use standard P' and FF' profiles to estimate the current distribution. P'
  % and FF' will be scaled in order to match the target Wmhd and target Ip. 

  
  psi = eq.psizr;
  %   psi_app = reshape(mpc * eq.ic + mpv*eq.iv, nr, nz);
  %   psi_pla = psi - psi_app;

  mu0 = pi*4e-7;
  psin = linspace(0,1,nr);

  % find boundary
  rz = solveTSP([target.rcp(:) target.zcp(:)], 0); % sort via traveling salesman
  rz = [rz; rz(1,:)];
  bry = fnplt(cscvn(rz'));  % spline interpolation
  rbbbs = bry(1,:);
  zbbbs = bry(2,:);
  % rbbbs = eq.rbbbs;
  % zbbbs = eq.zbbbs;

  in = inpolygon(rgg, zgg, eq.rbbbs, eq.zbbbs);
  psin_grid = (psi - eq.psimag) / (eq.psibry - eq.psimag);
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

  % make some plots
  if opts.plotit
    figure
    hold on
    plot(pprime)
    plot(eq.pprime)
    legend('Estimate', 'EFIT', 'fontsize', 18)

    figure
    hold on
    plot(pres)
    plot(eq.pres)
    legend('Estimate', 'EFIT', 'fontsize', 18)
  end


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
  
  % make some plots
  if opts.plotit
    figure
    hold on
    plot(ffprim)
    plot(eq.ffprim)  
    legend('Estimate', 'EFIT', 'fontsize', 18)

    figure
    contourf(jphi)
    colorbar
    title('Estimate')

    figure
    contourf(eq.jphi)
    colorbar
    title('EFIT')
  end

end
%%

% write to output
psizr_pla = mpp * pcurrt(:);
psizr_pla = reshape(psizr_pla, nr, nz);

pla = variables2struct(psizr_pla, pcurrt, jphi);

% make some plots
if opts.plotit
  figure
  contourf(rg,zg,jphi, 10)
  colorbar
  title('Estimate')
   
  figure
  subplot(121)
  [~,cs] = contourf(rg, zg, psizr_pla, 20);
  colorbar  
  caxis([min(cs.LevelList) max(cs.LevelList)])
  axis equal
  try 
    subplot(122)
    psi_app = reshape(mpc * eq.ic + mpv*eq.iv, nr, nz);
    contourf(rg, zg, eq.psizr - psi_app, cs.LevelList)
    colorbar
    axis equal
    caxis([min(cs.LevelList) max(cs.LevelList)])
  catch
  end
 
  
  drawnow
end



















