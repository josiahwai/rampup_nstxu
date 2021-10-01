

function eq = eq_analysis(psizr, pla, tok_data_struct, opts)

% ccc
% load('matlab.mat')

if ~exist('opts','var'), opts = struct; end
if ~isfield(opts, 'plotit'), opts.plotit = 1; end

xlim = tok_data_struct.limdata(2,:)';
ylim = tok_data_struct.limdata(1,:)';
[xlim, ylim] = interparc(xlim, ylim, 500, true, true);
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;

rbbbs0 = pla.rbbbs; % previous iteration/estimate of boundary
zbbbs0 = pla.zbbbs;

P = polyshape(rbbbs0, zbbbs0);
[rmaxis, zmaxis] = centroid(P);
[rmaxis, zmaxis, psimag] = isoflux_xpFinder(psizr, rmaxis, zmaxis, rg, zg);

psibry0 = bicubicHermite(rg, zg, psizr, min(xlim), 0);

[xy_touch,dist] = distance2curve([xlim ylim], [rbbbs0' zbbbs0']);

[~,i] = sort(dist);
i = i(1:floor(end/3));

psilim_touch = bicubicHermite(rg, zg, psizr, xy_touch(i,1), xy_touch(i,2));
psibarlim_touch = (psilim_touch - psimag) / (psibry0 - psimag);
[psibar_touch,i] = min(psibarlim_touch);
psi_touch = psilim_touch(i); 
r_touch = xy_touch(i,1);
z_touch = xy_touch(i,2);


% x-points
try
    
  % find best upper and lower guesses for x-point (pts on boundary where
  % grad(psi) closest to zero)
  iup = zbbbs0 > (max(zbbbs0) + min(zbbbs0))/2;
  ilo = find(~iup);
  iup = find(iup);
  [~,psi_r,psi_z] = bicubicHermite(rg,zg,psizr,rbbbs0,zbbbs0);

  [~,i] = min( psi_r(iup).^2 + psi_z(iup).^2);
  [~,j] = min( psi_r(ilo).^2 + psi_z(ilo).^2);  

  rx_up = rbbbs0(iup(i));
  zx_up = zbbbs0(iup(i));

  rx_lo = rbbbs0(ilo(j));
  zx_lo = zbbbs0(ilo(j));
  
  % zoom in on x-point
  [rx_up,zx_up,psix_up] = isoflux_xpFinder(psizr,rx_up,zx_up,rg,zg);
  [rx_lo,zx_lo,psix_lo] = isoflux_xpFinder(psizr,rx_lo,zx_lo,rg,zg);

  % check if inside limiter
  if ~inpolygon(rx_lo, zx_lo, xlim, ylim)
    rx_lo = nan;
    zx_lo = nan;
  end
  if ~inpolygon(rx_up, zx_up, xlim, ylim)
    rx_up = nan;
    zx_up = nan;
  end
  
catch
  rx_lo = nan;
  zx_lo = nan;
  rx_up = nan;
  zx_up = nan;
end

psibar_xlo = (psix_lo - psimag) / (psibry0 - psimag);
psibar_xup = (psix_up - psimag) / (psibry0 - psimag);

[~,i] = nanmin([psibar_touch psibar_xlo psibar_xup]);

if i == 1
  psibry = psi_touch;
  rbdef = r_touch; 
  zbdef = z_touch;
  islimited = 1;
elseif i == 2
  psibry = psix_lo;
  rbdef = rx_lo;
  zbdef = zx_lo;
  islimited = 0;
elseif i == 3
  psibry = psix_up;
  rbdef = rx_up;
  zbdef = zx_up;
  islimited = 0;
end

% trace boundary
if opts.robust_tracing
  [~,~,~, rbbbs, zbbbs] = trace_contour(...
    rg,zg,psizr,rbdef,zbdef, rmaxis, zmaxis,xlim,ylim,0,0);
  rbbbs = rbbbs{1};
  zbbbs = zbbbs{1};
else
  crz = contourc(rg, zg, psizr, [psibry psibry]);
  r = crz(1,:);
  z = crz(2,:);
  i = z > zx_lo & z < zx_up;
  rbbbs = r(i);
  zbbbs = z(i);
end


eq = variables2struct(psizr, psimag, psibry, rbdef, zbdef, rbbbs, zbbbs, islimited, ... 
  r_touch, z_touch, rx_lo, zx_lo, rx_up, zx_up);


%%

if opts.plotit
	figure
  hold on
  contour(rg, zg, psizr, 30)
  % contour(rg, zg, psizr, [psibry psibry], 'linewidth', 3)
  plot(rbbbs, zbbbs, 'r', 'linewidth', 3)
  scatter(rmaxis, zmaxis, 'filled')
  scatter(rx_lo, zx_lo, 100, 'x', 'linewidth', 3)
  scatter(rx_up, zx_up, 100, 'x', 'linewidth', 3)
  scatter(r_touch, z_touch, 100, 'o', 'linewidth', 3)
  axis equal
  set(gcf, 'Position', [76 118 466 649])
end
























































