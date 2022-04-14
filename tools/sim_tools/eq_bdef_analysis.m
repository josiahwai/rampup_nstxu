% Get information on the boundary-definition: whether limited or not, 
% touch points, x-points
function boundary = eq_bdef_analysis(eq, tok_data_struct, opts)

if ~exist('opts','var'), opts = struct; end
if ~isfield(opts, 'plotit'), opts.plotit = 1; end

i = eq.rbbbs == 0 & eq.zbbbs == 0;
eq.rbbbs(i) = [];
eq.zbbbs(i) = [];
struct_to_ws(eq);

xlim = tok_data_struct.limdata(2,:);
ylim = tok_data_struct.limdata(1,:);
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;
rgg = tok_data_struct.rgg;
zgg = tok_data_struct.zgg;

% limited or diverted
[rzbbbs_fine] = contourc(rg,zg,psizr,[psibry psibry]);
rbbbs_fine = rzbbbs_fine(1,:);
zbbbs_fine = rzbbbs_fine(2,:);

e = 0.002;

i = rbbbs_fine < min(rbbbs) - e | rbbbs_fine > max(rbbbs) + e | ...
    zbbbs_fine < min(zbbbs) - e | zbbbs_fine > max(zbbbs) + e;

rbbbs_fine(i) = [];
zbbbs_fine(i) = [];

[xy,dist] = distance2curve([xlim(:) ylim(:)], [rbbbs_fine(:) zbbbs_fine(:)]);
  
[mindist, i] = min(dist);
islimited = double(mindist < .005);
if islimited
  rtouch = xy(i,1);
  ztouch = xy(i,2);
  psitouch = bicubicHermite(rg, zg, psizr, rtouch, ztouch);
else
  rtouch = nan;
  ztouch = nan;
  psitouch = nan;
end

if isempty(islimited), islimited = nan; end
  
% x-points
try
    
  % find best upper and lower guesses for x-point (pts on boundary where
  % grad(psi) closest to zero)
  iup = zbbbs > (max(zbbbs) + min(zbbbs))/2;
  ilo = find(~iup);
  iup = find(iup);
  [~,psi_r,psi_z] = bicubicHermite(rg,zg,psizr,rbbbs,zbbbs);

  [~,i] = min( psi_r(iup).^2 + psi_z(iup).^2);
  [~,j] = min( psi_r(ilo).^2 + psi_z(ilo).^2);  

  rx_up = rbbbs(iup(i));
  zx_up = zbbbs(iup(i));

  rx_lo = rbbbs(ilo(j));
  zx_lo = zbbbs(ilo(j));
  
  % zoom in on x-point
  [rx_up,zx_up,psix_up] = isoflux_xpFinder(psizr,rx_up,zx_up,rg,zg);
  [rx_lo,zx_lo,psix_lo] = isoflux_xpFinder(psizr,rx_lo,zx_lo,rg,zg);

  % check if inside limiter
  if ~all(inpolygon([rx_lo rx_up], [zx_lo zx_up], xlim, ylim))
    rx_lo = nan;
    zx_lo = nan;
    rx_up = nan;
    zx_up = nan;
  end
  
catch
  rx_lo = nan;
  zx_lo = nan;
  rx_up = nan;
  zx_up = nan;
end

if islimited
  rbdef = rtouch;
  zbdef = ztouch;
  psibry = psitouch;
elseif abs(psix_lo - psibry) < abs(psix_up - psibry)
  rbdef = rx_lo;
  zbdef = zx_lo;
  psibry = psix_lo;
else
  rbdef = rx_up;
  zbdef = zx_up;
  psibry = psix_up;
end

if opts.plotit
  plot_nstxu_geo(tok_data_struct)
  scatter(rx_lo, zx_lo, 80, 'bx')
  scatter(rx_up, zx_up, 80, 'bx') 
  scatter(rtouch, ztouch, 80, 'b', 'filled')
  scatter(rbdef, zbdef, 200, 'r', 'linewidth', 2)
  contour(rg,zg,psizr,[psibry psibry])
  drawnow
end


% some geometric parameters
R0 = (max(rbbbs) + min(rbbbs)) / 2;
% Z0 = (median(zbbbs(zbbbs>0)) + median(zbbbs(zbbbs<0)))/2;
Z0 = (max(zbbbs) + min(zbbbs)) / 2;

a = (max(rbbbs) - min(rbbbs)) / 2;
epsilon = a / R0;
aspect_ratio = 1/epsilon;
kappa = (max(zbbbs) - min(zbbbs)) / (2*a);

[~,i] = max(zbbbs);
tritop = (R0 - rbbbs(i)) / a;

[~,i] = min(zbbbs);
tribot = (R0 - rbbbs(i)) / a;

tri = (tritop + tribot) / 2;

[rout, iout] = max(rbbbs);
zout = zbbbs(iout);

[rin, iin] = min(rbbbs);
zin = zbbbs(iin);

[ztop, itop] = max(zbbbs);
rtop = rbbbs(itop);

[zbot, ibot] = min(zbbbs);
rbot = rbbbs(ibot);

xi_ou = measure_squareness(rout, rtop, zout, ztop, rbbbs, zbbbs, false);
xi_iu = measure_squareness(rin, rtop, zin, ztop, rbbbs, zbbbs, false);
xi_ol = measure_squareness(rout, rbot, zout, zbot, rbbbs, zbbbs, false);
xi_il = measure_squareness(rin, rbot, zin, zbot, rbbbs, zbbbs, false);


boundary = variables2struct(islimited, psibry, rbdef, zbdef, rtouch, ztouch, psitouch, ...
  rx_lo, zx_lo, psix_lo, rx_up, zx_up, psix_up, R0, Z0, a, epsilon, aspect_ratio, ...
  kappa, tritop, tribot, tri, xi_ou, xi_iu, xi_ol, xi_il, rbbbs, zbbbs);












