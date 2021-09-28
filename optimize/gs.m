% first pass: estimate psin from the target boundary, fit pprime and
%             ffprim to get jphi distribution
%
% subsequent: use psi at previous equilibrium to iterate profiles


% %%
% function psi_pla  = gs(
% 
% rcp, zcp, ip, wmhd, psi_app, psi_pla

%%
ccc

debug = 0;

idx = 10;
load('matlab.mat')

eq = eqs.gdata(idx);
fns = fieldnames(targets);
for i = 1:length(fns)
  targets.(fns{i}) = targets.(fns{i})(idx,:);
end

targets.Wmhd = 7.5532e+04; 

struct_to_ws(tok_data_struct);
dr = mean(diff(rg));
dz = mean(diff(zg));
dA = dr*dz;



% ----------------------------------------------------------
% estimate a psin distribution, assuming the target boundary
% ----------------------------------------------------------

% sort and interpolate target boundary vertices, for well-defined polygon
r = [targets.rcp targets.rbdef];
z = [targets.zcp targets.zbdef];
rz = solveTSP([r' z'], false); % sort via traveling salesman, Cl := contour length of bry
rz = [rz; rz(1,:)];
n = length(rz);
rbbbs = interp1(1:n, rz(:,1), linspace(1,n), 'spline');
zbbbs = interp1(1:n, rz(:,2), linspace(1,n), 'spline');

P = polyshape(rbbbs,zbbbs);


shafranov_shift = .08;
[~, zmaxis] = centroid(P);
rmaxis = (min(rbbbs)+max(rbbbs))/2 + shafranov_shift;


% % Find distance from each grid point to axis & boundary: 
% % Draw a line from magnetic axis to each grid point to the magnetic axis.
% % Extend the line. Find where the line intersects target bry. 
% vec = [rgg(:) - rmaxis, zgg(:) - zmaxis];
% norms = vecnorm(vec')';
% vecn = vec ./ repmat(norms,1,2);
% for i = 1:(nr*nz)  
%   curve = [rmaxis zmaxis] + [0 0; vecn(i,:) * 2];
%   [xy, dist] = distance2curve(curve, [rbbbs(:) zbbbs(:)]);
%   [~,j] = min(dist);  
%   xy = xy(j,:);
%   
%   % distance ratio gridpoint to axis versus boundary to axis
%   xbar(i,1) = norm([rgg(i) - rmaxis, zgg(i) - zmaxis]) / norm(xy - [rmaxis zmaxis]); 
% end
% 
% xbar = reshape(xbar, nr, nz);
% in = inpolygon(rgg, zgg, rbbbs, zbbbs);
% xbar(~in) = nan;

dist2cent = sqrt((rgg(:)-rmaxis).^2 + (zgg(:)-zmaxis).^2);
[xy,dist2bry] = distance2curve([rbbbs(:) zbbbs(:)], [rgg(:) zgg(:)]);
xbar = dist2cent ./ (dist2cent + dist2bry);
xbar = reshape(xbar, nr, nz);
in = inpolygon(rgg, zgg, rbbbs, zbbbs);
xbar(~in) = nan;

psin_grid = xbar.^1.8;

if debug
  figure
  contourf(psin_grid.^1.8)
  colorbar

  psi = (eq.psizr - eq.psimag) / (eq.psibry - eq.psimag);
  psi(~in) = nan;
  figure
  contourf(psi)
  colorbar
  set(gcf, 'Position', [452 872 560 420])
end

%%

% estimate psibry-psimag
mu0 = pi*4e-7;
Rmax = max(rbbbs);
a =  (max(rbbbs) + min(rbbbs))/2 - shafranov_shift;
ip = eq.cpasma;
Cl = sum(sqrt( (diff(rbbbs).^2 + diff(zbbbs).^2)));
psibry_minus_psimag = -a*Rmax*mu0*ip*pi/Cl; 

%%

psin = linspace(0,1,nr);

% scale P' profile to match Wmhd 
[pprime0, ffprim0] = load_standard_efit_profiles; % profiles to be scaled


pres0 = cumtrapz(psin, pprime0) * psibry_minus_psimag / (2*pi);  % pressure profile to be scaled
pres0 = pres0-pres0(end);
pres0_grid = interp1(psin, pres0, psin_grid(:));

tmp_wmhd = nansum(3 * pi * rgg(:) .* pres0_grid * dA);

pprime_coef = targets.Wmhd / tmp_wmhd;

pres = pprime_coef * pres0;
pprime = pprime_coef * pprime0;

pprime_grid = interp1(psin, pprime, psin_grid(:));
pprime_grid = reshape(pprime_grid, nr, nz);


if debug
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

%%

% scale FF' to match Ip

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

ip_ffprim = targets.ip - ip_pressure;

ffprim_coeff = ip_ffprim / ip0_ffprim;

ffprim = ffprim0 * ffprim_coeff;
ffprim_grid = ffprim0_grid * ffprim_coeff;

% current from scaled profiles
jphi = rgg .* pprime_grid + ffprim_grid ./ (rgg * mu0);
pcurrt = jphi * dr * dz;
pcurrt(isnan(pcurrt)) = 0;
  
  
 
if debug
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




%%
% pres = cumtrapz(psin, eq.pprime) * (eq.psibry-eq.psimag) / (2*pi);
% pres = pres-pres(end);
% pres_grid = interp1(psin, pres, psin_grid(:));
% 
% Wmhd = nansum(3 * pi * rgg(:) .* pres_grid * dr * dz)






%%



% figure
% hold on
% plot(P)
% plot([rmaxis rmaxis+vecn(i,1)*2], [zmaxis zmaxis+vecn(i,2)*2])
% scatter(xy(j,1), xy(j,2))
% 
% % at each grid point, find normalized distance between centroid and
% % boundary
% dr = mean(diff(rg));
% dz = mean(diff(zg));
% rgridc = rgg + dr/2;
% zgridc = zgg + dz/2;
% 
% dist2cent = sqrt((rgridc(:)-rc).^2 + (zgridc(:)-zc).^2);
% 
% [xy,dist2bry] = distance2curve([rbbbs(:) zbbbs(:)], [rgridc(:) zgridc(:)]);
% 
% % j =1000;
% % scatter(rgridc(j), zgridc(j), 100, 'filled')
% % scatter(xy(j,1), xy(j,2), 100, 'filled')
% % dist(j)
% 
% 
% x = dist2cent ./ (dist2cent + dist2bry);
% 
% in = inpolygon(rgridc(:), zgridc(:), P.Vertices(:,1), P.Vertices(:,2));
% x(~in) = nan;
% x = reshape(x, nr, nz);






























