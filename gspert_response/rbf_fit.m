clear; clc; close all

isample = 80;
icoil = 1;
npts = 15;

targs = load('/Users/jwai/Research/rampup_nstxu/gspert_response/old/train_response_203008.mat').targs;
load('nstxu_obj_config2016_6565.mat')
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;

dpsidix = reshape(targs.dpsidix(isample,:,icoil), 65, 65);
psizr = squeeze(targs.psirz(isample,:,:))';

xg = rg;
yg = zg;
dx = mean(diff(xg));
dy = mean(diff(yg));
[xgg, ygg] = meshgrid(xg, yg);


nx = length(xg);
ny = length(yg);

xc = linspace(min(xg), max(xg), npts);
yc = linspace(min(yg), max(yg), npts);
[xc, yc] = meshgrid(xc, yc);
nc = numel(xc);  

terms = zeros(nx*ny,nc);
for i = 1:nc
  dist = sqrt((xc(i)-xgg(:)).^2 + (yc(i)-ygg(:)).^2);  
  terms(:,i) = rbf(dist);
end

z = dpsidix;
%%
coeff = terms \ z(:);

coeff = coeff + 0.3*median(abs(coeff))*randn(size(coeff));

zfit = terms*coeff;
zfit = reshape(zfit, ny, nx);


ix = floor(linspace(1,nx,npts));
iy = floor(linspace(1,ny,npts));
z_interp = griddata(xgg(iy,ix), ygg(iy,ix), z(iy,ix), xgg, ygg);



figure
subplot(121)
hold on
[~,cs] = contour(xg, yg, z, 20, '--r');
% contour(xg, yg, z_interp, cs.LevelList, 'b')
contour(xg(ix), yg(iy), z(iy,ix), cs.LevelList, 'b')

subplot(122)
hold on
[~,cs] = contour(xg, yg, z, 20, '--r');
contour(xg, yg, zfit, cs.LevelList, '-b')
set(gcf, 'Position', [700 658 609 554])


function out = rbf(dist, eps)
  if ~exist('eps','var'), eps=1; end
  out = exp(-(eps*dist).*2);
end
  















































