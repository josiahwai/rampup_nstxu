clear; clc; close all

isample = 50;
icoil = 1;
targs = load('/Users/jwai/Research/rampup_nstxu/gspert_response/old/train_response_203008.mat').targs;
load('nstxu_obj_config2016_6565.mat')
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;

dpsidix = reshape(targs.dpsidix(isample,:,icoil), 65, 65);
psizr = squeeze(targs.psirz(isample,:,:))';




degree_x = 30;
degree_y = 30;
xg = rg;
yg = zg;


nx = length(xg);
ny = length(yg);

xgn = xg;
ygn = yg;

xmax = max(xg);
xmin = min(xg);
xgn = 2*(xg - (xmin+xmax)/2) / (xmax-xmin);

ymax = max(yg);
ymin = min(yg);
ygn = 2*(yg - (ymin+ymax)/2) / (ymax-ymin);

cheby_x = {};
for i = 1:degree_x
  cheby_x{i} = chebyshevT(i, xgn(:));
end

cheby_y = {};
for j = 1:degree_y
  cheby_y{j} = chebyshevT(j, ygn(:));
end

terms = zeros(nx, ny, degree_x, degree_y);
for i = 1:degree_x
  for j = 1:degree_y
    terms(:,:,i,j) = cheby_x{i} * cheby_y{j}';
  end
end

terms = reshape(terms, nx*ny, degree_x*degree_y);

z = dpsidix;

coeff = terms \ z(:);

zfit = terms*coeff;
zfit = reshape(zfit, ny, nx);



figure
hold on
[~,cs] = contour(rg, zg, z, 20, '--r');
contour(rg, zg, zfit, cs.LevelList, '-b')
set(gcf, 'Position', [907 834 305 463])

















































