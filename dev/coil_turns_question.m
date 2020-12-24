ccc
load('eq204660_300.mat')
load('nstxu_obj_config2016_6565.mat')
struct_to_ws(tok_data_struct);
struct_to_ws(eq);

plot_nstxu_geo(tok_data_struct)
contour(rg,zg,psizr_app,100)
set(gcf,'position',[1047 16 365 655])

[r,z] = ginput;

[~, dpsidr, dpsidz] = bicubicHermite(rg, zg, psizr_app, r(1), z(1));
Br = -1/(2*pi*r(1))*dpsidz;
Bz =  1/(2*pi*r(1))*dpsidr;
B = sqrt(Br^2 + Bz^2);

d = norm([r(2) - r(1), z(2) - z(1)]);
mu0 = 4*pi*1e-7;


cccirc = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 11 12 13];
ccnames = {'OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', ...
  'PF5', 'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL'};

coil = 'PF5';
icirc = find(strcmp(cellstr(ccnames), coil));
icoils = find(cccirc == icirc);

I = 2*pi*d*B / mu0
eq.ic(icoils(1))
sum(eq.ic(icoils) .* ccnturn(icoils))







