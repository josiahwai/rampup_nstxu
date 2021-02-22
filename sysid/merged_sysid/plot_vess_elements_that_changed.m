ccc
load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);


% find elements that had resistance change > thresh
load('merge_sys_est.mat')
thresh = 2; % mOhm
Rxx = diag( load('NSTXU_vaccum_system.mat').NSTXU_vacuum_system.Rxx);
Rvv0 = Rxx(circ.iivx) * 1000;
Rvv = sys_est.Structure.Parameters(4).Value;

dr = Rvv0 - Rvv;
k = find(abs(dr) > thresh);
[~,isort] = sort(abs(dr(k)), 'descend');
k = k(isort);

disp('Coils that changed the most:')
for j = 1:length(k)
  i = k(j);
  disp( [circ.vvnames{i} ': ' num2str(dr(i)) ' mOhm']);
end

vvdata = tok_data_struct.vvdata;

dum = zeros(circ.nvx, 1);
dum(k) = 1;
z = circ.Pvv * dum;
ivess = find(z ~= 0);

rvess = vvdata(2,ivess);
zvess = vvdata(1,ivess);


figure
vvlabels = categorical(circ.vvnames);
bar(vvlabels, [Rvv0 Rvv])
ylabel('Resistance [mOhm]')
legend('Original', 'Fit', 'fontsize', 14)

figure
opt.icolor = 0;
plot_nstxu_geo(tok_data_struct, opt)
scatter(rvess, zvess, 100, 'filled', 'r')




































