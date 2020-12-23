load('fit_Rp.mat')

t1 = linspace(0,0.4);
rp1 = fit_Rp(t1)';

t2 = linspace(0.4,2);
rp2 = rp1(end) * ones(size(t2));

rp = [rp1 rp2];
t = [t1 t2];

cftool(t,rp)












