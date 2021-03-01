


file_args = {Mxx, Rxx, fit_coils, circ, Rext_mOhm, Lext_mH, enforce_stability};
parameters = {'Mvv', Mvv; 'Mvc', Mvc; 'Mvp', Mvp; 'Rvv_mOhm', Rvv0; 'Lvv', Lvv0; 'Mcp' Mcp};
odefun = 'coil_plus_vessel_dynamics';
sys = idgrey(odefun, parameters, 'd', file_args, Ts, 'InputDelay', 3);

tsim = 0:0.01:1;
load('xsim.mat')

ic = timeseries(xsim(:,circ.iicx), tsim);
ip = timeseries(xsim(:,circ.iipx), tsim);
icts = resample(icts,tsample);   
ipts = resample(ipts,tsample);

Tsmooth = 10;  % [ms]
nsmooth = floor(Tsmooth/1000/Ts);

ic = smoothdata(icts.Data, 1, 'movmean', nsmooth);
ip = smoothdata(ipts.Data, 1, 'movmean', nsmooth);

icdot = gradient(ic', Ts)';
ipdot = gradient(ip', Ts)';

icdot = smoothdata(icdot,1,'movmean',nsmooth);
ipdot = smoothdata(ipdot,1,'movmean',nsmooth);

u = double([vts.Data icdot ipdot]);

shotdata = iddata(y, u, Ts);







[yest,t,xest] = lsim(sys_est, u, tsample, x0);

figure
hold on
plot(t, xest(:,circ.iicx), 'b')
plot(t, y(:,circ.iicx), '--r')
title([num2str(shot) ' Coil Currents'], 'fontsize', 14)
ylabel('[A]', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
mylegend({'True', 'Simulated'}, {'--','-'}, [], {'r','b'}, [], 'Northeast', 14);

load('coils_greybox.mat')
figure
hold on
plot(t,xest(:,circ.iivx),'b')
plot(coils.t, coils.iv,'--r')
xlim([0 0.9])
title([num2str(shot) ' Vessel Currents'], 'fontsize', 14)
ylabel('[A]', 'fontsize', 14)
xlabel('Time [s]', 'fontsize', 14)
mylegend({'True', 'Simulated'}, {'--','-'}, [], {'r','b'}, [], 'Northeast', 14);
















