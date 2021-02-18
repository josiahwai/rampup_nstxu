function [x,y] = greymodel_oneattime(t, x, u, rcoil, lcoil, m_plasma_fac, voltage_fac, coil_voltage, file_args)

Ts = file_args{1};
circ = file_args{2};
lstar_invs = file_args{3};
timebase = file_args{4};
traj = file_args{5};
icoil_to_estimate = file_args{6};
rxx = file_args{7};


[~,k] = min(abs(timebase-t));
xtrue = traj.x(k,:)';


rxx(icoil_to_estimate) = rcoil;
lstari = squeeze(lstar_invs(k,:,:));

Ac = -lstari * diag(rxx);
Bc = lstari(:,1:circ.ncx);

Ac(circ.iicx_remove,:) = 0;
Ac(:,circ.iicx_remove) = 0;
Bc(circ.iicx_remove,:) = 0;
Bc(:,circ.iicx_remove) = 0;

Astab = removeFirstEigval(Ac);

[A, B] = c2d(Astab, Bc, Ts);

% u(icoil_to_estimate) = u(icoil_to_estimate) * voltage_fac;
u(icoil_to_estimate) = coil_voltage(k);

xmod = xtrue;
xmod(icoil_to_estimate) = x(icoil_to_estimate);

x = A*xmod + B*u(:);
y = x(icoil_to_estimate);



















