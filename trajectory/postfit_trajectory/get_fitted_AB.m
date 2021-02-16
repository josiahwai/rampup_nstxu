function [Ad,Bd] = get_fitted_AB(traj, circ, fitted_resistances, Ts, i)

struct_to_ws(fitted_resistances);

lstari = squeeze(traj.lstari(i,:,:));
rxx = [rcc; rvv; rp_t(i)];

Ac = -lstari * diag(rxx);
Bc = lstari(:,1:circ.ncx);

Ac(circ.iicx_remove,:) = 0;
Ac(:,circ.iicx_remove) = 0;
Bc(circ.iicx_remove,:) = 0;
Bc(:,circ.iicx_remove) = 0;

Astab = numerically_stabilize(Ac, 1e5);
[Ad, Bd] = c2d(Astab, Bc, Ts);

Bd = Bd * diag(voltage_scale);





