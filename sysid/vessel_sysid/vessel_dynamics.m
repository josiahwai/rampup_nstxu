% xdot = A*x + B*u
% x = iv(ivess), u = [icdot_true ivdot_true ipdot_true]

function [Ad, Bd, C, D] = vessel_dynamics(Rv, Lv, Mvc, Mvv, Mvp, Ts, ivess)

Ac = -Rv / Lv;

Mvv(ivess) = 0;

Bc = - [Mvc Mvv Mvp] / Lv;

[Ad,Bd] = c2d(Ac,Bc,Ts);

C = eye(size(Ad));

D = zeros(size(C,1),size(Bd,2));







































