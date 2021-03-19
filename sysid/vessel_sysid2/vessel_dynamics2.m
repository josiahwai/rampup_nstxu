% xdot = A*x + B*u
% x = iv(ivess), u = [icdot_true ivdot_true ipdot_true]

function [Ad, Bd, C, D] = vessel_dynamics2(Rv, Mvc, Mvv, Mvp, Ts, enforce_stability)

Ac = -inv(Mvv)*diag(Rv);

if enforce_stability
  Ac = numerically_stabilize(Ac, 1e3); 
end

Bc = -inv(Mvv) * [Mvc Mvp];

[Ad,Bd] = c2d(Ac,Bc,Ts);

C = eye(size(Ad));

D = zeros(size(C,1),size(Bd,2));







































