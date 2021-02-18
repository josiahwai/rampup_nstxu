function [A, B, C, D] = vessel_dynamics(Mvv, Mvc, Mvp, Rvv, Lvv, Ts)

Mvv = Mvv - diag(diag(Mvv)) + diag(Lvv); 

Ac = -inv(Mvv)*diag(Rvv/1000);

Bc = -inv(Mvv)*[Mvc Mvp];

nv = length(Rvv);

A = inv(eye(nv)-Ts*Ac);

B = A*Bc*Ts;

C = eye(size(A));

D = zeros(size(C,1),size(B,2));

end
