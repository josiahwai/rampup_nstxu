function [x, Ad, Bd, C, D, iv] = full_dynamics_vessmod(x, u, Mxx, Rxx, Rp, Lp, Mpc, Mpv,...
  Ts, circ, enforce_stability, iv, icdot, ipdot)

Mpcv = [Mpc Mpv];
Mss = [Mxx Mpcv'; Mpcv Lp];
Rss = [Rxx; Rp];

invMss = inv(Mss);
Ac = -invMss*diag(Rss);
Bc = invMss(:,circ.iicx);

if enforce_stability
  Ac = numerically_stabilize(Ac, 1e3); 
end


[Ad,Bd] = c2d(Ac,Bc,Ts);

C = eye(size(Ad));
D = zeros(size(C,1),size(Bd,2));

xprev = x;
x = Ad*x + Bd * u;



ii_unipolar = [2 5 10 13];
x(ii_unipolar) = max(x(ii_unipolar), 0);






% Find A and B matrix portions corresponding to vessel
Mvv = Mxx(circ.iivx, circ.iivx);
Mvc = Mxx(circ.iivx, circ.iicx);
Rvv = Rxx(circ.iivx);

Avess = -inv(Mvv)*diag(Rvv);
Bvess = -inv(Mvv)*[Mvc Mpv'];

[Avess_d, Bvess_d] = c2d(Avess, Bvess, Ts);

dx = x - xprev; 
icdot = dx(circ.iicx) / Ts;
ipdot = dx(circ.iipx) / Ts;

iv = Avess_d * iv + Bvess_d * [icdot; ipdot];


















































