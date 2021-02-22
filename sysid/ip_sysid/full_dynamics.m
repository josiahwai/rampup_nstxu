function [x, Ad, Bd, C, D] = full_dynamics(x, u, Mxx, Rxx, Rp, Lp, Mpc, Mpv, Ts, circ, enforce_stability)

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

x = Ad*x + Bd * u;



ii_unipolar = [2 5 10 13];
x(ii_unipolar) = max(x(ii_unipolar), 0);











































