% xdot = A*x + B*u
% x = [ic iv]', u = [ps_voltages icdot_true ipdot_true]'

function [Ad, Bd, C, D] = vessel_dynamics3(Mvv, Mvc, Mvp, Rvv_mOhm, Lvv, Mcp, Mvc_scale, ...
  Ts, Mxx, Rxx, fit_coils, circ, Rext_mOhm, Lext_mH, enforce_stability)

Rvv = Rvv_mOhm / 1000;   


% if min(size(Rxx))==1, Rxx = diag(Rxx); end 
% 
% % inject the fitted Lext and Rext
% for i=1:length(fit_coils)
%     Mxx(fit_coils(i),fit_coils(i)) = Mxx(fit_coils(i),fit_coils(i)) + Lext_mH(i)/1000;
%     Rxx(fit_coils(i),fit_coils(i)) = Rxx(fit_coils(i),fit_coils(i)) + Rext_mOhm(i)/1000;
% end
% 
% Mvc = diag(Mvc_scale) * Mvc;
% 
% % inject the Mvc estimate
% Mxx(circ.iicx, circ.iivx) = Mvc';
% Mxx(circ.iivx, circ.iicx) = Mvc;
% 
% % inject the Rvv estimate
% Rxx(circ.iivx, circ.iivx) = diag(Rvv);
% 
% 
% 
% % Find A and B matrix portions corresponding to coils
% A = -inv(Mxx)*Rxx;
% Acoil = A(circ.iicx,:);
% 
% b1 = inv(Mxx) * [eye(circ.ncx); zeros(circ.nvx, circ.ncx)];
% b2 = zeros(circ.ncx+circ.nvx, circ.ncx);
% b3 = inv(Mxx) * [-Mcp; -Mvp];
% B = [b1 b2 b3];
% Bcoil = B(circ.iicx,:);

% Find A and B matrix portions corresponding to vessel
Avess = [ -inv(Mvv)*diag(Rvv)];
Bvess = [ -inv(Mvv)*[Mvc Mvp] ];

% % Merge, stabilize, and discretize
% Ac = [Acoil; Avess];
% Bc = [Bcoil; Bvess];
% 
if enforce_stability
  Avess = numerically_stabilize(Avess, 1e3); 
end

% Ad = inv(eye(size(Ac))-Ts*Ac);
% Bd = Ad*Bc*Ts;
[Ad,Bd] = c2d(Avess,Bvess,Ts);

C = eye(size(Ad));
D = zeros(size(C,1),size(Bd,2));







































