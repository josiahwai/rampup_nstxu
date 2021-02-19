function [Ad, Bd, C, D] = coil_and_vessel_dynamics(Mxx, Rxx, fit_coils, ...
  Rext_mOhm, Lext_mH, Mvc, Rvv, Ts, circ, Mcp, Mvp)

nx = length(Rxx);

% Plug in the fit Lext and Rext
Rxx = diag(Rxx);
for i=1:length(fit_coils)
    Mxx(fit_coils(i),fit_coils(i)) = Mxx(fit_coils(i),fit_coils(i)) + Lext_mH(i)/1000;
    Rxx(fit_coils(i),fit_coils(i)) = Rxx(fit_coils(i),fit_coils(i)) + Rext_mOhm(i)/1000;
end

% plug in the fit resistances
rxx = diag(Rxx);
rxx(circ.iivx) = Rvv/1000;
Rxx = diag(rxx);

Mcc = Mxx(circ.iicx, circ.iicx);
Mvv = Mxx(circ.iivx, circ.iivx);

Rc = Rxx(circ.iicx, circ.iicx);
Rv = Rxx(circ.iivx, circ.iivx);

% Coils
Mvc_old = Mxx(circ.iivx, circ.iicx);
Acoil = -inv(Mcc) * Rc;
Bcoil = -inv(Mcc) * [eye(circ.ncx) zeros(circ.ncx) Mvc_old' Mcp];

% Vessel
Avess = -inv(Mvv) * Rv;
Bvess = -inv(Mvv) * [zeros(circ.nvx, circ.ncx) Mvc zeros(circ.nvx) Mvp];

% Merged
A = blkdiag(Acoil, Avess);
B = [Bcoil; Bvess];

% Discretize
Ad = inv(eye(size(A))-Ts*A);
Bd = Ad*B*Ts;
C = eye(size(A));
D = zeros(size(C,1),size(Bd,2));


% 
% % plug in the fit Mvc
% Mxx(circ.iicx, circ.iivx) = Mvc';
% Mxx(circ.iivx, circ.iicx) = Mvc;
% 
% % Form the continuous-time matrices
% invM = inv(Mxx);
% A = -invM * Rxx;
% B1 = invM(:,circ.iicx);
% B2 = -invM*[Mcp; Mvp];
% B = [B1 B2];
% 
% A = numerically_stabilize(A, 1e5);




% A0 = -inv(Mxx) * Rxx;
% B0 = -inv(Mxx);
% B0 = B0(:,circ.iicx);
% 
% Ac = A0(circ.iicx,:);
% Ap = A0(circ.iipx,:);
% Avv = -inv(Mvv)*diag(Rvv/1000);
% Avc = zeros(circ.nvx,circ.ncx);
% Avp = zeros(circ.nvx,1);
% A = [Ac; Avc Avv Avp; Ap];
% 
% Bc1 = B(circ.iicx,:);
% Bp1 = B(circ.iipx,:);
% Bv1 = zeros(circ.nvx,circ.ncx);
% Bc2 = zeros(circ.ncx, circ.ncx + circ.np);
% Bp2 = zeros(1, circ.ncx + circ.np);
% Bv2 = -inv(Mvv)*[Mvc Mvp];
% B = [Bc1 Bc2; Bv1 Bv2; Bp1 Bp2];
% 
% % Discretize
% Ad = inv(eye(circ.nx)-Ts*Ac);
% Bd = Ad*B*Ts;
% C = eye(circ.nx);
% D = zeros(size(C,1),size(Bd,2));







