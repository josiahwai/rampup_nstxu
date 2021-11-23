


% [x,y] = vessel_dynamics(t, x, u, rc, rv, lc, lv, mcc_triu, mvv_triu, mcv, file_args)


function [Ad, Bd, C, D] = coil_vessel_dynamics(Mvv, Mvc, Mvp, Rvv_mOhm, Lvv, Mcp, Mvc_scale, ...
  Ts, Mxx, Rxx, fit_coils, circ, Rext_mOhm, Lext_mH, enforce_stability)




function [Ad, Bd, C, D] = coil_vessel_dynamics(x, u, rc, rv, lc, lv, mcc_triu, mvv_triu, mcv, file_args)
% 
% args = file_args{1};
% circ = args.circ;
% ts = args.ts;
% 
% % extract components from x
% ic = x(circ.iicx);
% iv = x(circ.iivx);
% 
% % reconstruct mutual inductances from parameters
% mcc = mcc_triu + mcc_triu' + diag(lc);
% mvv = mvv_triu + mvv_triu' + diag(lv);
% M = [mcc mcv; mcv' mvv];
% R = diag([rc; rv]);
% 
% Minv = inv(M);
% 
% A = -Minv*R;
% B = [Minv(:,circ.iicx) Minv];
% 
% A = numerically_stabilize(A, 100);
% 
% % [Ad,Bd] = c2d(A,B,ts);
% 
% Ad = eye(size(A)) + A*ts;
% Bd = B*ts;
% 
% x = Ad*x(:) + Bd*u(:);
% 
% i = x(circ.ii_unipolar) < 0;
% x(circ.ii_unipolar(i)) = 0;
% 
% y = x;














% Rvv = Rvv_mOhm / 1000;   
% 
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
% % Ac = -inv(Mxx)*Rxx;
% % Bc = inv(Mxx) * [eye(circ.ncx); zeros(circ.nvx, circ.ncx)];
% % 
% % Ad = inv(eye(size(Ac))-Ts*Ac);
% % Bd = Ad*Bc*Ts;
% % 
% % C = eye(size(Ac));
% % D = zeros(size(C,1),size(Bd,2));
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
% 
% % Find A and B matrix portions corresponding to vessel
% Avess = [zeros(circ.nvx, circ.ncx) -inv(Mvv)*diag(Rvv)];
% Bvess = [zeros(circ.nvx, circ.ncx) -inv(Mvv)*[Mvc Mvp] ];
% 
% % Merge, stabilize, and discretize
% Ac = [Acoil; Avess];
% Bc = [Bcoil; Bvess];
% 
% if enforce_stability
%   Ac = numerically_stabilize(Ac, 1e3); 
% end
% 
% % Ad = inv(eye(size(Ac))-Ts*Ac);
% % Bd = Ad*Bc*Ts;
% [Ad,Bd] = c2d(Ac,Bc,Ts);
% 
% C = eye(size(Ad));
% D = zeros(size(C,1),size(Bd,2));







































