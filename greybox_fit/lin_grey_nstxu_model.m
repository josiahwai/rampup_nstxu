function [A,B,C,D] = lin_greybox_model_res(rcc, rvv, rp, lp, Ts, circ, sys)


rxx = [rcc; rvv; rp];
lstar = sys.lstar;
lstar(end,end) = lp;
lstari = inv(lstar);

Ac = -lstari * diag(rxx);
Bc = lstari(:,1:circ.ncx);

Ac(circ.iremove,:) = 0;
Ac(:,circ.iremove) = 0;
Bc(circ.iremove,:) = 0;
Bc(:,circ.iremove) = 0;

Astab = removeFirstEigval(Ac);

[A, B] = c2d(Astab, Bc, Ts);
C = circ.Pxx_keep;
D = zeros(size(C,1), size(B,2));



















