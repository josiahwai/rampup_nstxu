% xdot = A*x + B*u
% x = [ic iv]
% u = [vc, icdot, ivdot, ipdot, mcIp, mvIp]
%
% Note that icdot, ivdot, ipdot are the experimental values. 


function [x,y] = vessel_dynamics(t, x, u, rc, rv, lc, lv, mcc_triu, mvv_triu, mcv, file_args)

args = file_args{1};
circ = args.circ;
ts = args.ts;

% extract components from x
ic = x(circ.iicx);
iv = x(circ.iivx);

% reconstruct mutual inductances from parameters
mcc = mcc_triu + mcc_triu' + diag(lc);
mvv = mvv_triu + mvv_triu' + diag(lv);
M = [mcc mcv; mcv' mvv];
R = diag([rc; rv]);

Minv = inv(M);

A = -Minv*R;
B = [Minv(:,circ.iicx) Minv];

A = numerically_stabilize(A, 100);

% [Ad,Bd] = c2d(A,B,ts);

Ad = eye(size(A)) + A*ts;
Bd = B*ts;

x = Ad*x(:) + Bd*u(:);

i = x(circ.ii_unipolar) < 0;
x(circ.ii_unipolar(i)) = 0;

y = x;








































