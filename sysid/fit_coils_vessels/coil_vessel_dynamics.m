% xdot = A*x + B*u
% x = [ic iv]
% u = [vc; 0] - [mcIp; mvIp]*ipdot 
%
% Definitions: ic=coil currents, iv=vessel currents, vc=power supply
% voltages, [mcIp; mvIp] = mutual inductances of coils and vessels with
% plasma circuit, ipdot=d/dt(plasma current)


function [x,y] = coil_vessel_dynamics(t, x, u, rc, rv, lc, lv, mcc_triu, mvv_triu, mcv, file_args)

args = file_args{1};
circ = args.circ;
ts = args.ts;

% reconstruct mutual inductances from parameters
mcc = mcc_triu + mcc_triu' + diag(lc);
mvv = mvv_triu + mvv_triu' + diag(lv);
M = [mcc mcv; mcv' mvv];
R = diag([rc; rv]);

Minv = inv(M);

A = -Minv*R;
B = [Minv(:,circ.iicx) Minv];

% reduce the speed of the super-fast stable poles, which causes numerical
% (not physical) instabilities
A = numerically_stabilize(A, 100);  

% [Ad,Bd] = c2d(A,B,ts);

Ad = eye(size(A)) + A*ts;
Bd = B*ts;

x = Ad*x(:) + Bd*u(:);

i = x(circ.ii_unipolar) < 0;
x(circ.ii_unipolar(i)) = 0;

y = x;








































