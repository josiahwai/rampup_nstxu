
function [A,B,C,D] = LinearPendulumDisc2(mglb, Ts)

m = mglb(1);
g = mglb(2);
l = mglb(3);
b = mglb(4);

Acont = [0 1; -g/l, -b/m/l^2];
Bcont = zeros(2,0);
C = [1 0];
D = zeros(1,0);

[A,B] = c2d(Acont,Bcont,Ts);
end 