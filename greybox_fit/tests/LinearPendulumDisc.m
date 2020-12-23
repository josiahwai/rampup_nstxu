
function [A,B,C,D] = LinearPendulumDisc(m,g,l,b,t,Ts)
Acont = [0 1; -g/l, -b/m/l^2];
Bcont = zeros(2,0);
C = [1 0];
D = zeros(1,0);

[A,B] = c2d(Acont,Bcont,Ts);
end 