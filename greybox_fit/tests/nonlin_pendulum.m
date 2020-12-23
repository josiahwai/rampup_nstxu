
function [A,B,C,D] = nonlin_pendulum(m,g,l,bvec,t,sim_timebase,Ts)

b = bvec(sim_timebase==t);

Acont = [0 1; -g/l, -b/m/l^2];
Bcont = zeros(2,0);
C = [1 0];
D = zeros(1,0);

[A,B] = c2d(Acont,Bcont,Ts);

end 