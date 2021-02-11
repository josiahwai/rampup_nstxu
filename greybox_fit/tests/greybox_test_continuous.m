
ccc
odefun = 'LinearPendulumCont';

m = 1;
g = 9.81;
l = 1;
b = 0.2;
parameters = {'mass',m;'gravity',g;'length',l;'friction',b};
fcn_type = 'c';
init_sys = idgrey(odefun,parameters,fcn_type); 
init_sys.Structure.Parameters(1).Free = false;
init_sys.Structure.Parameters(2).Free = false;
init_sys.Structure.Parameters(3).Minimum = 0;  % estimate l 
init_sys.Structure.Parameters(4).Minimum = 0;  % estimate l 

% create some data
dt = 1e-3;
t = 0:dt:10;
N = length(t);
b_true = 0.3;
l_true = 0.5;
x0 = [1 1]';
x = x0;
for i = 1:N
  [A,B,C,D] = LinearPendulum2(m,g,l_true,b_true);
  xdot = A*x;
  x = x + xdot*dt;
  y(i) = C*x + randn*0.05;
end
data = iddata(y(:), zeros(N,0), dt);   

sys = greyest(data,init_sys);

l_est = sys.Structure.Parameters(3).Value;
b_est = sys.Structure.Parameters(4).Value;

disp(['l_est: ' num2str(l_est)])
disp(['l_true: ' num2str(l_true)])
disp(['b_est: ' num2str(b_est)])
disp(['b_true: ' num2str(b_true)])



function [A,B,C,D] = LinearPendulum2(m,g,l,b,Ts)
A = [0 1; -g/l, -b/m/l^2];
B = zeros(2,0);
C = [1 0];
D = zeros(1,0);
end 
















