ccc
odefun = 'LinearPendulumDisc';

Ts = 0.05;
t = 0:Ts:10;
m = 1;
g = 9.81;
l = 1;
b = 0.2;
parameters = {'mass',m;'gravity',g;'length',l;'friction',b; 'time', t};
fcn_type = 'd';
init_sys = idgrey(odefun, parameters, fcn_type, {}, Ts); 
init_sys.Structure.Parameters(1).Free = false;
init_sys.Structure.Parameters(2).Free = false;
init_sys.Structure.Parameters(3).Minimum = 0;  % estimate l 
init_sys.Structure.Parameters(4).Minimum = 0;  % estimate b
init_sys.Structure.Parameters(5).Free = false;

% create some data
N = length(t);
b_true = 0.3;
l_true = 0.5;
x0 = [1 1]';
x = x0;
for i = 1:N
  [A,B,C,D] = LinearPendulumDisc(m,g,l_true,b_true,t(i),Ts);
  x = A*x;
  y(i) = C*x + 0.05*randn;
end
data = iddata(y(:), zeros(N,0), Ts);   

sys = greyest(data,init_sys);

l_est = sys.Structure.Parameters(3).Value;
b_est = sys.Structure.Parameters(4).Value;

disp(['l_est: ' num2str(l_est)])
disp(['l_true: ' num2str(l_true)])
disp(['b_est: ' num2str(b_est)])
disp(['b_true: ' num2str(b_true)])


scatter(t,y)




















