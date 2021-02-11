ccc
odefun = 'LinearPendulumDisc2';

Ts = 0.05;
t = 0:Ts:10;
m = 1;
g = 9.81;
l = 1;
b = 0.2;
parameters = {'mglb', [m g l b]};
fcn_type = 'd';
init_sys = idgrey(odefun, parameters, fcn_type, {}, Ts); 
init_sys.Structure.Parameters.Free(1:3) = false;


% create some data
N = length(t);
b_true = 0.3;
l_true = 0.5;
x0 = [1 1]';
x = x0;
for i = 1:N
  [A,B,C,D] = LinearPendulumDisc2([m g l b_true], Ts);
  x = A*x;
  y(i) = C*x + 0.05*randn;
end
data = iddata(y(:), zeros(N,0), Ts);   

sys = greyest(data,init_sys);

b_est = sys.Structure.Parameters.Value(4);

disp(['b_est: ' num2str(b_est)])
disp(['b_true: ' num2str(b_true)])


scatter(t,y)




















