ccc
load('matlab.mat')

uhat = gradient(xhat(circ.iicx,:)) / Ts;

E = [];
npv = circ.nvx + circ.np;
Apow = eye(npv);
F = [];
F_row = zeros(npv, N*circ.ncx);
ivp = ivp0;

for i = 1:N
  A = Ad_list{i};
  B = Bd_list{i};    
  
  idx = (circ.ncx*(i-1)+1):(circ.ncx*i);
  F_row = A * F_row;
  F_row(:,idx) = B;
  F = [F; F_row];
  
  Apow = A*Apow;
  E = [E; Apow];
  
  ivp = A*ivp + B*uhat(:,i);
  ivp_hat2(:,i) = ivp;
end

uhat = reshape(uhat, [], 1);
ivp_hat = E*ivp0 + F*uhat;
ivp_hat = reshape(ivp_hat, [], N);

figure
plot(t, ivp_hat)

figure
plot(t, ivp_hat2)

%%
nx = circ.nvx + circ.np;
nu = circ.ncx;


E = [];
Apow = eye(nx);
F_row = zeros(nx, N*nu);
F = zeros(N*nx, N*nu);

for i = 1:N
  A = Ad_list{1};
  B = Bd_list{1}; 
  
  idx = (nu*(i-1)+1):(nu*i);
  F_row = A * F_row;
  F_row(:,idx) = B;
  F = [F; F_row];
  
  Apow = A*Apow;
  E = [E; Apow];
end
  

%%
E = [];
Apow = eye(2);
F = zeros(N*2, N*1);
x0 = [1; 1];
x = x0;
uhat = 0.1*ones(1,N);

for i = 1:N
  A = [0.95 -0.7; 0.7 0.5];
  B = [1; 0];
  
  m = diag(ones(N-i+1,1), 1-i);  
  F = F + kron(m, Apow*B); 
  Apow = A*Apow;
  E = [E; Apow];
  
  u = uhat(:,i);
  x = A*x + B*u;
  xsim(:,i) = x;
end


uhat = reshape(uhat, [], 1);
xsim2 = E*x0 + F*uhat;
xsim2 = reshape(xsim2, [], N);


figure
plot(t, xsim)

figure
plot(t, xsim2)


%%
nx = 2;
nu = 1;


E = [];
Apow = eye(nx);
F_row = zeros(nx, N*nu);
F = [];

for i = 1:N
  A = Ad_list{1};
  B = Bd_list{1}; 
  
  idx = (nu*(i-1)+1):(nu*i);
  F_row = A * F_row;
  F_row(:,idx) = B;
  
  F = [F; F_row];
  
  
  Apow = A*Apow;
  E = [E; Apow];
end






























































