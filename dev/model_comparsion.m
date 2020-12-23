ccc
t = 100:10:960;
N = length(t);

md = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/std/';

Aref = load([md '204660_' num2str(t(50)) '_sys.mat']).sys.As;
A{1} = load([md '204660_' num2str(t(1)) '_sys.mat']).sys.As;

for i = 2:N
  t(i)
  A{i} = load([md '204660_' num2str(t(i)) '_sys.mat']).sys.As;
  % d(i) = norm(A{i} - A{i-1}); 
  d(i) = norm(A{i} - Aref);
  f(i) = A{i}(1,1);
end


figure
plot(t,d)

figure
plot(t, f)
















