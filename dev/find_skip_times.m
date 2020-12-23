ccc
t = 100:10:960;
N = length(t);

md = '/Users/jwai/Research/rampup_nstxu/buildmodel/built_models/std/';

tref = 500;
Aref = load([md '204660_' num2str(tref) '_sys.mat']).sys.As;
A{1} = load([md '204660_' num2str(t(1)) '_sys.mat']).sys.As;

for i = 2:N  
  A{i} = load([md '204660_' num2str(t(i)) '_sys.mat']).sys.As;
  % d(i) = norm(A{i} - A{i-1}); 
  % d(i) = norm(A{i} - Aref);
  f(i,:) = reshape(A{i} - Aref, 1, []);
  Af(i,:,:) =  A{i} - Aref; 
end

close all
fsmooth = smoothdata(f, 'movmedian', 5);
figure
plot(t,fsmooth)

fsmooth = smoothdata(fsmooth, 'movmean', 5);
figure
plot(t,fsmooth)

% figure
% plot(t,d)

% figure
% plot(t, f)
















