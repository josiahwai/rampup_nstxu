coils1 = load('coils204660_1.mat').coils204660;
coils2 = load('coils204660_2.mat').coils204660;

t = [coils1.t coils2.t(2:end)];
fit_options = fitoptions('Method','Smooth','SmoothingParam',0.99999);
  
fn = fieldnames(coils1);

for i = 1:length(fn)  
  y = [coils1.(fn{i}) coils2.(fn{i})(:, 2:end)];  
  
  yfit = 0*y;
  for j = 1:size(y,1)    
    f = fit(t(:), y(j,:)', 'smoothingspline', fit_options);
    yfit(j,:) = f(t);
  end
  

  coils.(fn{i}) = yfit;
  
  figure 
  hold on
  plot(coils.t, y', '--r')
  plot(coils.t, yfit', 'b')
  title(fn{i})
end



save('coils204660', 'coils')


