coils1 = load('coils204660.mat').coils;
t1 = coils1.t;


shot = 204660;
coils2 = mds_fetch_current_voltage(shot, 0);
t2 = coils2.t;

load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);

ic1 = coils1.ic(:, circ.ikeep);
ic2 = coils2.ic;

close all
figure
for i = 1:8
  subplot(2,4,i)
  hold on
  plot(t1,ic1(:,i),'r')
  plot(t2,ic2(i,:),'--b')
  title(circ.keep_coils{i})
  xlim([0 1])
end












