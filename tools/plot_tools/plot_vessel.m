
load('nstxu_obj_config2016_6565.mat')
circ = nstxu2016_circ(tok_data_struct);

plot_nstxu_geo(tok_data_struct);


for i = 1:40
  
  j = find(circ.vvgroup==i);
  r = tok_data_struct.vvdata(2,j);
  z = tok_data_struct.vvdata(1,j);
  
  scatter(r,z,100,'filled')
  text(r(1), z(1), [num2str(i) ' ' circ.vvnames{i}], 'fontsize', 18)
end









