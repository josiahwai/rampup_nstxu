
% create a hybrid of tok_data_structs

% coneqt_nstxu_obj_config2016_6565.mat is the version from Dan Boyer's coneqt
% repository. nstxu_obj_config2016_6565_original.mat was created using the
% info in its 'make_tok_inputs' field. Copy some of the fields from one to another 
% because GA codes don't always play nicely with each version. 

% The final version is: "nstxu_obj_config2016_6565.mat" and will be the 
% version used by all scripts in this repository. 

saveit = 1;

t1 = load('coneqt_nstxu_obj_config2016_6565.mat').tok_data_struct;
t2 = load('nstxu_obj_config2016_6565_original.mat').tok_data_struct;

fields2copy = {'mcc', 'mcv', 'mpc', 'mpv', 'mvv', 'resc', 'resv'};   % most of these are almost identical, some small numerical differences. resv has a few larger differences. 
tok_data_struct = copyfields(t2, t1, fields2copy, true);

tok_data_struct.mpp_full = t1.mpp;

if saveit
  save('nstxu_obj_config2016_6565.mat', 'tok_data_struct');
end










