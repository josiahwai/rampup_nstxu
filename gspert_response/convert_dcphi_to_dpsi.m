clear; clc; close all

saveit = 1;

builds_dir = '/Users/jwai/Desktop/final-jellyfish/josiahwai/pertnet/data/raw_data/pertbuilds_train';

d = dir(builds_dir);
d(1:2) = [];

vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
mpp = tok_data_struct.mpp;

for i = 1:length(d)
  
  disp(i / length(d))
  
  targs = load([d(i).folder '/' d(i).name]).targs;
  
  nsamples = length(targs.actualtime);
  ncoils = size(targs.dcphidix, 3);
  
  targs.dpsidix = 0*targs.dcphidix;
  for isample = 1:nsamples    
    targs.dpsidix(isample,:,:) = mpp * squeeze(targs.dcphidix(isample,:,:));    
  end
  
  if saveit
    save([d(i).folder '/' d(i).name], 'targs')
  end
end
  
%%
% Proof that both mpp methods are the same, even though mpp is slightly different
% (mpp2 method is faster, but tok_data_struct1 was used in gspert)
close all

tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
circ = nstxu2016_circ(tok_data_struct);
mpp1 = tok_data_struct.mpp;

vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
mpp2 = tok_data_struct.mpp;

X = targs.dcphidix;
Y = X*0;

icoil = 6;
isample = 50;

dpsi1 = mpp_x_vec(mpp1, X(isample,:,icoil));
dpsi1 = reshape(dpsi1, 65, 65);

dpsi2 = mpp2 * X(isample,:,icoil)';
dpsi2 = reshape(dpsi2, 65, 65);

figure
hold on
[~,cs] = contour(dpsi1, 20, '--r');
contour(dpsi2, cs.LevelList, 'b')
colorbar



%%
% Proof that inverting dpsidix to obtain dcphidix produces same result as having
% dcphidix. This is useful since  dpsidix has lower dimensional space than dcphidix

X = targs.dcphidix;
Y = targs.dpsidix;

isample = 70;
icoil = 5;

vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
tok_data_struct = vac_sys.build_inputs.tok_data_struct;
mpp = tok_data_struct.mpp;

mppi = inv(mpp);

X2 = mppi * squeeze(targs.dpsidix(isample,:,icoil))';
X2 = reshape(X2, 65, 65);

X1 = reshape(X(isample,:,icoil), 65, 65);

figure
hold on
[~,cs] = contourf(X1, 20);
colorbar

figure
hold on
contourf(X2, 20);
colorbar

figure
hold on
[~,cs] = contour(X1, 20, '--r');
contour(X2, cs.LevelList, 'b')
colorbar



























