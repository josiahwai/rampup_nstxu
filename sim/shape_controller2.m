clear all; clc; close all
warning('off','MATLAB:polyshape:repairedBySimplify')
warning('off', 'curvefit:fit:nonDoubleYData')
warning('off', 'stats:pca:ColRankDefX')

ROOT = getenv('RAMPROOT');

% we will try to match shape of this shot at the shape_time
shot = 204660;  
shape_time = 0.7; 
enforce_stability = 0;

t0 = 0.07;
tf = 0.9;
N = 50;
t = linspace(t0, tf, N);

% load vacuum geometry, inductances, resistances
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
struct_to_ws(tok_data_struct);
circ = nstxu2016_circ(tok_data_struct);
sys = load([ROOT 'sysid/fit_coils_vessels/coil_vessel_fit.mat']).coil_vessel_fit;
res = load([ROOT 'sysid/fit_plasma_resistance/fits_all/res' num2str(shot) '.mat']).res;

% fetch efits - will use efit profiles and shape as targets
[targets, targets_array, efit01_eqs] = read_target(shot, t, tok_data_struct);
[~,i] = min(abs(targets.time - shape_time));
target = targets_array(i);

response = vacuum_response3(target, tok_data_struct);

eq = efit01_eqs.gdata(i-10);
eq_target = efit01_eqs.gdata(i);

y = read_isoflux(eq, target, tok_data_struct);

target.icx = y.icx;
target.psicp_err = 0*target.rcp(:);
target.psi_r = 0;
target.psi_z = 0;

targetvec = [target.icx; 
             target.psicp_err; 
             target.psi_r; 
             target.psi_z];

yvec = [y.icx; 
        y.psicp - y.target_bdef_psi; 
        y.target_bdef_psi_r; 
        y.target_bdef_psi_z];

C = [response.disdis(circ.iicx,:); 
     response.dpsicpdis - response.dpsibrydis; 
     response.dpsibrydis_r; response.dpsibrydis_z];

C = C(:,circ.iicx);
   
wt.icx = ones(circ.ncx,1) * 1e-6;
wt.icx(circ.iicx_remove) = 1e8;
wt.psicp = ones(size(target.rcp(:))) * 1;
wt.psi_r = 1;
wt.psi_z = 1;

W = diag([wt.icx; wt.psicp; wt.psi_r; wt.psi_z]);
   
dy_target = targetvec - yvec;
di_target = pinv(W*C)*W*dy_target;

%%
mpcx = mpc * circ.Pcc;

dpsizr_app = reshape(mpcx*di_target, nz, nr);

psizr_new = eq.psizr + dpsizr_app;

[rx,zx,psix] = isoflux_xpFinder(psizr_new,0.6,-1,rg,zg);

figure
plot_eq(eq_target)
contour(rg,zg,psizr_new,[psix psix], 'b')
contour(rg,zg,eq.psizr,[eq.psibry eq.psibry], 'g')

%%
sys = load([ROOT 'sysid/fit_coils_vessels/coil_vessel_fit.mat']).coil_vessel_fit;

Mxxi = inv(sys.Mxx);
A = -Mxxi*diag(sys.Rxx);
B = Mxxi(:,circ.ikeep);

A = A(circ.ikeep, circ.ikeep);
B = B(circ.ikeep, :);

wt.vc = ones(circ.ncx_keep,1) * 500;
wt.ic_tracking = ones(circ.ncx_keep,1) * 1;

Q = diag(wt.ic_tracking);
R = diag(wt.vc);

Klqr = lqr(A,B,Q,R);


I_target = y.icx + di_target;
dv = Klqr * di_target(circ.ikeep);
v0 = sys.rc(circ.ikeep) .* I_target(circ.ikeep);

v = v0 + dv;
















































