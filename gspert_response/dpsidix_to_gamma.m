clear all; clc; close all
% Reconstruct gamma from psi-response

shot = 204660;
isample = 44;   % 25 90

% load stuff
build_dir = '/Users/jwai/Research/rampup_nstxu/gspert_response/pertbuilds_train/';
targs = load([build_dir 'train_response_' num2str(shot) '.mat']).targs;
vac_sys = load('NSTXU_vacuum_system_fit.mat').NSTXU_vacuum_system_fit;
mpp = vac_sys.build_inputs.tok_data_struct.mpp;
tok_data_struct = load('nstxu_obj_config2016_6565.mat').tok_data_struct;
mpc = tok_data_struct.mpc;
mpv = tok_data_struct.mpv;
circ = nstxu2016_circ(tok_data_struct);
mppi = inv(mpp);

% tree = 'EFIT01';
% tokamak = 'nstxu';
% server = 'skylark.pppl.gov:8501';
% eq = read_eq(shot, targs.actualtime(isample), tree, tokamak, server);

load('eq204660_250.mat')


eq.ecturn = tok_data_struct.ecnturn;
eq.ecid   = ones(size(eq.ecturn));
eq.turnfc = tok_data_struct.fcnturn';
eq.fcturn = circ.fcfrac;
eq.fcid = circ.fccirc';      

% response
build_inputs.tokamak = 'NSTXU';
build_inputs.vacuum_objs = tok_data_struct;
build_inputs.ichooseq = 4;
build_inputs.irzresp_dynamic = 5;
build_inputs.irzresp_output = 5;
build_inputs.iplcirc = 1;
build_inputs.cccirc = circ.cccirc(:);
build_inputs.vvcirc = circ.vvcirc(:);
build_inputs.vvgroup = circ.vvgroup(:);  
build_inputs.equil_data = eq;

sys = build_tokamak_system(build_inputs);
delete('NSTXU_netlist.dat')
      

Rxx = sys.Pxx'*sys.rxx*sys.Pxx;
rxx = diag(Rxx);
rxx(circ.iremove) = 10;
Rxx = diag(rxx);

amat = -inv(sys.lstar) * Rxx;

e = esort(eig(amat));
e(1:5)


%%
% convert dpsi to dcphi
dcphidix = squeeze(targs.dcphidix(isample,:,:));
dpsidix = squeeze(targs.dpsidix(isample,:,:));
dcphidix2 = mppi * dpsidix;


% get xmat from dcphi
P = circ.Pxx(1:end-1, 1:end-1);
dcphidis = sys.gspert_data.dcphidis;

xmat = P' * [mpc'; mpv'] * dcphidis * P;
mmat = P' * sys.mxx(1:end-1,1:end-1) * P;
r = P' * sys.rxx(1:end-1,1:end-1) * P;
r = diag(r);
r(circ.iremove) = 10;
r = diag(r);
amat = -inv(xmat + mmat) * r;

e = esort(eig(amat));
e(1:5)

targs.gamma(isample)

%%
isample = 90
dcphidicx = squeeze(targs.dcphidix(isample,:,:)) * circ.Pcc_keep;
dcphidiv = sys.gspert_data.dcphidis(:,circ.iiv);
dcphidivx = dcphidiv * circ.Pvv;

dcphidix = [dcphidicx dcphidivx];

xmat = P' * [mpc'; mpv'] * dcphidix;
mmat = P' * sys.mxx(1:end-1,1:end-1) * P;
r = P' * sys.rxx(1:end-1,1:end-1) * P;
amat = -inv(xmat + mmat) * r;

e = real(esort(eig(amat)));
e(1:5)

%%
dcphidix =inv(Mpp) * dpsidix;
X = Pxx' * MpcMpv' * dcphidix;
amat = -inv(M + X) * Rxx;
e = esort(eig(amat));
e(1:5)







































