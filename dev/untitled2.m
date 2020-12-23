clear; clc; warning('off', 'all'); close all
root = '/u/jwai/iter_advanced_divertor_2019/current/';
addpath(genpath('/u/jwai/toksys/current/matlab'));
addpath(genpath(root));
noise_dir =[root 'inputs/noise/'];

%=========
% SETTINGS
%=========

use_previous_noise = 0;


%==============================================
% Load some stuff & build the state space model
%==============================================
load('eq15MA')
load('iter_obj_2010v3')
tok_data_struct.resc(1:12) = tok_data_struct.resc(1:12)/1e6; % modify for s/c coils

% Dimensions of coil and vessel elements, outputs, grid
% .....................................................

npfc = 12;  % 12 Poloidal Field coils
npfx = 11;  % 11 PF circuits (9 and 10 connected in series)
nvsc = 2;   %  2 Vertical Stabilization coils
nvsx = 1;   %  1 VS circuit (vs1 and vs2 connected in antiseries)
nvvc = 148; % 148 vacuum vessel coils 
nvvx = 148; % 148 vacuum vessel circuits
nip = 1;    %  1 plasma circuit
nc = npfc + nvsc + nvvc + nip;  % total # coils
nx = npfx + nvsx + nvvx + nip;  % total # circuits

npfu = npfx; % # control voltages for PF coils
nvsu = nvsx; % # control voltages for VS coils
nu   = npfu + nvsu; % total # control inputs

nxr = npfx + nip;  % reduced set of circuits for mpc prediction model
ixr = [1:npfx nx]; % index of reduced set of circuits

Pcx = blkdiag(eye(8), [1 1]', eye(2), [1 -1]'); % mapping from circuit currents to coil currents
Pxc = blkdiag(eye(8), [1 0], eye(2), [1 0]);    % mapping from coil currents to circuit currents

% state vector (of circuits) is organized as [Ipf Ivs Ivv Ip]. Index for
% each subset of circuits
ipfx = 1:npfx;
ivsx = npfx + 1: npfx + nvsx;
ivvx = npfx + nvsx + 1: npfx + nvsx + nvvx;
iipx = nx;

% Outputs are:
% Ip, 3*rx, 3*zx, 2*spR, 2*spZ, psibry, fexp, psi @ xp2, dpsidr @ xp2
% psi @ nSegs, current at nxr circuits
nSegs = 31;
ne = nSegs + 15; 
ny = ne + nxr;

% Grid size
rg = eq.rg';
zg = eq.zg;
nr = length(rg);
nz = length(zg);

control_dims = struct('npfc',npfc,'npfx',npfx,'nvsx',nvsx,'nvsc',nvsc,...
  'nvvc',nvvc,'nvvx',nvvx,'nip',nip,'nc',nc,'nx',nx,'npfu',npfu,...
  'nvsu',nvsu,'nu',nu,'nxr',nxr,'ixr',ixr,'Pcx',Pcx,'Pxc',Pxc,'ipfx',ipfx,...
  'ivsx',ivsx,'ivvx',ivvx,'iipx',iipx,'nSegs',nSegs,'ne',ne,'ny',ny,...
  'rg',rg,'zg',zg,'nr',nr,'nz',nz);


% Initialize Model
% ................

%Define Control Segments
segData = find_isosegs(tok_data_struct, nSegs);
cpRZD = isoflux_cpFinder(eq.psizr, eq.psibry, rg, zg, segData);
mag0 = cp_rz2mag(cpRZD, segData);
magLim = max(.08, mag0 - .4);
cpRZLim = cp_mag2rz(magLim,segData);

% initialize the grad-shafranov solver
IlimAll = [48 55 55 55 52 52  45  45  45  45  45  45 10 10]'*1000;
Ilim = IlimAll;
Ilim([9,13,14]) = [];
[spec, init, config] = initialize_gsdesign(tok_data_struct, eq, 1.4*IlimAll, segData);

%Integration Timing
tstart = 0.00;
tend   = 150;
dt = 1e-3;
times = linspace(tstart, tend, (tend-tstart)/dt);
nTimes = length(times);

%Timing for control updates & discrete time model
controlStep = 1000;
ts = controlStep*dt;
nControlSteps = floor(nTimes/controlStep);
controlTimes = linspace(tstart, tend, nControlSteps);
cc = 0;

%Timing for plasma model updates
plasmaStep = 1000;
nPlasmaSteps = floor(nTimes/plasmaStep);
plasmaTimes = linspace(tstart, tend, nPlasmaSteps);
pp = 0;


% Measurement noise standard deviations
ipNoise                = 50000;  % A
pfcoilMsmtNoise(1:npfx)  = 250;    % A
vscoilMsmtNoise      = 50; % A
xpPosNoise(1:6)      = 0.02;   % m
spRNoise(1:2)        = 0.02;   % m
spZNoise(1:2)        = 0.02;   % m
psiBryNoise          = 1;      % wb
fexpNoise            = 0.4;    
gradpsix2Noise       = 4;      % wb/m
dpsibrypsix2Noise    = 1;      % wb
cpPsiNoise(1:nSegs)  = 1;      % wb

stateMsmtNoise  = [pfcoilMsmtNoise vscoilMsmtNoise ipNoise]';
outputMsmtNoise = [ipNoise xpPosNoise spRNoise spZNoise psiBryNoise ...
  fexpNoise gradpsix2Noise dpsibrypsix2Noise cpPsiNoise]';
stdV = [stateMsmtNoise; outputMsmtNoise];
Rv = diag(stdV.*stdV); % covariance


% Process noise standard deviations
coilProcNoise(1:npfx)  = 100;
ipProcNoise            = 2000;
vsProcNoise            = 10;
vesselProcNoise(1:nvvx) = 10;

stdW  = 0.3*[coilProcNoise vsProcNoise vesselProcNoise  ipProcNoise]'; % std dev
Qw = diag(stdW.*stdW);                                     % process covariance
Pk = 10 * Qw;   % error covariance matrix (for Kalman filter)

if use_previous_noise
  load([noise_dir 'allMsmtNoise.mat']);
  load([noise_dir 'allProcessNoise.mat']);
else  
  allMsmtNoise = randn(nControlSteps,length(stdV)) * diag(stdV) * 0;
  allProcessNoise = randn(nControlSteps,length(stdW)) * diag(stdW) * 0;
end

% Initialize plant and estimator models
I_init = [Pxc * eq.ic; zeros(nvvx, 1); eq.cpasma];  
if use_previous_noise
  load([noise_dir 'initialCoilNoise.mat']);
else
  initialCoilNoise = randn(nx,1) .* [pfcoilMsmtNoise'; vscoilMsmtNoise; ones(nvvx,1); ipNoise] * 0;
end
estI_init = I_init + initialCoilNoise;
plt = initializeModel(I_init, nx, nxr, eq, Pcx, Pxc);
est = initializeModel(estI_init,nx, nxr, eq, Pcx, Pxc);
estsim = [];
pltsim = [];

% ===============
% MPC formulation
% ===============

% Move blocking matrix with geometrically increasing step sizes
tLookAhead = 9.2;
nb = 10;
N = tLookAhead/dt/controlStep;

% find geometric scaling factor b
syms b i
eqn = symsum(b^i,i,1,nb) == N;
bvals = vpasolve(eqn, b);
clear b i
b = bvals(abs(imag(bvals)) < 2*eps & bvals > 0); % real positive soln

% create move block matrix
blks = [];
for i = 1:nb, blks(i) = round(b^i); end
N = sum(blks);
addblks = [0 cumsum(blks)];
Tb = zeros(N, nb);
for j = 1:nb
  i = addblks(j)+1:addblks(j+1);
  Tb(i,j) = 1;
end
mb = kron(Tb, eye(npfu));


% cost matrices
% pos wts: [Ip, rx1, rx2, rx3, zx1, zx2, zx3, rsp_in, rsp_out, zsp_in, zsp_out]
% psi wts: [psibry psibry-psix3 cpPsi]
load('wt');

R = 1e-4 * diag([1 1 1 1 1 1 1 1 1 1 1]);
Rhat = kron(eye(N),R);
Q = diag(wt.^2);

% Initialize
opt    = mpcqpsolverOptions;
simErrorFlag = false;
simEndFlag   = true;
Iviol  = [];
pp = 0;
cc = 0;
vprev = zeros(npfu,1);

tt0 = 0;

%====================
% Integrate the Model
%====================
for tt = tt0:nTimes-1
  % try
  %========================
  % Update the Plasma Model
  %========================
  if mod(tt,plasmaStep) == 0
    
    pp = pp + 0;
    
    % Relinearize plant/estimator plasma models
    plt = updateModel_fexp1(plt, tok_data_struct, control_dims, segData,...
      tt, dt, spec, init, config, cpRZD, ts);
    
    est = plt;
%     est = updateModel_fexp1(est, tok_data_struct, ncx, nvx, nxr, nu, ny, nSegs, segData,...
%       tt, dt, spec, init, config, cpRZD, Pcx);
    
    % ===============
    % MPC Formulation
    % ===============
    % discrete augmented matrices
    A = est.amat_sdr;
    B = est.bmat_pf_sdr;
    C = est.cmat_r;
    
    ChatAug = kron(eye(N), est.cmat_aug_r);
    Apow  = eye(nxr);
    
    % prediction model
    F  = zeros(N*size(B));
    OB  = [];
    for i = 1:N
      F  = F  + kron(diag(ones(N-i+1,1), -i+1), Apow*B);
      Apow  = A*Apow;
      OB  = [OB; Apow];
    end
    
    
    % updated cost & targ
    Atemp = [A zeros(nxr,ne); C*A zeros(ne)];
    Btemp = [B; C*B];
    Qbar = dare(Atemp,Btemp,Q,R);
    Qhat = blkdiag(kron(eye(N-1),Q), Qbar);
    
    
    % quadprog
    H = mb' * (F'*ChatAug'*Qhat*ChatAug*F + Rhat) * mb;
    if max(max(H-H')) < 1e-10, H = (H+H')/2; end
    Linv = inv(chol(H,'lower'));
    Linv(abs(Linv) < 1e-10) = 0; % set small values to zero
    
  end
  
  
  %=====================
  % Update control input
  %=====================
  if mod(tt, controlStep) == 0
    
    cc = cc + 1;
    
    fprintf('Iteration %d of %d...\n', tt/controlStep+1, max(nControlSteps,1))
    
    
    % INPUT RATE CONSTRAINTS
    uCirc = [1 1 1 1 1 1 1 1 2 1 1]';   % 2x power for CS1L/CS1U
    dvMax =  200*uCirc;
    dvMin = -200*uCirc;
    G0 = [eye(npfx); -eye(npfx)];
    G0hat = kron(eye(N), G0);
    g0 = [dvMax; -dvMin];
    g0hat = repmat(g0,N,1);
    alpha = diag(ones(N,1)) + diag(-1*ones(N-1,1), -1);
    Sv = kron(alpha, eye(npfx));
    v0hat = repmat(est.v0(ipfx),N,1);
    Aineq0 = G0hat*Sv*mb;
    bineq0 = g0hat - G0hat*Sv*v0hat + G0hat*[vprev; zeros((N-1)*npfu,1)];
    
    
    % INPUT CONSTRAINTS
    vmax =  1500*uCirc;
    vmin = -1500*uCirc;
    G1    = [eye(npfu); -eye(npfu)];
    G1hat = kron(eye(N), G1);
    Aineq1 = G1hat*mb;
    g1 = [vmax; -vmin];
    g1hat = repmat(g1,N,1);
    bineq1 = g1hat - G1hat*v0hat;
    
    
    % STATE AND OUTPUT CONSTRAINTS
    Imax = [ Ilim;  15.5e6];
    Imin = [-Ilim;  14.5e6];
    dymax =  [Imax;  inf(ne,1)];
    dymin =  [Imin; -inf(ne,1)];
    
    % min wall gapsf
    psizr = est.eq.psizr;
    for ii = 1:nSegs
      cpPsiLim(ii) = bicubicHermite(rg, zg, psizr, cpRZLim(ii,1), cpRZLim(ii,2));
    end
    dymax(end-nSegs+1:end) = est.cpPsiD - cpPsiLim';
    
    % outer strike point on plate
    dymin(nxr+11) = -4.55 - est.zD(11);
    dymin(nxr+9)  =  5.55 - est.zD(9);
    
    
    % form the constraint
    y0 = [est.I0(ixr); est.e0];
    y0hat = repmat(y0,N,1);
    g2 = [dymax; -dymin];
    g2hat = repmat(g2,N,1);
    G2 = [eye(ny); -eye(ny)];
    G2hat = kron(eye(N), G2);
    Aineq2 = G2hat*ChatAug*F*mb;
    bineq2 = g2hat - G2hat*y0hat - G2hat*ChatAug*OB*est.xk;
    
    
    % POWER DERIVATIVE CONSTRAINT (takes same form as dv constraint G4 * dv <= g4)
    g4 = [200; 200]*1e6;  % MW/s
    Ifactor = Ilim' .* sign(est.I(ipfx)');
    Ifactor = est.I(ipfx)' * 1.04;
    G4 = Ifactor / (dt * controlStep);
    G4 = [G4; -G4];
    G4hat = kron(eye(N),G4);
    g4hat = repmat(g4,N,1);
    Aineq4 = G4hat*Sv*mb;
    bineq4 = g4hat - G4hat * (Sv*v0hat + [vprev; zeros((N-1)*npfu,1)]);
    
    
    % IGNORE constraint sets (debugging)
    %     bineq0 = inf(size(bineq0));  % dv
    %     bineq1 = inf(size(bineq1));  % v
    %     bineq4 = inf(size(bineq4));  % pdot
    %     bineq2 = inf(size(bineq2));  % y
    
    
    Aineq = [Aineq0; Aineq1; Aineq2; Aineq4];
    bineq = [bineq0; bineq1; bineq2; bineq4];
    
    
    % future targets
    rhat = [];
    for k = 1:N     
      new_t = (tt+k*controlStep)*dt;      
      
      zD_old = targets_fexp1(tt*dt, est.eq.psibry);
      zD_new = targets_fexp1(new_t, est.eq.psibry);           
      
      r = -est.e0 + zD_new - zD_old;
      
      % encourage a steady state solution
      tss = 120; 
      if tt*dt > tss
        if ~exist('Iss','var'), Iss = [est.I0(ipfx); 0]; end
        xtarg = Iss - [est.I0(ipfx); 0];
      else
        xtarg = zeros(nxr,1);
      end           
      
      rhat = [rhat; xtarg; r];
    end
    
    
    % solve Quadratic Program
    f = mb'*F'*ChatAug'*Qhat*(ChatAug*OB*est.xk-rhat);
    if ~exist('iA','var')
      iA = false(size(bineq));
      iA_all = iA;
    end
    [Uhat,status,iA] = mpcqpsolver(Linv,f,-Aineq,-bineq,[],zeros(0,1),iA,opt);
    
    iA_all = iA_all | iA;
    
    % Error Handling
    if status > 0   % valid solution - save results
      lastValid.cc = cc;
      lastValid.V = mb*Uhat + v0hat;
      u = Uhat(ipfx);
    else    % invalid solution - use last valid results
      warning('on','all')
      warning('No QP solution.');
      warning('off', 'all')
      Uhat = lastValid.V - v0hat;
      dcc = cc - lastValid.cc;          % how many steps since last valid
      i = dcc*npfx + 1: (dcc+1) * npfx;
      u = Uhat(i);
    end
    
    est.u = u;
    plt.u = est.u + est.v0 - plt.v0;
    v = est.u + est.v0;
    plt.v = v;
    est.v = v;
    
    % =============
    % Kalman filter 
    % =============
%     [est.amatD, est.bmatD] = c2d(est.amat,est.bmat,ts);
%     [plt.amatD, plt.bmatD] = c2d(plt.amat,plt.bmat,ts);
%     
%     % measure
%     cmat_states = [eye(nu) zeros(nu,nx-nu); zeros(1,nx-1) 1];
%     plt.cmat_aug = [cmat_states; plt.cmat];
%     est.cmat_aug = [cmat_states; est.cmat];
%     
%     vk     = allMsmtNoise(cc,:)';
%     plt.dy = plt.cmat_aug * plt.dI + vk;
%     est.dy = est.cmat_aug * est.dI;
%     
%     plt.y = plt.dy + [plt.I0(ixr); plt.e0];
%     est.y = est.dy + [est.I0(ixr); est.e0];
%     
%     % predict
%     Pk_kprev = est.amatD*Pk*est.amatD' + Qw;
%     
%     % optimal filter gain
%     L = Pk_kprev * est.cmat_aug' * inv(est.cmat_aug * Pk_kprev * est.cmat_aug' + Rv);
%     
%     % modify residual
%     residual = plt.y - est.y;
%     residual([15 16 18 19]) = 0; % dont use xp2/xp3 (jumps between diff xps)
%     residual([21 23]) = 0; % dont use outer sp (nonlinear corner)
%     
%     % update estimator
%     est.dI = est.amatD*est.dI + est.bmatD*est.u + L*residual;
%     est.I  = est.I0 + est.dI;
%     est.xk = est.dI(ixr);
%     
%     % update innovation
%     Pk = (eye(nx) - L*est.cmat_aug)*Pk_kprev;
%     
%     % update plant
%     wk = allProcessNoise(cc,:)';
%     plt.dI = plt.amatD*plt.dI + plt.bmatD*plt.u + wk;
%     plt.I = plt.I0 + plt.dI;
%     plt.xk = plt.dI(ixr);
    
    est.dI = est.amat_sd*est.dI + est.bmat_pf_sd*est.u;
    est.I  = est.I0 + est.dI;
    est.xk = est.dI(ixr);
    plt = est;
    
    % ==============
    % Save iteration
    % ==============
    
    % write stuff to sim structs
    Iviol = unique([Iviol; find(abs(plt.I(ipfx)) > Imax(ipfx))]);
    
    power(cc) = plt.I(ipfx)' * v / 1e6;       % <250 MW
    powerDot(cc) = plt.I(ipfx)' * (v-vprev) / (controlStep*dt) / 1e6;  % < 200 MW/s
    vprev = v;
    
    pltsim = writeToSims(pltsim,plt,ne,tt,controlStep,plasmaStep,tok_data_struct,...
      segData, controlTimes, plasmaTimes, Iviol, simEndFlag);
    
    estsim = writeToSims(estsim,est,ne,tt,controlStep,plasmaStep,tok_data_struct,...
      segData, controlTimes, plasmaTimes, Iviol, simEndFlag);
    
  end
  
  % catch err
  %      fprintf('Error occurred on line %d. Message was: \n\n %s', err.stack(:).line, err.message);
  %      break
  % end
end


% Make some plots
simEndFlag = true;
pltsim.power = power;
pltsim.powerDot = powerDot;
pltsim.iA = iA_all;
pltsim.allMsmtNoise = allMsmtNoise;
pltsim.allProcessNoise = allProcessNoise;

pltsim = writeToSims(pltsim,plt,ne,tt,controlStep,plasmaStep,tok_data_struct,...
  segData, controlTimes, plasmaTimes, Iviol, simEndFlag);

estsim = writeToSims(estsim,est,ne,tt,controlStep,plasmaStep,tok_data_struct,...
  segData, controlTimes, plasmaTimes, Iviol, simEndFlag);


if cc > 20, save('output/pltsim.mat', 'pltsim'), end
if cc > 20, save('output/estsim.mat', 'estsim'), end
delete('ITER_netlist.dat')











