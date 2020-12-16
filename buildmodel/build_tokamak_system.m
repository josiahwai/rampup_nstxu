function tok_system = build_tokamak_system(build_inputs)
%
% BUILD_TOKAMAK_SYSTEM 
%
%   Function to build linear dynamics and output equations for 
%   electromagnetic model of axisymmetric tokamak system.
%
% USAGE: build_tokamak_system.m
%
% INPUTS:
%
%   build_inputs...structure containing the following fields
%
%       - tok_data_struct...TokSys struct describing the tokamak geometry
%
%       - cccirc............connection vector for the coils
%
%       - ccxnames..........character array containing coil circuit names
%                           (optional)
%
%       - vvgroup...........grouping vector combining individual vessel
%                           elements into grouped conductors
%
%       - vvcirc............connection vector for the vessel groups
%
%       - vvxnames..........character array containing vessel circuit names
%                           (optional)
%
%       - plasma_resp.......flag for selecting the plasma response model 
%                               [0] no plasma              (vacuum) 
%                               [1] linear rigid r -motion (rzrigid)
%                               [2] linear rigid z -motion (rzrigid)
%                               [3] linear rigid rz-motion (rzrigid)
%                               [4] linear nonrigid        (gspert)
%
%       - eq................struct containing equilibrium data 
%
% OUTPUTS: 
%
%   tok_system...structure containing dynamics and output equation objects
%
% AUTHOR: Patrick J. Vail
%
% DATE: 06/20/2017
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 06/20/2017
%
%..........................................................................

%--------------------------------------------------------------------------
% Extract data from input struct build_inputs
%--------------------------------------------------------------------------

if ~isfield(build_inputs, 'tok_data_struct')
    error('Input build_inputs must contain field ''tok_data_struct''')
end
if ~isfield(build_inputs, 'cccirc')
    error('Input build_inputs must contain field ''cccirc''')
end
if ~isfield(build_inputs, 'vvgroup')
    error('Input build_inputs must contain field ''vvgroup''')
end
if ~isfield(build_inputs, 'vvcirc')
    error('Input build_inputs must contain field ''vvcirc''')
end
if ~isfield(build_inputs, 'plasma_resp')
    error('Input build_inputs must contain field ''plasma_resp''')
end

tok_data_struct = build_inputs.tok_data_struct;
cccirc          = build_inputs.cccirc;
vvgroup         = build_inputs.vvgroup;
vvcirc          = build_inputs.vvcirc;
plasma_resp     = build_inputs.plasma_resp;

if isfield(build_inputs, 'ccxnames')
    ccxnames = build_inputs.ccxnames;
end
if isfield(build_inputs, 'vvxnames')
    vvxnames = build_inputs.vvxnames;
end

if plasma_resp ~= 0
    if ~isfield(build_inputs, 'eq')
        error('Input build_inputs must contain field ''eq''')
    else
        eq = build_inputs.eq;
    end
end

nc = tok_data_struct.nc;
nv = tok_data_struct.nv;

nz = tok_data_struct.nz;
nr = tok_data_struct.nr;

nfl = tok_data_struct.nfl;   % number of flux loops
nlv = tok_data_struct.nlv;   % number of voltage loops

%...................
% Mutual inductances

Mcc = tok_data_struct.mcc;   % mutuals coil-to-coil
Mcv = tok_data_struct.mcv;   % mutuals coil-to-vessel
Mvv = tok_data_struct.mvv;   % mutuals vessel-to-vessel

Mpc = tok_data_struct.mpc;   % mutuals plasma-to-coils
Mpv = tok_data_struct.mpv;   % mutuals plasma-to-vessel

Mlc = tok_data_struct.mlc;   % mutuals floops-to-coils
Mlv = tok_data_struct.mlv;   % mutuals floops-to-vessel

Mhc = tok_data_struct.mhc;   % mutuals vloops-to-coils
Mhv = tok_data_struct.mhv;   % mutuals vloops-to-vessel

%......................
% Conductor resistances

Resc = tok_data_struct.resc; % coil   resistances
Resv = tok_data_struct.resv; % vessel resistances

%.........................
% Output names and signals

flnames = tok_data_struct.flnames;
lvnames = tok_data_struct.lvnames;

%flsignals = tok_data_struct.flsignals; % flux loop diagnostic signals

%--------------------------------------------------------------------------
% Compute a linear plasma response
%--------------------------------------------------------------------------

switch plasma_resp
    
    case 0 % vacuum model
        
        response = [];
        
        nc = length(find(cccirc));
        
    case 1 % rigid r-response  (rzrigid)
        
        options.cccirc    = cccirc;
        options.iradial   = 1;
        options.ivertical = 0;
        
        response = rzrigid(eq, tok_data_struct, options);
        
    case 2 % rigid z-response  (rzrigid)
        
        options.cccirc    = cccirc;
        options.iradial   = 0;
        options.ivertical = 1;
        
        response = rzrigid(eq, tok_data_struct, options);
    
    case 3 % rigid rz-response (rzrigid)
        
        options.cccirc    = cccirc;
        options.iradial   = 1;
        options.ivertical = 1;
        
        response = rzrigid(eq, tok_data_struct, options);
    
    case 4 % nonrigid response (gspert)
        
end

if plasma_resp == 0
    
    djphiDIc = zeros(nz*nr,nc);
    djphiDIv = zeros(nz*nr,nv);
    
    Xcc = zeros(nc,nc);
    Xcv = zeros(nc,nv);
    Xvc = zeros(nv,nc);
    Xvv = zeros(nv,nv);
    Xpc = zeros(nz*nr,nc);
    Xpv = zeros(nz*nr,nv);
    Xfc = zeros(nfl,nc);
    Xfv = zeros(nfl,nv);
    
else
    
    djphiDIc = response.djphiDIc;
    djphiDIv = response.djphiDIv;
    
    Xcc = response.Xcc;
    Xcv = response.Xcv;
    Xvc = response.Xvc;
    Xvv = response.Xvv;
    Xpc = response.Xpc;
    Xpv = response.Xpv;
    Xfc = response.Xfc;
    Xfv = response.Xfv;
    
end

%--------------------------------------------------------------------------
% Define circuit connections
%--------------------------------------------------------------------------

% Connect the coil conductors and passive conductors into lumped circuits 
% (reduce the dimension of the state vector)

%..........................................................
% Remove open-circuited coils and selected vessel groupings

% convention used is to index the conductors to keep

% coils

idx = find(cccirc ~= 0);

cccirc = cccirc(idx);
nc = length(cccirc);

Resc = Resc(idx);

Mcc = Mcc(idx,idx);
Mcv = Mcv(idx,:);
Mpc = Mpc(:,idx);
Mlc = Mlc(:,idx);

% vessel

idx1 = find(vvcirc ~= 0);

vvcirc = vvcirc(idx1);

idx2 = ismember(vvgroup,idx1);
idx3 = find(idx2);

nv = length(idx3);

Resv = Resv(idx3);

Mvv = Mvv(idx3,idx3);
Mcv = Mcv(:,idx3);
Mpv = Mpv(:,idx3);
Mlv = Mlv(:,idx3);

djphiDIv = djphiDIv(:,idx3);

Xcv = Xcv(:,idx3);
Xvc = Xvc(idx3,:);
Xvv = Xvv(idx3,idx3);
Xpv = Xpv(:,idx3);
Xfv = Xfv(:,idx3);

% response of plasma current density for all (un-connected) conds

djphiDI = [djphiDIc djphiDIv];

% mutual inductance and resistance matrices for all (un-connected) conds

M = [(Mcc+Xcc) (Mcv+Xcv); (Mcv'+Xvc) (Mvv+Xvv)];

R = diag([Resc; Resv]);

MpcMpv = [Mpc Mpv];
XpcXpv = [Xpc Xpv];

% mutual inductance matrices from (un-connected) conds to diagnostics

Mlc = Mlc + Xfc;
Mlv = Mlv + Xfv;

%...................................................................
% Define matrices that map from connected to un-connected conductors

% COILS

ncx = max(cccirc); % number of coil circuits
Pcc = zeros(nc,ncx);

for ii = 1:ncx
    idx = find(cccirc == ii);
    Pcc(idx,ii) = ones(length(idx),1);
end
    
% VESSEL

nvx = max(vvcirc); % number of vessel circuits

if isempty(nvx)
    nvx = 0; 
end
    
Pvv = zeros(nv,nvx);

% Split currents as if vessel elements were connected in parallel

for ii = 1:nvx    
    idx = find(vvgroup == ii);     
    Resv_parallel = 1/sum(1./Resv(idx)); % parallel resistance of vvcircuit
    vvfrac = Resv_parallel./Resv(idx);
    Pvv(idx,ii) = vvfrac;
end
        
% TOTAL

ncvx = ncx + nvx;      % total number of circuits after connections
Pxx  = [Pcc zeros(nc,nvx); zeros(nv,ncx) Pvv];

%...................................................................
% Compute mutual inductance and resistance matrices for the circuits

Mxx = Pxx'*(M)*Pxx;
Rxx = Pxx'*(R)*Pxx;

%--------------------------------------------------------------------------
% Build the dynamics equation: xdot = A*x + B*v + F*w
%--------------------------------------------------------------------------

% Build the dynamics matrix A

amat = -inv(Mxx)*Rxx;

% Build the input voltage matrix B

vmat = [eye(ncx); zeros(nvx,ncx)];

bmat = Mxx\vmat;

% Build the disturbance matrix F (betap, li, noise)

fmat = zeros(ncvx,2);

%--------------------------------------------------------------------------
% Build the output equation: y = C*I + D*v + H*w
%--------------------------------------------------------------------------

% Build the output matrix C

cmat = [pinv(Pcc) zeros(ncx,nv); (Mlc+Xfc) (Mlv+Xfv)];

cmat = cmat*Pxx;

% Build the feedforward matrix D

dmat = zeros(ncx+nfl,ncx);

% Build the disturbance matrix H

hmat = zeros(size(cmat,1),2);

%--------------------------------------------------------------------------
% Make structure tok_system
%--------------------------------------------------------------------------

%...............................................
% Define character array of state variable names

if exist('ccxnames', 'var')
    states = ccxnames;
else
    cstates = cell(ncx,1);
    for ii = 1:ncx
        cstates{ii} = strcat('cc', int2str(ii));
    end
    states = char(cstates);
end

if nvx ~= 0
    if exist('vvxnames', 'var')
        states = char(states, vvxnames);
    else
        vstates = cell(nvx,1);
        for ii = 1:nvx
            vstates{ii} = strcat('vv', int2str(ii));
        end
        states = char(states, char(vstates));
    end
end

%................................
% Define output names and signals

outputs = char(states(1:ncx,:), char(flnames));

%...........................................
% Make structure with descriptions of fields

descriptions = struct( ...
    'ncx',     'number of coil circuits',                               ...
    'nvx',     'number of vessel circuits',                             ...
    'amat',    'state matrix for state vector I(ncx+nvx)',              ...
    'bmat',    'input matrix for coil voltages',                        ...
    'cmat',    'output matrix with outputs cc',                         ...
    'dmat',    'feedforward matrix',                                    ...
    'states',  'state variable names',                                  ...
    'outputs', 'output names',                                          ...
    'Pcc',     'map coil circuits to individual conductors',            ...
    'Pvv',     'map vessel circuits to indiviudal conductors',          ...
    'Pxx',     'map all circuits to individual conductors',             ...
    'M',       'mutual inductance matrix (unconnected conds)',          ...
    'R',       'resistance matrix (unconnected conds)',                 ...
    'Mxx',     'mutual inductance matrix (connected circuits)',         ...
    'Rxx',     'resistance matrix (connected circuits)',                ...
    'djphiDI', 'current density response matrix  (unconnected conds)',  ...
    'MpcMpv',  'mutual inductance matrix on grid (unconnected conds)',  ...
    'XpcXpv',  'linearized plasma response matrix for flux on grid',    ...
    'response','plasma response model objects (rzrigid or gspert)'      ...
);

tok_system = struct( ...
    'tokamak',        tok_data_struct.tokamak, ...
    'ncx',            ncx,                     ...
    'nvx',            nvx,                     ...
    'amat',           amat,                    ...
    'bmat',           bmat,                    ...
    'cmat',           cmat,                    ...
    'dmat',           dmat,                    ...
    'fmat',           fmat,                    ...
    'hmat',           hmat,                    ...
    'states',         states,                  ...
    'outputs',        outputs,                 ...
    'Pcc',            Pcc,                     ...
    'Pvv',            Pvv,                     ...
    'Pxx',            Pxx,                     ...
    'M',              M,                       ...
    'R',              R,                       ...      
    'Mxx',            Mxx,                     ...
    'Rxx',            Rxx,                     ...
    'djphiDI',        djphiDI,                 ...
    'MpcMpv',         MpcMpv,                  ...
    'XpcXpv',         XpcXpv,                  ...
    'response',       response,                ...
    'descriptions',   descriptions,            ...
    'build_inputs',   build_inputs             ...
);
 
end
