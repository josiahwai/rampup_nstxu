
function [spec, config] = make_gsdesign_inputs(x, tok_data_struct, eq, traj, t, tspan)
clear spec init config

init = eq;

config = tok_data_struct;
config.max_iterations = 30;
config.constraints = 1;
config.nkn = 10;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

spec.targets.rsep = init.rbbbs;
spec.targets.zsep = init.zbbbs;
spec.weights.sep = ones(length(spec.targets.rsep),1) * 1e-3;
   
spec.targets.psibry = init.psibry;
spec.weights.psibry = 1e-3;

gs_configure

spec.targets.li = interp1(tspan, traj.li, t);
spec.weights.li = 10;

spec.targets.betap = interp1(tspan, traj.betap, t);
spec.weights.betap = 10;

config.constraints = 1; % allow for scaling/peaking of profiles
config.pres0 = interp1(tspan, traj.pres, t)';
config.fpol0 = interp1(tspan, traj.fpol, t)';

spec.cccirc = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 ...
    11 12 13];

spec.limits.ic(1:tok_data_struct.nc,1) = [-20 0 0 -8 0 0 -13 -13 -13 ...
    -13 0 0 0 0 0 0 -24 -24 -24 -24 -13 -13 -13 -13 0 0 -8 0 0]'*1000;

spec.limits.ic(1:tok_data_struct.nc,2) = [20 15 0 13.5 15 15 8 8 8 8 ...
    13 13 13 13 13 13 0 0 0 0 8 8 8 8 15 15 13.5 0 15]'*1000;

% inner wall limited, specify a target for the boundary point so that solver 
% does not go vertically unstable
% if t*1000 < 200
%   spec.targets.rbdef = 0.315;
%   spec.targets.zbdef = 0.02;
%   spec.weights.bdef = 100;
% end
     
% indices in x for: coil currents, vessel currents, plasma current
iic = 1:8;
iiv = 9:48;
iip = 49;
  
% Lock vessel current
[~,~,~,Pvv] = nstxu2020_vvcirc;
iv = x(iiv);
spec.locks.iv = Pvv*iv;

% Lock plasma current
spec.locks.cpasma = x(iip);

% Lock coil currents
% Some coil currents were turned off for campaign, so lock these to zero

removed_coils = {'PF1BU', 'PF1BL', 'PF1CU', 'PF1CL', 'PF4'};  
coils = {'OH', 'PF1AU', 'PF1BU', 'PF1CU', 'PF2U', 'PF3U', 'PF4', 'PF5', ...
  'PF3L', 'PF2L', 'PF1CL', 'PF1BL', 'PF1AL'};

for i = 1:length(coils)
  iuse(i) = ~ismember(coils{i}, removed_coils);
end
ic = zeros(size(coils));
ic(iuse) = x(iic);
spec.locks.ic = ic(spec.cccirc);

















