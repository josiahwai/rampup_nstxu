% To summarize:
%
% PF5 positive current pulls on plasma, decreases outer gap
% PF5 negative current pushes on plasma, increases outer gap
% 
% PF3 positive current pulls on plasma, decreases upper gap, increases elongation
% PF3 negative current pushes on plasma, increases upper gap, decreases elongation
% 
% To increase elongation without changing width, 'a', need to pull with PF3 (+) but this
% increases plasma width so also need to push with PF5 (-)

clear all; clc; close all

load('matlab.mat')
load('eqs.mat')

i = 12;
eq0 = eqs{i};
x0 = pinv(circ.Pxx) * [eq0.ic; eq0.iv; eq0.cpasma];  % initial state

eq = eq0;
x = x0;

profiles.pres = efit01_eqs.gdata(i).pres;
profiles.fpol = efit01_eqs.gdata(i).fpol;

%%
[spec,init,config] = make_gs_inputs(x, profiles, eq0, tok_data_struct);
config.max_iterations = 8;

%   init = efit01_eqs.gdata(i);
spec.weights.sep = ones(size(spec.weights.sep)) * 1;
spec.targets.ic = spec.locks.ic;
spec.weights.ic = ones(size(spec.targets.ic)) * 3e-3;
spec.locks.ic = nan(size(spec.locks.ic));
spec.targets.iv = spec.locks.iv;
spec.weights.iv = ones(size(spec.targets.iv)) * 3e-3;
spec.locks.iv = nan(size(spec.locks.iv));


eq = gsdesign(spec, init, config);

%%
close all
x = x0;
x(6) = x(6) + 300; % positive to pull it out towards coil, negative to push in
x(9) = x(9) + 300; 
% x(8) = x(8) - 100;  

% x([6 9]) = x([6 9]) + [100 -100]';   % moves it up


[spec,init,config] = make_gs_inputs(x, profiles, eq0, tok_data_struct);
config.max_iterations = 8;

%   init = efit01_eqs.gdata(i);
spec.weights.sep = ones(size(spec.weights.sep)) * 1;
spec.targets.ic = spec.locks.ic;
spec.weights.ic = ones(size(spec.targets.ic)) * 3e-3;
spec.locks.ic = nan(size(spec.locks.ic));
spec.targets.iv = spec.locks.iv;
spec.weights.iv = ones(size(spec.targets.iv)) * 3e-3;
spec.locks.iv = nan(size(spec.locks.iv));


eq1 = gsdesign(spec, init, config);


figure
plot_eq(eq, 'r')
plot_eq(eq1, 'k')






















































