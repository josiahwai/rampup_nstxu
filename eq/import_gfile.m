% Import gfile and solve for the gsdesign equilibrium that best matches the 
% gfile equilibrium. 

% EXAMPLE
% shot = 204660;
% time_ms = 120;
% eq = import_gfile(shot, time_ms, 0, 1)

function eq = import_gfile(shot, time_ms, saveit, plotit)

if ~exist('plotit', 'var'), plotit = 0; end
if ~exist('saveit', 'var'), saveit = 0; end


% ========
% SETTINGS
% ========
% shotdir = '/u/jwai/rampup_nstxu/eq/geqdsk/';
shotdir = '/Users/jwai/Research/rampup_nstxu/eq/geqdsk/';
savedir = '/Users/jwai/Research/rampup_nstxu/eq/eq/';
load('nstxu_obj_config2016_6565.mat')

% ============================
% DEFINE SPEC, INIT, CONFIG
% ============================

efit_eq = read_eq(shot, time_ms/1000, shotdir);
init = efit_eq.gdata;

config = tok_data_struct;
config.max_iterations = 30;
config.constraints = 1;
config.nkn = 4;
config.no_edge_current = false;
config.no_edge_gradient = false;
config.plot_settings.SOL.n = 10;
config.plot_settings.SOL.d = 1e-3;

spec.targets.rsep = init.rbbbs;
spec.targets.zsep = init.zbbbs;
spec.weights.sep = 1*ones(length(spec.targets.rsep),1);
   
spec.targets.cpasma = init.cpasma;
spec.weights.cpasma = 1;

spec.targets.psibry = init.psibry;
spec.weights.psibry = 1;

config.constraints = 1;  % allow for scaling/peaking of profiles
config.pres0 = init.pres;
config.fpol0 = init.fpol;


gs_configure
gs_initialize
gs_eq_analysis
spec.targets.li = li;
spec.targets.betap = betap;
spec.weights.li = 10;
spec.weights.betap = 10;


spec.cccirc = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 ...
    11 12 13];

spec.limits.ic(1:tok_data_struct.nc,1) = [-20 0 0 -8 0 0 -13 -13 -13 ...
    -13 0 0 0 0 0 0 -24 -24 -24 -24 -13 -13 -13 -13 0 0 -8 0 0]'*1000;

spec.limits.ic(1:tok_data_struct.nc,2) = [20 15 0 13.5 15 15 8 8 8 8 ...
    13 13 13 13 13 13 0 0 0 0 8 8 8 8 15 15 13.5 0 15]'*1000;


% ============================
% MAP COIL CURRENTS TO TOKSYS
% ============================

% import coil currents
load(['coils' num2str(shot) '.mat'])
i = find( floor(coils.t*1000) == time_ms);
if isempty(i)
  warning('Could not import coil currents')
end
ic = coils.ic(i,:)';
iv = coils.iv(i,:)';

% vacuum vessel degrouping matrix
vvgroup = [1  1  2  2  3  4  5  5  5  6  6  6  6  6  6  6  6  6  6  6 ...
     6  6  6  6  7  7  7  7  7  8  8  8  8  9  9  9  9  9  9  9  9  9 ...
     9  9  9  9  9  9  9  9  9  9  9  9  9  9  9 10 10 10 11 11 11 11 ...
    11 12 12 13 13 13 14 15 16 17 18 18 18 19 19 20 20 20 20 20 21 21 ...
    21 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 ...
    22 22 22 23 23 23 23 24 24 24 24 24 25 25 25 25 25 25 25 25 25 25 ...
    25 25 25 25 25 26 26 26 27 28 29 29 30 30 31 31 32 32 33 34 35 36 ...
    37 38 39 40];
  
Pvv = zeros(length(vvgroup), max(vvgroup));  
for i = 1:max(vvgroup)
  Pvv(:,i) = (vvgroup==i) / sum(vvgroup==i);
end
spec.locks.iv = Pvv*iv;
spec.locks.ic = ic(spec.cccirc);

eq = gsdesign(spec, init, config);



% ==============
% PLOT AND SAVE
% ==============
if plotit
  figure
  hold on
  rg = tok_data_struct.rg;
  zg = tok_data_struct.zg;
  plot_nstxu_geo(tok_data_struct)

  [~,cs] = contour(rg,zg,init.psizr,30,'--r', 'linewidth', 0.5);
  contour(rg,zg,eq.psizr,cs.LevelList,'k', 'linewidth', 0.5);
  set(gcf,'Position',[204 38 312 533])
end

if saveit
  fn = [savedir 'eq' num2str(shot) '_' num2str(time_ms) '.mat'];
  save(fn, 'eq');
end


end















