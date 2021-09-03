% Inputs: eq
function plot_eq(eq, color)

if ~exist('color','var'), color = 'r'; end
if isfield(eq,'gdata'), eq = eq.gdata; end

% Load tok_data_struct
load('nstxu_obj_config2020_6565.mat')
rg = tok_data_struct.rg;
zg = tok_data_struct.zg;

plot_nstxu_geo(tok_data_struct)

psizr = eq.psizr;
psibry = eq.psibry;

% time = 1000*eq.time;
% shot = eq.shotnum;

contour(rg,zg,psizr,[psibry psibry], 'color', color, 'linewidth', 1.5);
contour(rg,zg,psizr, 10, 'color', [1 1 1]*0.8, 'linewidth', 0.5);

% set(gcf,'Position',[204 38 312 533])

% title([num2str(shot) ': ' num2str(time) 'ms'])















