%...................................
% Desired shot number, tree, and tag
close all
shot = 204660;

MDSPLUS=getenv('MDSPLUS_DIR');
addpath(genpath(MDSPLUS))
mdshost = 'skylark.pppl.gov:8501';
mdsclose;
mdsdisconnect;
mdsconnect(mdshost);

tree = 'ENGINEERING';
mdsopen(tree, shot);
tag_prefix = '.EPICS.FCPC.DIGITIZERS:';
tag_suffix =  {'PF1AU_P1SV', 'PF1AU_P2SV'};
               
tags = strcat(tag_prefix, tag_suffix);



figure
hold on
for i = 1:length(tags)
  tag = tags{i};
  mdsTimes = mdsvalue(strcat('dim_of(', tag, ')'));
  mdsTagValue = mdsvalue(tag);
  plot(mdsTimes, smooth(mdsTagValue, 100), 'linewidth', 2)
end
xlim([0 0.3])



% SIMULATION
xnames = {'OH', 'PF1AU', 'PF2U', 'PF3U', 'PF5', 'PF3L', 'PF2L', 'PF1AL'};
load('sim_inputs204660.mat')

t = sim_inputs.tspan(2:end);
v = sim_inputs.Uhat(2,:);

plot(t, sim_inputs.Uhat', 'color', [1 1 1] * 0.8)
plot(t, v, 'linewidth', 2, 'color', 'r')
legend(tag_suffix, 'fontsize', 12)




























