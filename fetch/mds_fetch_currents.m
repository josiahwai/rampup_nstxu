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
tags = {'.ANALYSIS:IOH'};

figure
hold on
for i = 1:length(tags)
  tag = tags{i};
  mdsTimes = mdsvalue(strcat('dim_of(', tag, ')'));
  mdsTagValue = mdsvalue(tag);
  plot(mdsTimes, smooth(mdsTagValue, 100))
end
% xlim([0 0.3])
legend(tags, 'fontsize', 12)




























