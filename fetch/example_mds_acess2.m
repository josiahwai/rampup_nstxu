MDSPLUS_DIR=getenv('MDSPLUS_DIR');
addpath(genpath(MDSPLUS_DIR))

shot = 204069;
tree = 'efit01';
tag = '.RESULTS.AEQDSK:VLOOPMHD'; 
times = [];
plotit = 1;

signal = mds_fetch_signal(shot, tree, times, tag, plotit);


