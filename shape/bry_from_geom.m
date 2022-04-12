

clc; close all

load('eq204653_950.mat')


% target boundary
gap_opts.plotit = 1;
gap_opts.use_out_up_lo = 1;

gaps = get_nstxu_gaps(eq, gap_opts);
target.rcp = gaps.r;
target.zcp = gaps.z;












