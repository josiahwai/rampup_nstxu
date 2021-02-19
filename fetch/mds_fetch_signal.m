function signal = mds_fetch_signal(shot, tree, tag, plotit)

if ~exist('plotit','var'), plotit = 0; end

mdshost = 'skylark.pppl.gov:8501';
mdsconnect(mdshost);
mdsopen(tree, shot);
times = mdsvalue(strcat('dim_of(', tag, ')'));
sigs = mdsvalue(tag);

signal = variables2struct(shot,tag,times,sigs);

mdsclose;
mdsdisconnect;

if plotit
  figure
  plot(times, sigs)
end