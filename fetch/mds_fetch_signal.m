function signal = mds_fetch_signal(shot, tree, tag)


mdshost = 'skylark.pppl.gov:8501';
mdsconnect(mdshost);
mdsopen(tree, shot);
times = mdsvalue(strcat('dim_of(', tag, ')'));
sigs = mdsvalue(tag);

signal = variables2struct(shot,tag,times,sigs);

mdsclose;
mdsdisconnect;

