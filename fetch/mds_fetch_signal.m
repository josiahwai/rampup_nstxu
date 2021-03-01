% use times=[] to fetch all available times

function signal = mds_fetch_signal(shot, tree, times, tag, plotit)

if ~exist('plotit','var'), plotit = 0; end

mdshost = 'skylark.pppl.gov:8501';
mdsconnect(mdshost);
mdsopen(tree, shot);
mdstimes = mdsvalue(strcat('dim_of(', tag, ')'));
mdssigs = mdsvalue(tag);


if isempty(times)
  times = mdstimes;
  sigs = mdssigs;
else
  if length(mdstimes) ~= size(mdssigs,1)
    mdssigs = mdssigs';
  end
  sigs = interp1(mdstimes, mdssigs, times)';
end


signal = variables2struct(shot,tag,times,sigs);

mdsclose;
mdsdisconnect;

if plotit
  figure
  plot(mdstimes, mdssigs)
  if length(times) == 1
    xline(times, '--k', 'linewidth', 2);
  end
end