function ivsignals = get_ivsignals(shot, plotit)

if ~exist('plotit','var'), plotit = 0; end

tree = 'efit01';

mdshost = 'skylark.pppl.gov:8501';
mdsconnect(mdshost);
mdsopen(tree, shot);
tag = '.RESULTS.AEQDSK:CCBRSP';

times = mdsvalue(strcat('dim_of(', tag, ')'));
mdsTagValue = mdsvalue(tag);

nvx = 40;
iv = mdsTagValue(end-nvx+1:end,:)';

ivsignals = struct('shot', shot, 'times', times, 'sigs', iv);

if plotit
  figure
  plot(times, iv)
  xlabel('Time')
  ylabel('Vessel Currents [A]')
end