MDSPLUS=getenv('MDSPLUS_DIR');
addpath(genpath(MDSPLUS))

%...................................
% Desired shot number, tree, and tag

shot = 204660;
tree = 'efit01';
tag = '.RESULTS.AEQDSK:ECCURT'; %LI, IPMEAS

tree = 'WF';
tag = '.PNB';


%.........................................
% Open connection to NSTX-U MDSplus server

mdshost = 'skylark.pppl.gov:8501';
mdsconnect(mdshost);

%..........
% Open tree

mdsopen(tree, shot);

%.............
% Get the data

mdsTimes = mdsvalue(strcat('dim_of(', tag, ')'));
mdsTagValue = mdsvalue(tag);

%......
% Plot
figure(2)
clf
plot(mdsTimes,mdsTagValue)
hold on
plot(mdsTimes, smoothdata(mdsTagValue, 'movmedian', 5), 'linewidth', 2)
xlabel('Time [s]')
ylabel(tag)
xlim([0 1])
