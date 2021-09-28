MDSPLUS_DIR=getenv('MDSPLUS_DIR');
addpath(genpath(MDSPLUS_DIR))

%...................................
% Desired shot number, tree, and tag

shot = 204069;
tree = 'efit01';
tag = '.RESULTS.AEQDSK:WMHD'; %LI, IPMEAS

% tree = 'WF';
% tag = '.PNB';


%.........................................
% Open connection to NSTX-U MDSplus server

mdshost = '127.0.0.1:8000';
mdsconnect(mdshost)

%..........
% Open tree

mdsopen(tree, shot)

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
