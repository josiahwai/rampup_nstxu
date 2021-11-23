%...................................
% Desired shot number, tree, and tag

shot = 204118;
tree = 'EFIT01';
tag = '.RESULTS.AEQDSK:IPMEAS'; %LI, IPMEAS

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
plot(mdsTimes,mdsTagValue)
xlabel('Time [s]')
ylabel(tag)

