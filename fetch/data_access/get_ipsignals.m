function ipsignals = get_ipsignals(shot, times, tree, tags)
%
% GET_IPSIGNALS 
%
% SYNTAX: get_ipsignals
% 
% PURPOSE: Read Ip signals from MDSplus.
% 
% INPUTS:
%
%   shot...............NSTX-U shot number
%
%   times..............time interval with format [tstart tend]; if all
%                      available times are desired, then input []
%
%   tree...............MDSplus tree from which to obtain data; if no tree
%                      is specified and shot >= 200000, then use the 
%                      engineering tree
%
%   tags...............cell array containing MDSplus tag names for ip
%                      signals
%
%
% OUTPUTS:
%
%   ipsignals..........structure with the following fields
%
%                         shot:   NSTX-U shot number
%
%                         nodes:  MDSplus node names for ip
%
%                         times:  array of times
%
%                         sigs:   matrix of waveforms with dimensions 
%                                 [length(times), length(nodes)]
%
% AUTHOR: Dan Boyer
%
% DATE: 04/14/20
%
%..........................................................................

if nargin == 0
    error('ERROR get_ipsignals: must specify NSTX-U shot number')
end

if isempty(times)
end

if isempty(tree) && shot >= 200000
    tree = 'engineering';
end

%.........................................
% Open connection to NSTX-U MDSplus server

mdshost = 'skylark.pppl.gov:8501';
mdsconnect(mdshost);

%..........
% Open tree

mdsopen(tree, shot);

%.............
% Get the data


mdsTimes = double(mdsvalue(['dim_of(' tags{1} ')']));
mdsSigs = zeros(length(mdsTimes),length(tags));
for i=1:length(tags)
    mdsSigs(:,i) = double(mdsvalue(tags{i}));
end

mdsdisconnect;

%................
% Output the data

ipsignals = struct(       ...
    'shot',  shot,        ...
    'times', mdsTimes,    ...
    'sigs',  mdsSigs,     ...
    'nodes', {tags}        ...
);

end
