function lvsignals = get_lvsignals(shot, times, tree, tags)
%
% GET_LVSIGNALS
%
% SYNTAX: get_lvsignals
% 
% PURPOSE: Read loop voltage signals from MDSplus.
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
%                      operations tree
%
%   tags...............cell array containing MDSplus tag names for loop 
%                      voltage signals
%
% OUTPUTS:
%
%   lvsignals..........structure with the following fields
%
%                         shot:   NSTX-U shot number
%
%                         tags:   MDSplus tag names for loop voltages
%
%                         times:  array of times
%
%                         sigs:   matrix of lvsignals with dimensions 
%                                 [length(times), length(tags)]
%
% AUTHOR: Patrick J. Vail
%
% DATE: 06/08/17
%
% MODIFICATION HISTORY
%   Patrick J. Vail: Original File 06/08/17
%
%..........................................................................

if nargin == 0
    error('ERROR get_lvsignals: must specify NSTX-U shot number')
end

if isempty(tree) && shot >= 200000
    tree = 'operations';
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

mdsTimes = mdsvalue(strcat('dim_of(', tags{1}, ')'));

% Find times between tstart and tend

if ~isempty(times)
    
    tstart = times(1);
    tend   = times(2);
    
    idx = intersect(find(mdsTimes >= tstart), find(mdsTimes <= tend));
    mdsTimes = mdsTimes(idx);
    
end

% Load the signals

mdsSigs = zeros(length(mdsTimes), length(tags));

for ii = 1:length(tags)
    
    tarray = strcat(num2str(mdsTimes(1)),':', num2str(mdsTimes(end)));
    
    try
        mdsSigs(:,ii) = mdsvalue(strcat(tags{ii}, '[', tarray, ']'));
    catch
    end
        
end

mdsdisconnect;

%................
% Output the data

lvsignals = struct(      ...
    'shot',  shot,       ...
    'tags',  char(tags), ...
    'times', mdsTimes,   ...
    'sigs',  mdsSigs     ...
);

end
