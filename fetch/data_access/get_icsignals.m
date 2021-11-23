function icsignals = get_icsignals(shot, times, tree, include_coils)
%
% GET_ICSIGNALS 
%
% SYNTAX: get_icsignals
% 
% PURPOSE: Read auctioneered coil currents from MDSplus.
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
%   include_coils......cell array containing names of coil circuits for
%                      which to retrieve currents
%
% OUTPUTS:
%
%   icsignals..........structure with the following fields
%
%                         shot:   NSTX-U shot number
%
%                         nodes:  MDSplus node names for coils
%
%                         times:  array of times
%
%                         sigs:   matrix of waveforms with dimensions 
%                                 [length(times), length(nodes)]
%
% AUTHOR: Patrick J. Vail
%
% DATE: 06/09/17
%
% MODIFICATION HISTORY
%   Patrick J. Vail: Original File 06/09/17
%
%..........................................................................

if nargin == 0
    error('ERROR get_icsignals: must specify NSTX-U shot number')
end

if isempty(times)
end

if isempty(tree) && shot >= 200000
    tree = 'engineering';
end

%...............................
% Select nodes for fetching data

nodes = {'RA_AUC_OH',  'RA_AUC_1AU', ...
         'RA_AUC_1BU', 'RA_AUC_1CU', ...
         'RA_AUC_2U',  'RA_AUC_3U',  ...
         'RA_AUC_4',   'RA_AUC_5',   ...
         'RA_AUC_3L',  'RA_AUC_2L',  ...
         'RA_AUC_1CL', 'RA_AUC_1BL', ...
         'RA_AUC_1AL'                ...
};

remove_idx = [];

if sum(ismember(include_coils, 'OH')) == 0
    remove_idx = [remove_idx 1];
end

if sum(ismember(include_coils, 'PF1aU')) == 0
    remove_idx = [remove_idx 2];
end

if sum(ismember(include_coils, 'PF1bU')) == 0
    remove_idx = [remove_idx 3];
end

if sum(ismember(include_coils, 'PF1cU')) == 0
    remove_idx = [remove_idx 4];
end

if sum(ismember(include_coils, 'PF2U')) == 0
    remove_idx = [remove_idx 5];
end

if sum(ismember(include_coils, 'PF3U')) == 0
    remove_idx = [remove_idx 6];
end

if sum(ismember(include_coils, 'PF4')) == 0
    remove_idx = [remove_idx 7];
end

if sum(ismember(include_coils, 'PF5')) == 0
    remove_idx = [remove_idx 8];
end

if sum(ismember(include_coils, 'PF3L')) == 0
    remove_idx = [remove_idx 9];
end

if sum(ismember(include_coils, 'PF2L')) == 0
    remove_idx = [remove_idx 10];
end

if sum(ismember(include_coils, 'PF1cL')) == 0
    remove_idx = [remove_idx 11];
end

if sum(ismember(include_coils, 'PF1bL')) == 0
    remove_idx = [remove_idx 12];
end

if sum(ismember(include_coils, 'PF1aL')) == 0
    remove_idx = [remove_idx 13];
end 

keep_idx = setdiff(1:13,remove_idx);

nodes = nodes(keep_idx);

%.........................................
% Open connection to NSTX-U MDSplus server

mdshost = 'skylark.pppl.gov:8501';
mdsconnect(mdshost);

%..........
% Open tree

mdsopen(tree, shot);

%.............
% Get the data

mdspath  = '.PPCC.PCS.RA:';
mdsTimes = mdsvalue(strcat('dim_of(', mdspath, char(nodes(1)), ')'));

% Load the signals

mdsSigs = zeros(length(mdsTimes), length(nodes));

for ii = 1:length(nodes)
    
    mdsSigs(:,ii) = mdsvalue(strcat(mdspath, char(nodes(ii))));
    
end
    
mdsdisconnect;

%................
% Output the data

icsignals = struct(       ...
    'shot',  shot,        ...
    'nodes', char(nodes), ...
    'times', mdsTimes,    ...
    'sigs',  mdsSigs      ...
);

end
