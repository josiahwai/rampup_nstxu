function icsignals = get_vobjcsignals(shot, times, tree, include_coils)
%
% GET_VOBJCSIGNALS 
%
% SYNTAX: get_vobjcsignals
% 
% PURPOSE: Read Vobjective from MDSplus
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
%   vobjsignals..........structure with the following fields
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
% AUTHOR: Dan Boyer
%
% DATE: 09/07/19
%
% MODIFICATION HISTORY
%   Dan Boyer: Original File 09/07/19
%
%..........................................................................

if nargin == 0
    error('ERROR get_vobjsignals: must specify NSTX-U shot number')
end

if isempty(times)
end

if isempty(tree) && shot >= 200000
    tree = 'engineering';
end

%...............................
% Select nodes for fetching data

alpha1_nodes = {'RTF_A1_OH',  'RTF_A1_PF1AU', ...
         'RTF_A1_PF1BU', 'RTF_A1_PF1CU', ...
         'RTF_A1_PF2U',  'RTF_A1_PF3U',  ...
         'RTF_A1_PF4',   'RTF_A1_PF5',   ...
         'RTF_A1_PF3L',  'RTF_A1_PF2L',  ...
         'RTF_A1_PF1CL', 'RTF_A1_PF1BL', ...
         'RTF_A1_PF1AL'                ...
};
alpha2_nodes = {'RTF_A2_OH',  'RTF_A2_PF1AU', ...
         'RTF_A2_PF1BU', 'RTF_A2_PF1CU', ...
         'RTF_A2_PF2U',  'RTF_A2_PF3U',  ...
         'RTF_A2_PF4',   'RTF_A2_PF5',   ...
         'RTF_A2_PF3L',  'RTF_A2_PF2L',  ...
         'RTF_A2_PF1CL', 'RTF_A2_PF1BL', ...
         'RTF_A2_PF1AL'                ...
};
cv1_nodes = {'RTC_CV1OH',  'RTC_CV1PF1AU', ...
         'RTC_CV1PF1BU', 'RTC_CV1PF1CU', ...
         'RTC_CV1PF2U',  'RTC_CV1PF3U',  ...
         'RTC_CV1PF4',   'RTC_CV1PF5',   ...
         'RTC_CV1PF3L',  'RTC_CV1PF2L',  ...
         'RTC_CV1PF1CL', 'RTC_CV1PF1BL', ...
         'RTC_CV1PF1AL'                ...
};
cv2_nodes = {'RTC_CV2OH',  'RTC_CV2PF1AU', ...
         'RTC_CV2PF1BU', 'RTC_CV2PF1CU', ...
         'RTC_CV2PF2U',  'RTC_CV2PF3U',  ...
         'RTC_CV2PF4',   'RTC_CV2PF5',   ...
         'RTC_CV2PF3L',  'RTC_CV2PF2L',  ...
         'RTC_CV2PF1CL', 'RTC_CV2PF1BL', ...
         'RTC_CV2PF1AL'                ...
};
polarity_nodes = {'SETUP.RECTIFIERS.OH.POLARITY', 
    'SETUP.RECTIFIERS.PF1AU.POLARITY',
    'SETUP.RECTIFIERS.PF1BU.POLARITY',
    'SETUP.RECTIFIERS.PF1CU.POLARITY',
    'SETUP.RECTIFIERS.PF2U.POLARITY',
    'SETUP.RECTIFIERS.PF3U.POLARITY',
    'SETUP.RECTIFIERS.PF4.POLARITY',
    'SETUP.RECTIFIERS.PF5.POLARITY',
    'SETUP.RECTIFIERS.PF3L.POLARITY',
    'SETUP.RECTIFIERS.PF2L.POLARITY',
    'SETUP.RECTIFIERS.PF1CL.POLARITY',
    'SETUP.RECTIFIERS.PF1BL.POLARITY',
    'SETUP.RECTIFIERS.PF1AL.POLARITY'};

nseries_nodes = {'SETUP.RECTIFIERS.OH.NSERIES', 
    'SETUP.RECTIFIERS.PF1AU.NSERIES',
    'SETUP.RECTIFIERS.PF1BU.NSERIES',
    'SETUP.RECTIFIERS.PF1CU.NSERIES',
    'SETUP.RECTIFIERS.PF2U.NSERIES',
    'SETUP.RECTIFIERS.PF3U.NSERIES',
    'SETUP.RECTIFIERS.PF4.NSERIES',
    'SETUP.RECTIFIERS.PF5.NSERIES',
    'SETUP.RECTIFIERS.PF3L.NSERIES',
    'SETUP.RECTIFIERS.PF2L.NSERIES',
    'SETUP.RECTIFIERS.PF1CL.NSERIES',
    'SETUP.RECTIFIERS.PF1BL.NSERIES',
    'SETUP.RECTIFIERS.PF1AL.NSERIES'};

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

alpha1_nodes = alpha1_nodes(keep_idx);
alpha2_nodes = alpha2_nodes(keep_idx);
cv1_nodes = cv1_nodes(keep_idx);
cv2_nodes = cv2_nodes(keep_idx);
polarity_nodes = polarity_nodes(keep_idx);
nseries_nodes = nseries_nodes(keep_idx);

%.........................................
% Open connection to NSTX-U MDSplus server

mdshost = 'skylark.pppl.gov:8501';
mdsconnect(mdshost);

%..........
% Open tree

mdsopen(tree, shot);

%.............
% Get the data

mdspath  = '.PPCC.PCS.RT:';
mdsTimes = mdsvalue(strcat('dim_of(', mdspath, char(alpha1_nodes(1)), ')'));

% Load the signals

alpha1Sigs = zeros(length(mdsTimes), length(alpha1_nodes));
alpha2Sigs = zeros(length(mdsTimes), length(alpha2_nodes));
cv1Sigs = zeros(length(mdsTimes), length(cv1_nodes));
cv2Sigs = zeros(length(mdsTimes), length(cv2_nodes));

for ii = 1:length(alpha1_nodes)   
    alpha1Sigs(:,ii) = mdsvalue(strcat(mdspath, char(alpha1_nodes(ii))));
    alpha2Sigs(:,ii) = mdsvalue(strcat(mdspath, char(alpha2_nodes(ii))));
    cv1Sigs(:,ii) = mdsvalue(strcat(mdspath, char(cv1_nodes(ii))));
    cv2Sigs(:,ii) = mdsvalue(strcat(mdspath, char(cv2_nodes(ii))));
    polarity = mdsvalue(strcat('.PPCC.PCS.RT.', char(polarity_nodes(ii))));
    nseries = mdsvalue(strcat('.PPCC.PCS.RT.', char(nseries_nodes(ii))));
    vcoil0 = 0*alpha1Sigs(:,ii);
    for jj=1:length(mdsTimes)
        if cv1Sigs(jj,ii) == 1
            vcoil0(jj) = cos(alpha1Sigs(jj,ii)*pi/180);
        elseif cv2Sigs(jj,ii) == 1
            vcoil0(jj) = -cos(alpha2Sigs(jj,ii)*pi/180);
        end
    end
    mdsSigs(:,ii) = double(polarity*int32(nseries)*1012)*double(vcoil0);
end
    
mdsdisconnect;

%................
% Output the data

icsignals = struct(       ...
    'shot',  shot,        ...
    'times', mdsTimes,    ...
    'sigs',  mdsSigs      ...
);

end
