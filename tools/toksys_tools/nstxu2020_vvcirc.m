function [vvgroup, vvcirc, vvxnames, Pvv] = nstxu2020_vvcirc

% NSTXU_VVCIRC
%
%   Output the grouping vector, connection vector, and circuit names array 
%   for the desired vessel configuration.
%
% USAGE: nstxu_vvcirc.m
%
% OUTPUTS: 
%
%   vvgroup.............grouping vector combining individual vessel
%                       elements into grouped conductors
%
%   vvcirc..............connection vector for the vessel groups
%
%   vvxnames............character array containing vessel circuit names
%                       
% AUTHOR: Patrick J. Vail
%
% DATE: 11/14/2017
%
% MODIFICATION HISTORY:
%   Patrick J. Vail: Original File 11/14/2017
%   JOSIAH WAI: 40-element version for 2020 version of nstxu
%   tok_data_struct
%
%..........................................................................

% Connection vector and circuit names with all vessel circuits included

vvgroup = [1  1  2  2  3  4  5  5  5  6  6  6  6  6  6  6  6  6  6  6 ...
     6  6  6  6  7  7  7  7  7  8  8  8  8  9  9  9  9  9  9  9  9  9 ...
     9  9  9  9  9  9  9  9  9  9  9  9  9  9  9 10 10 10 11 11 11 11 ...
    11 12 12 13 13 13 14 15 16 17 18 18 18 19 19 20 20 20 20 20 21 21 ...
    21 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 22 ...
    22 22 22 23 23 23 23 24 24 24 24 24 25 25 25 25 25 25 25 25 25 25 ...
    25 25 25 25 25 26 26 26 27 28 29 29 30 30 31 31 32 32 33 34 35 36 ...
    37 38 39 40];
    
vvcirc = 1:40;


Pvv = zeros(length(vvgroup), max(vvgroup));  
for i = 1:max(vvgroup)
  Pvv(:,i) = (vvgroup==i) / sum(vvgroup==i);
end

vvxnames = char({  ...
    'VS1U',   ...
    'VS2U',   ...
    'VS3U',   ...
    'VS4U',   ...
    'VS5U',   ...
    'VS6U',   ...
    'VS7U',   ...
    'VS8U',   ...
    'VS9U',   ...
    'VS10U',  ...
    'VS11U',  ...
    'VS12U',  ...
    'VS13U',  ...
    'VS14U',  ...
    'VS15U',  ...
    'VS15L',  ...
    'VS14L',  ...
    'VS13L',  ...
    'VS12L',  ...
    'VS11L',  ...
    'VS10L',  ...
    'VS9L',   ...
    'VS8L',   ...
    'VS7L',   ...
    'VS6L',   ...
    'VS5L',   ...
    'VS4L',   ...
    'VS3L',   ...
    'VS2L',   ...
    'VS1L',   ...
    'DPU1',   ...
    'DPL1',   ...
    'PPSIUU', ...
    'PPSIUL', ...
    'PPPOUU', ...
    'PPPOUL', ...
    'PPPOLU', ...
    'PPPOLL', ...
    'PPSILL', ...
    'PPSILU'
});




























