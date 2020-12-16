  function  eq= eq_mds(shot, tree, server, toupper, verbose)
%
%  SYNTAX:
%         eq= eq_mds(shot, tree, server, toupper, verbose); % full call
%         eq= eq_mds(shot);                        % defaults to DIII-D server
%         eq= eq_mds(shot, 'EFIT01', 'NSTX');      % for NSTX
%         eq= eq_mds(shot, 'NB', 'DIII-D');        % Other MDS+ Tree Data -beam
%  Should also Run:
%         eq= eq_mk_time_last(eq); % puts time at end (i.e. (130,1) => (1,130)
%
%  PURPOSE: Get EFIT equilibrium (GEQDSK AEQDSK MEQDSK) from mdsplus database. 
%
%  INPUT: <default>
%    shot   = shot number
%    tree   = which efit tree to use <'EFIT01'>
%    server = which MDS+ database to use: defaults to <'DIII-D'>
%             'DIII-D' opens 'atlas.gat.com'
%             'NSTX'   opens 'europa.pppl.gov:8501'
%             'OPEN'   assumes mdsconnect already called
%              server  arbitrary server is called using mdsconnect(server)
%    toupper= 1= all variables made upper case, =-1 all var. made lower case
%             [0]= no change, variables made depending on mds case (typical UC)
%    verbose = set to 1 to get diagnostic prints during execution
%
%  OUTPUT:
%    eq    = structure containing EFIT eqdsk variables for all times
%            structure follows MDS tree structure with 'TOP' replaced with 'eq'
%            See: eq.allnames
%
%  EXAMPLE Variables (Examples below are for DIII-D -- NSTX is different)
%
%   Ex:  eq.RESULTS.AEQDSK.IPMEAS, eq.RESULTS.GEQDSK.PSIRZ 
%   Use: plot(eq.RESULTS.GEQDSK.GTIME,eq.RESULTS.GEQDSK.CPASMA)
%        contour(eq.RESULTS.GEQDSK.R(:,1), eq.RESULTS.GEQDSK.Z(:,1), ...
%        eq.RESULTS.GEQDSK.PSIRZ(:,:,round(length(eq.RESULTS.GEQDSK.GTIME)/2))')
%
%    Some Extra items added to structure (all lower case)
%      eq.id    =   sting array of important data identifyer enf
%      eq.shot  =   shot number
%      eq.server=   MDS+ server
%      eq.allnames= list of all variables in structure with full struc. address
%      eq.mdsnames= list of all variables in structure with mds struc. address
% 
%  NOTE: 1)Function returns all TIME data for shot in large arrays.
%          Example:
%          eq.RESULTS.AEQDSK.BETAN(255,1)
%          eq.RESULTS.AEQDSK.ATIME(255,1); % Time data vector for AEQDSK data
%          plot(eq.RESULTS.AEQDSK.ATIME,eq.RESULTS.AEQDSK.BETAN)
%          eq.RESULTS.GEQDSK.GTIME(130,1)
%          eq.RESULTS.GEQDSK.PPRIME(65,130)
%          eq.RESULTS.GEQDSK.PSIRZ(65,65,130) 
%        2) Time is in vector eq...GTIME and is in ms. All other units are as
%           they come from EFIT. psi(Vs/r) I(A) R(m) ...
%
%  SEE ALSO: eq= eq_mk_time_last(eq); (should be run on eq)
%            eq_time_lim eq_ga_env str_to_ws
%
%  NOTE: Written very generally and should return any subtree of MDS+
%        Example: ions= eq_mds(shot,'IONS');

%  WRITTEN BY:  Jim Leuer    ON      3/1/05
%  taken from mds_eq.m uses sub-structure to store
%  USES:   eq_mod
% To see MDS structure on HYDRA run traverser
% tested on DIII-D and NSTX data and should work for JET data but not tested
% ==========================================================================
%  @(#)eq_mds.m	1.7 06/23/10

wait('WARNING eq_mds: eq_mds is OBSOLETE - use get_mds_tree instead')

if(nargin<2)
   wait('ERROR eq_mds: requires at least 2 inputs')
elseif(nargin==2)
   eq = get_mds_tree(shot, tree);
elseif(nargin==3)
   eq = get_mds_tree(shot, tree, server);
elseif(nargin==4)
   eq = get_mds_tree(shot, tree, server, toupper);
elseif(nargin==5)
   eq = get_mds_tree(shot, tree, server, toupper, verbose);
end
