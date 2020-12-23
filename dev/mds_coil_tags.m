xnames = {'OH', 'PF1AU', 'PF2U', 'PF3U', 'PF5', 'PF3L', 'PF2L', 'PF1AL'};

tree = 'ENGINEERING';
tag_prefix = '.EPICS.FCPC.DIGITIZERS:';
tag_prefix = '.ANALYSIS:';


% ( ) OH:     {'OH_P1S_1PV', 'OH_P1S_2PV', ... 
%              'OH_P2S_1PV', 'OH_P2S_2PV', ...
%              'OH_P3S_1PV', 'OH_P3S_2PV', ...
%              'OH_P4S_1PV', 'OH_P4S_2PV', ...
%              'OH_P5S_1PV', 'OH_P5S_2PV', ...
%              'OH_P6S_1PV', 'OH_P6S_2PV', ...
%              'OH_P7S_1PV', 'OH_P7S_2PV', ...
%              'OH_P8S_1PV', 'OH_P8S_2PV'};
%               && IOH
%               k1 = [1 3 9 11];     % same sign as IOH
%               k2 = [6 8 10 12];    % opposite IOH and k1
%                
% (+) PF1AU: PF1AU_P1SV, PF1AU_P2SV, PF1AU_P2SI && IPF1AU
% (+) PF2U: PF2U_P1SV, PF2U_P2SV,  PF2U_CUR2 && IPF2U
% (-) PF3U: PF3U_P1S_1PV, PF3U_P1S_2PV, PF3U_P2S_1PV, PF3U_P2S_2PV, PF3U_CUR2 && IPF3U 
% (-) PF5: PF5_P1SV, PF5_P2SV, PF5_P3SV, PF5_CUR1 && IPF5
% (-) PF3L: PF3L_P1S_1PV, PF3L_P1S_2PV, PF3L_P2S_1PV, PF3L_P2S_2PV, PF3L_CUR2 && IPF3L
% (+) PF2L: PF2L_P1SV, PF2L_P2SV,  PF2L_CUR2 && IPF2L
% (+) PF1AL: PF1AL_P1SV, PF1AL_P2SV, PF1AL_CUR1 && IPF1AL

% note for PF3U/L: the 2PV voltages are opposite 1PV voltages. sign 2PV goes
%  with sign IPF3U/L













