load('nstxu_obj_config2016_6565.mat')
load('/Users/jwai/Research/snowSim/eq/eq_snow204118_509ms.mat')

% eq.turnfc = tok_data_struct.fcnturn';
% eq.fcturn = ones(size(eq.turnfc));
% eq.fcid = [2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 11 12 13] - 1;
% eq.ecid = [1 1 1 1 1 1 1 1];
% eq.ecturn = [112 110 109.5 108.5 108.5 109.5 110 112];


% eq.turnfc = [880 tok_data_struct.fcnturn'];
% eq.fcturn = ones(size(eq.turnfc));
% eq.fcid = [1 2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 11 12 13];
% eq.ecid = [];
% eq.ecturn = [];



fcnturn = tok_data_struct.fcnturn';
Ifrac = zeros(size(fcnturn));
cccirc = [2 3 4 5 5 6 6 6 6 7 7 7 7 7 7 8 8 8 8 9 9 9 9 10 10 11 12 13] - 1;
for icirc = 1:max(cccirc)
  icoils = find(cccirc == icirc);
  Ifrac(icoils) = fcnturn(icoils) / sum(fcnturn(icoils));
end

eq.turnfc = tok_data_struct.fcnturn';
eq.fcturn = Ifrac;
eq.fcid = cccirc;
eq.ecid = [1 1 1 1 1 1 1 1];
eq.ecturn = [112 110 109.5 108.5 108.5 109.5 110 112];



[rzrig_data, cc0] = rzrig(eq,'NSTXU',tok_data_struct, 0, 0, [], 0, 0);































