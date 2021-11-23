% xdot = A*x + B*u
% x = [ic iv]
% u = [vc, icdot, ivdot, ipdot, mcIp, mvIp]
%
% Note that icdot, ivdot, ipdot are the experimental values. 


function [x,y] = vessel_dynamics(t, x, u, rc, rv, lc, lv, mcc_triu, mvv_triu, mcv, file_args)

args = file_args{1};

circ = args.circ;
ts = args.ts;

% extract components from x
ic = x(circ.iicx);
iv = x(circ.iivx);

% extract components from u
vc = u(1:circ.ncx)';
icdot = u(circ.ncx+1:2*circ.ncx)';
ivdot = u(2*circ.ncx+1:2*circ.ncx+circ.nvx)';
ipdot = u(2*circ.ncx+circ.nvx+1:2*circ.ncx+circ.nvx+1)';
mcIp = u(end-circ.nvx-circ.ncx+1:end-circ.nvx)';
mvIp = u(end-circ.nvx+1:end)';

% u = u';
% 
% vc = u(1:circ.ncx);
% u(1:circ.ncx) = [];
% 
% icdot = u(1:circ.ncx);
% u(1:circ.ncx) = [];
% 
% ivdot = u(1:circ.nvx);
% u(1:circ.nvx) = [];
% 
% ipdot = u(1);
% u(1) = [];
% 
% mcIp = u(1:circ.ncx);
% u(1:circ.ncx) = [];
% 
% mvIp = u(1:circ.nvx);
% u(1:circ.nvx) = [];


% reconstruct mutual inductances from parameters
mcc = mcc_triu + mcc_triu' + diag(lc);
mvv = mvv_triu + mvv_triu' + diag(lv);

mcci = inv(mcc);
mvvi = inv(mvv);

dicdt = -mcci*diag(rc)*ic + mcci*vc - mcv*ivdot - mcIp*ipdot;
divdt = -mvvi*diag(rv)*iv - mcv'*icdot - mvIp*ipdot;


ic = ic + dicdt*ts;
i = ic(circ.ii_unipolar) < 0;
ic(circ.ii_unipolar(i)) = 0;

iv = iv + divdt*ts;

x = [ic; iv];
y = x;








































