function [x,y] = nl_grey_nstxu_model(t, x, u, rcc, rvv, rp_t, lp_t, voltage_scale, file_args)

Ts = file_args{1};
circ = file_args{2};
lstar_invs = file_args{3};
sim_timebase = file_args{4};

[~,k] = min(abs(sim_timebase-t));
rp = rp_t(k);
rxx = [rcc; rvv; rp];

lstari = squeeze(lstar_invs(k,:,:));
% lp = lp_t(k);
% lstar(end,end) = lp;
% lstari = inv(lstar);

Ac = -lstari * diag(rxx);
Bc = lstari(:,1:circ.ncx);

Ac(circ.iicx_remove,:) = 0;
Ac(:,circ.iicx_remove) = 0;
Bc(circ.iicx_remove,:) = 0;
Bc(:,circ.iicx_remove) = 0;

Astab = removeFirstEigval(Ac);

[A, B] = c2d(Astab, Bc, Ts);
C = circ.Pxx_keep;

u = u(:) .* voltage_scale;
x = A*x + B*u;
y = (C*x)';



















