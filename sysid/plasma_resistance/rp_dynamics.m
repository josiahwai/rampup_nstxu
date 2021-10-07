function [x,y] = rp_dynamics(t, x, u, Rp_t, Lp_t, file_args)


% Unpack file_args
traj = file_args{1};
parameter_times = file_args{2};
Ts = file_args{3};

Mpc = squeeze( interp1(traj.t, traj.mcIp, t, 'linear', 'extrap'));
Mpv = squeeze( interp1(traj.t, traj.mvIp, t, 'linear', 'extrap'));
Rp = interp1(parameter_times, Rp_t, t, 'linear', 'extrap');
Lp = interp1(parameter_times, Lp_t, t, 'linear', 'extrap');

Ac = -Rp / Lp;
Bc = -1/Lp * [Mpc Mpv];

[Ad,Bd] = c2d(Ac,Bc,Ts);

x = Ad*x(:) + Bd*u(:);
y = x;







































