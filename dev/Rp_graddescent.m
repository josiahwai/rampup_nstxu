eps = 2e-10;
Rp = 3e-6 * ones(15,1);

for kk = 1:5
  build_nstxu_rz_system
  mpc_trajectory_rz
  
  ip_sim = x_all(end,:)';
  ip_true = traj.x(:,end);
  
%   dip_sim = gradient(ip_sim);
%   dip_true = gradient(ip_true);
  
  Rp(1:end-1) = Rp(1:end-1) + eps * (diff(ip_sim) - diff(ip_true));
  Rp(end) = Rp(end-1);
  Rp = max(Rp, 0);
end

cftool(t,r)