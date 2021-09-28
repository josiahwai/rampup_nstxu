
function gaps = get_nstxu_gaps(eq, opts)

if ~exist('opts','var'), opts = struct; end
if ~isfield(opts, 'plotit'), opts.plotit = 1; end
  
struct_to_ws(eq);

% The cols of segs are := [r_start r_end z_start z_end]
% (r_start, z_start) should be the end point closest to the limiter, and
% (r_end, z_end) should be the point closer to core. 
% segs = [1.5747    1.0371         0         0; 
%         0.2950    0.7778         0         0;
%         0.2950    0.7778   -0.3813   -0.2402;
%         0.2950    0.7587   -0.9342   -0.5071;
%         0.8045    0.8045   -1.5520   -0.6122;
%         1.2621    0.9616   -1.1173   -0.4595;
%         1.4908    1.0867   -0.4385   -0.2593;
%         0.2950    0.7778    0.3813    0.2402;
%         0.2950    0.7587    0.9342    0.5071;
%         0.8045    0.8045    1.5520    0.6122;
%         1.2621    0.9616    1.1173    0.4595;
%         1.4908    1.0867    0.4385    0.2593];

segs = [0.2950    0.7778         0         0;
        0.2950    0.7778   -0.3813   -0.2402;
        0.2950    0.7587   -0.9342   -0.5071;
        0.8045    0.8045   -1.5520   -0.6122;
        1.2621    0.9616   -1.1173   -0.4595;
        1.4908    1.0867   -0.4385   -0.2593;
        1.5747    1.0371         0         0; 
        1.4908    1.0867    0.4385    0.2593;
        1.2621    0.9616    1.1173    0.4595;
        0.8045    0.8045    1.5520    0.6122;
        0.2950    0.7587    0.9342    0.5071;
        0.2950    0.7778    0.3813    0.2402];

for iseg = 1:size(segs,1)
  
  seg = segs(iseg,:);
  
  [r0,rf,z0,zf] = unpack(seg);

  N = 20;
  rq = linspace(r0, rf, N);
  zq = linspace(z0, zf, N);
  psiq = bicubicHermite(rg, zg, psizr, rq, zq);

  % index points closest to intersect locations
  idx = find(diff(sign(psiq - psibry)));
    
  if isempty(idx)
    warning(['No intersect location found on segment ' num2str(iseg)])
    [r, z, dist, psi_r, psi_z] = unpack(nan(5,1));
    
  else
    
    if numel(idx) > 1
      warning('Found multiple intersect locations. Using intersect point closest to core.')
      idx = idx(end);
    end
    
    % Newton's method to find psi=psibry
    r = rq(idx);
    z = zq(idx);
    seghat = [rf-r0 zf-z0];
    seghat = seghat / norm(seghat);

    for dum = 1:10
      [psi, psi_r, psi_z] = bicubicHermite(rg, zg, psizr, r, z);
      directional_gradient = ([psi_r psi_z] * seghat');  
      ds = (psibry - psi) / directional_gradient;
      r = r + seghat(1) * ds;
      z = z + seghat(2) * ds;    
    end    
    dist = sqrt( (r-r0)^2 + (z-z0)^2);
  end
  
  gaps.r(iseg,:) = r;
  gaps.z(iseg,:) = z;
  gaps.dist(iseg,:) = dist;  
  gaps.psi_r(iseg,:) = psi_r;
  gaps.psi_z(iseg,:) = psi_z;
  gaps.segs(iseg,:) = seg;
    
end


if opts.plotit

  r0 = segs(:,1);
  rf = segs(:,2);
  z0 = segs(:,3);
  zf = segs(:,4);
  
  plot_eq(eq)
  plot([r0 rf]', [z0 zf]', 'color', [1 1 1] * 0.8, 'linewidth', 3)
  scatter(gaps.r, gaps.z, 60, 'b', 'filled')
end


















