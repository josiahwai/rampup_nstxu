
function gaps = get_nstxu_gaps_from_bbbs(rbbbs, zbbbs, opts)

if ~exist('opts','var'), opts = struct; end
if ~isfield(opts, 'plotit'), opts.plotit = 0; end
if ~isfield(opts, 'verbose'), opts.verbose = 0; end
if ~isfield(opts, 'use_out_up_lo'), opts.use_out_up_lo = 0; end

% The cols of segs are := [r_start r_end z_start z_end]
% (r_start, z_start) should be the end point closest to the limiter, and
% (r_end, z_end) should be the point closer to core. 


segs = [0.2950    0.7778         0         0;
        0.2950    0.7778   -0.3813   -0.2402;
        0.2950    0.7587   -0.9342   -0.5071;
        0.8045    0.8045   -1.5520   -0.5122;
        1.2621    0.9616   -1.1173   -0.4595;
        1.4908    1.0867   -0.4385   -0.2593;
        1.5747    1.0371         0         0; 
        1.4908    1.0867    0.4385    0.2593;
        1.2621    0.9616    1.1173    0.4595;
        0.8045    0.8045    1.5520    0.5122;
        0.2950    0.7587    0.9342    0.5071;
        0.2950    0.7778    0.3813    0.2402];
      
if opts.use_out_up_lo
  segs = [0.2950    0.7778         0         0;
          1.5747    1.0371         0         0;
          1.2621    0.9616    1.1173    0.4595;
          1.2621    0.9616   -1.1173   -0.4595];        
end
      
for i = 1:size(segs,1)  
  [gaps.r(i,1), gaps.z(i,1)] = intersections(rbbbs, zbbbs, segs(i,1:2), segs(i,3:4));   
end


if opts.plotit

  r0 = segs(:,1);
  rf = segs(:,2);
  z0 = segs(:,3);
  zf = segs(:,4);
  
  plot(rbbbs, zbbbs)
  hold on
  plot([r0 rf]', [z0 zf]', 'color', [1 1 1] * 0.8, 'linewidth', 3)
  scatter(gaps.r, gaps.z, 60, 'b', 'filled')
end


















