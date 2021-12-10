% Ref: https://iopscience.iop.org/article/10.1088/0741-3335/55/9/095009/meta
%
% EXAMPLE: let rbbbs, zbbbs be the boundary:
% 
% [rout, iout] = max(rbbbs);
% zout = zbbbs(iout);
% 
% [rin, iin] = min(rbbbs);
% zin = zbbbs(iin);
% 
% [ztop, itop] = max(zbbbs);
% rtop = rbbbs(itop);
% 
% [zbot, ibot] = min(zbbbs);
% rbot = rbbbs(ibot);
% 
% THEN: 
%
% xi_ou = measure_squareness(rout, rtop, zout, ztop, rbbbs, zbbbs)
% xi_iu = measure_squareness(rin, rtop, zin, ztop, rbbbs, zbbbs)
% xi_ol = measure_squareness(rout, rbot, zout, zbot, rbbbs, zbbbs)
% xi_il = measure_squareness(rin, rbot, zin, zbot, rbbbs, zbbbs)

function xi = measure_squareness(r1, r2, z1, z2, rbbbs, zbbbs, plotit)

  if ~exist('plotit', 'var'), plotit = 1; end

  A = r1 - r2;
  B = z2 - z1;
  rellipse = linspace(r2, r1);
  zellipse = z1 + sign(B)*sqrt(B^2 - ((B/A)*(rellipse - r2)).^2);

  [rc, zc] = intersections([r2 r1], [z1 z2], rellipse, zellipse);
  [rd, zd] = intersections([r2 r1], [z1 z2], rbbbs, zbbbs);

  LOD = norm([rd - r2, zd - z1]);
  LOC = norm([rc - r2, zc - z1]);
  LCE = norm([rc - r1, zc - z2]);

  xi = (LOD - LOC) / LCE;


  if plotit 
    figure
    hold on
    axis equal
    plot(rbbbs, zbbbs)
    plot(rellipse, zellipse)
    scatter([rc rd r2 r1], [zc zd z1 z2], 30, 'k', 'filled')
    plot([r2 r1], [z1 z2])
  end
end
