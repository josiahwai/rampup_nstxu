
% test_fine.m
%  generic_startup

  rc= 0.5;
  zc= 0.1;
  dr= 0.01;
  dz= 0.01;
  r1= rc-.5*dr;
  r2= rc+.5*dr;
  z1= zc-.5*dz;
  z2= zc+.5*dz;

%  zp=  (-.099:-.0001:-.1)';
%  zp=  (1:-1:-100)';
%  zp=   -1.091642;
%  rp=  2.072899*ones(size(zp));

  rp= rc;
  zp= [-zc-1 -zc -zc+1];
  drp= dr;
  dzp= dz;

  ri= r1*ones(size(zp));
  ro= r2*ones(size(zp));
  zl= z1*ones(size(zp));
  zu= z2*ones(size(zp));
  cur= 1*ones(size(zp));
  rp =rp*ones(size(zp));

  [hr,hz,fl] =fine(ri,ro,zl,zu,cur,rp,zp);

  amu= 4*pi*1e-7;
  disp(['Fine: flux (H)= ', num2str(fl*amu)])


  [mut,br,bz]=mindbf(rc,zc,rp,zp);

  disp(['Green: flux (H)= ', num2str(mut)])

%  fl2= mut.*cur/amu;

%  disp([zp,fl,fl2,(fl2-fl)./fl2*100])

%  err= (fl2-fl)./fl2*100

%  plot(err)

