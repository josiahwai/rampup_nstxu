  function [bx,by,bz]= cfil_3d_b(x,y,z,a)

% PURPOSE:
% Generates B-field in carthesian 3-space [X] == [x,y,z] 
% produced by a circular filament with radius a located with center at 0,0,0
% and center line in +Z direction. Unit current (1A) in counter clockwise
% The coordinate system is centered about the geometric center of filament.
% Note: All Units are MKS (meters, Amps, Tesla)
%
% SYNTAX:
%        [bx,by,bz]= cfil_3d_b(x,y,z,a)
%
% INPUT:
%       [x,y,z] == [X]  vectors of coordinates, one point per row
%       [a] ==          vector  of filament radius one filament per row
%
%  Note: If x,y,z & a are matrix or vectors they must be same size.
%        Any of the components can be a scalar. 
%
% OUTPUT:
%       [bx,by,bz] == [B] vectors of B-fields at point X
%
% NOTE: use cfilr_3b_b if [X] is inside the filament
%

% Jim Leuer, General Atomics, 2-98

% ---------------
% initialization stuff


  if nargin <= 3
    prog_name= 'cfil_3d_b';
    disp(['%ERROR ',prog_name,': Must have at least 4 arguments']);
    eval(['help ',prog_name,]) % print out help from prog_name
    return
  end

% ---------------------------------------------------------------------------
% compute magnetics using the least memory and fastest computation

% Constants

  r2=  x.^2 + y.^2;
  
  r=   sqrt(r2);
  z2=  z.^2;
  rma= (r-a).^2 + z2;
  rpa= (r+a).^2 + z2;
  k2=  4*a.*r./rpa;

% ---------------
% elliptic integral

  [k,e]= ellipke(k2);
  
% ---------------
% calculate fields


  r= sqrt(rpa);
  bx= -2e-7*z./(r2.*r).*(k-0.5*(rma+rpa)./rma.*e); %sign error 020403  
  by= y.*bx;
  bx= x.*bx; % sign error 020403
  bz= 2e-7./r.*(k - (r2 - a.^2 + z2)./rma.*e);  

% eliminate problem with points on centerline which are NAN

  id= find(isnan(bx));

  if ~isempty(id)
    disp('cfil_3d_b: correcting centerline NAN problem')
    bx(id)= zeros(id,1);  
    by(id)= bx(id);
    bz(id)= 2e-7*pi*a(id).^2./r(id).^3;
  end

  return

% ------------------------------------------------------------
% TEST Stuff Below:
% ------------------------------------------------------------
% sign problem found 4-3-02 testing below:
%rc= 0.8608
%zc= 0.1685
rc= 4
zc= .2

rpt= 3
zpt= 1
ppt= 30

xpt= rpt*cos(pi/180*ppt)
ypt= rpt*sin(pi/180*ppt)

% check field at center of dipole

x= xpt;
y= ypt;
z= zpt-zc;
a= rc;

 [bx,by,bz]= cfil_3d_b(x,y,z,a);

%
  del= .001;  
  ri= a*(1-del);
  ro= a*(1+del);
  zl= zc - del*a;
  zu= zc + del*a;
  cur= 4*pi*1e-7*ones(size(x));
  rp= sqrt(x.^2+y.^2);
  zp= zpt;
  [br1,bz1,fl1]= fine(ri,ro,zl,zu,cur,rp,zp); % note cur contains muo
  b1= sqrt(br1.^2 + bz1.^2);
  r= sqrt(x.^2+y.^2);
  bx1= x./r.*br1;
  by1= y./r.*br1;

  [bx,bx1]
  [by,by1] % this one is opposite sign??
  [bz,bz1]

% =======================================================  
% Next Test

% test input parameters

  if exist('fil_3d_b.diary') == 2 
     ! rm fil_3d_b.diary
  end
  format compact
  diary fil_3d_b.diary
  nargin= 4;

%  x= 0:.01:1; y= 0:.02:2; z= 0:.02:2; c=ones(size(x)); d=ones(size(x));
%  x= 0:.01:2; y= zeros(size(x)); z= 0:.01:2; c=ones(size(x)); d=ones(size(x));
%  x= 0:.01:2; y= zeros(size(x)); z= zeros(size(x));
%  x= 0; y= zeros(size(x)); z= zeros(size(x));
  zo= +2
  x= (0:.1:1.0001)'; y= zeros(size(x)); z= zo*ones(size(x));
  a=ones(size(x));
  nargin=4;
  [bx,by,bz]= cfil_3d_b(x,y,z,a);
%  cfil_3d_b;
  br= sqrt(bx.^2+by.^2);
  b= sqrt(bx.^2+by.^2+bz.^2);

%
  del= .001;  
  ri= a*(1-del);
  ro= a*(1+del);
  zl= -del*a;
  zu= +del*a;
  cur= 4*pi*1e-7*ones(size(x));
  rp= sqrt(x.^2+y.^2);
  zp= z;
  [br1,bz1,fl1]= fine(ri,ro,zl,zu,cur,rp,zp); % note cur contains muo
  b1= sqrt(br1.^2 + bz1.^2);

  r= sqrt(x.^2+y.^2);
  bx1= x./r.*br1;
  by1= y./r.*br1;
  id= find(isnan(bx1));
  if ~isempty(id)
    bx1(id)= zeros(length(id),1);
    by1(id)= zeros(length(id),1);
  end
  disp('[bx,bx1,delta]')
  disp([bx,bx1,bx-bx1])

  disp('[by,by1,delta]')
  disp([by,by1,by-by1])

  disp('[bz,bz1,delta]')
  disp([bz,bz1,bz-bz1])

  err=  (b-b1)./b

  figure(1)
  clf
  plot(r,(b-b1)./b)
  hold on
  grid on
%  axis([0,2,0,.4])


  format short
  disp(' Test Output of fil_3d_b:')
  disp('     x         y         z         c')
  disp([x,y,z,c]);

% ---------------
 print results below:
 
  b= sqrt(bx.^2+by.^2);
  format long
  disp('           bx                by                  bz              B')
  disp([bx,by,bz,b]);
  format short

  diary off
