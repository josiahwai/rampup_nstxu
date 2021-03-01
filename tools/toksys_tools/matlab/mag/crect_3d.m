  function [bx,by,bz]= crect_3d(rc,zc,dr,dz,cur,xp,yp,zp)

% PURPOSE:
% Generates B-field in carthesian 3-space [X] == [x,y,z] at point xp,yp,zp 
% produced by a uniform current density coil located at rc,zc with axis in +Z
% and having a rectangular cross section dr x dz and carrying current cur.
% Note: All Units are MKS (meters, Amps, Tesla)
%
% SYNTAX:
%        [bx,by,bz]= crect_3d(rc,zc,dr,dz,cur,xp,yp,zp)
%
%       rc,zc   = center of rectangle
%       dr,dz   = radial width and height of rectangle
%       cur     = current in rectangle
%       xp,yp,zp= location of field point to determine field
%
%  Note: Works for vectors all of the same size
%
% OUTPUT:
%       [bx,by,bz] == [B] vectors of B-fields at point X
%
% NOTE: use fine so see limintations in fine.m
%

% Jim Leuer, General Atomics, 2-98

% ---------------
% initialization stuff


  if nargin <= 7
    prog_name= 'crect_3d';
    disp(['%ERROR ',prog_name,': Must have at least 8 arguments']);
    eval(['help ',prog_name,]) % print out help from prog_name
    return
  end

% ---------------------------------------------------------------------------
% compute magnetics using the least memory and fastest computation

% Constants
%
  ri= rc - 0.5*dr;
  ro= ri + dr;
  zl= zc - 0.5*dz;
  zu= zl + dz;  
  rp= sqrt(xp.^2+yp.^2);
  cur= 4e-7*pi*cur; % fine returns H change from H to B)
  [br,bz,fl]= fine(ri,ro,zl,zu,cur,rp,zp); % note cur contains muo
  
% new fix for centerline problem

  id= find(rp > eps);
  bx= zeros(size(br));
  by= bx;
  bx(id)= xp(id)./rp(id).*br(id);
  by(id)= yp(id)./rp(id).*br(id);

%  id= find(isnan(bx));
%  if ~isempty(id)
%    disp('crect_3d: correcting centerline NAN problem')
%    bx(id)= zeros(length(id),1);
%    by(id)= zeros(length(id),1);
%  end

  return

% ------------------------------------------------------------
% TEST Stuff Below:

% test input parameters

%  format long
%  nargin= 8;

%  x= 0:.01:1; y= 0:.02:2; z= 0:.02:2; c=ones(size(x)); d=ones(size(x));
%  x= 0:.01:2; y= zeros(size(x)); z= 0:.01:2; c=ones(size(x)); d=ones(size(x));
%  x= 0:.01:2; y= zeros(size(x)); z= zeros(size(x));
%  x= 0; y= zeros(size(x)); z= zeros(size(x));
%  zo= +1
%  x= (0:.1:1.0001)'; y= zeros(size(x)); z= zo*ones(size(x));
%  a=ones(size(x));
%  nargin=4;
%  [bx,by,bz]= cfil_3d_b(x,y,z,a);
%  cfil_3d_b;
%  br= sqrt(bx.^2+by.^2);
%  b= sqrt(bx.^2+by.^2+bz.^2);

%
%  rc= a;
%  zc= zeros(size(rc));
%  del= .001;  
%  dr= del*ones(size(rc));
%  dz= del*ones(size(rc));
%  cur= ones(size(rc));

%  cur= 4*pi*1e-7*ones(size(x));
%  xp= x; yp= y; zp= z;
%  [bx1,by1,bz1]= crect_3d(rc,zc,dr,dz,cur,xp,yp,zp);
%  br1= sqrt(bx.^2+by.^2);
%  b1= sqrt(bx.^2+by.^2+bz.^2);

% disp('[bx,bx1,delta]')
%  disp([bx,bx1,bx-bx1])

%  disp('[by,by1,delta]')
%  disp([by,by1,by-by1])

%  disp('[bz,bz1,delta]')
%  disp([bz,bz1,bz-bz1])

%  err=  (b-b1)./b

%  figure(1)
%  clf
%  plot(r,(b-b1)./b)
%  hold on
%  grid on
%  axis([0,2,0,.4])


%  format short
%  disp(' Test Output of fil_3d_b:')
%  disp('     x         y         z         c')
%  disp([x,y,z,c]);

% ---------------
% print results below:
 
%  b= sqrt(bx.^2+by.^2);
%  format long
%  disp('           bx                by                  bz              B')
%  disp([bx,by,bz,b]);
%  format short

%  diary off
