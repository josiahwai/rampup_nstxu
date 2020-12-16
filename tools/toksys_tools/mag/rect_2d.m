  function [bx,by]= rect_2d(a,b,cur,x,y)

% PURPOSE:
% Generates B-field in 2-d planar space [X] == [x,y] at point x,y
% produced by a uniform current density straight rectangular coil 
% of inifinte length in Z direction having dimensions a and b in the x and y
% directions and centered at x=y=0 & carrying current cur.
% Note: All Units are MKS (meters, Amps, Tesla)
%
% SYNTAX:
%        [bx,by]= rect_2d(a,b,cur,x,y)
%
%       a,b    = x and y dimensions of infinite length rectangular block
%       cur    = current in rectangle
%       x,y    = location of field point to determine field
%
%  Note: Should work for vectors all of the same size
%
% OUTPUT:
%       [bx,by] == [B] vectors of B-fields at point X
%
% NOTE: Z is the ignorable dimension 
%

% Jim Leuer, General Atomics, 3-98

% ---------------
% initialization stuff


  if nargin <= 4
    prog_name= 'rect_2d';
    disp(['%ERROR ',prog_name,': Must have at least 5 arguments']);
    eval(['help ',prog_name,]) % print out help from prog_name
    return
  end

% ---------------------------------------------------------------------------
% compute magnetics using the least memory and fastest computation

% Constants (use symmetry)
%

  xm= abs(x) - a/2;     % abs() => calculate in 1st quad and use symmetry
  xp= xm + a;
  ym= abs(y) - b/2;
  yp= ym + b;
  xpyp= xp.^2 + yp.^2;
  xpym= xp.^2 + ym.^2;
  xmyp= xm.^2 + yp.^2;
  xmym= xm.^2 + ym.^2;
  
  cur= cur./(a.*b);   % cur is now J
  cur= 2.0e-7*cur;    % mu*J/(2*pi)

% ---------------------------------------------------------------------------
% formula

  bx=   yp.*( atan(xp./yp) - atan(xm./yp) )...
      - ym.*( atan(xp./ym) - atan(xm./ym) )...
      + 0.5*( xp.*log(xpyp./xpym) - xm.*log(xmyp./xmym) );

  by=   xp.*( atan(yp./xp) - atan(ym./xp) )...
      - xm.*( atan(yp./xm) - atan(ym./xm) )...
      + 0.5*( yp.*log(xpyp./xmyp) - ym.*log(xpym./xmym) );

  id= find(isnan(bx));
  if ~isempty(id)
    disp('rect_2d: correcting centerline NAN problem')
%           using symmetry xm & ym are zero
    bx(id)=   yp(id).*( atan(xp(id)./yp(id)) - atan(xm(id)./yp(id)) )...
            - ym(id).*( atan(xp(id)./ym(id)) )...
            + 0.5*( xp(id).*log(xpyp(id)./xpym(id)) ); 

    by(id)=   xp(id).*( atan(yp(id)./xp(id)) - atan(ym(id)./xp(id)) )...
            - xm(id).*( atan(yp(id)./xm(id)) )...
            + 0.5*( yp(id).*log(xpyp(id)./xmyp(id)) );

  end

% cur=1; % test to make numbers around 1

  bx= -cur*sign(y).*bx; % bring back all 4 quadrans
  by=  cur*sign(x).*by;

  return

% ------------------------------------------------------------
% TEST Stuff Below:

% test input parameters

  format long
  nargin= 5;

  xp= (0:.1000:1.0001)'; yp= zeros(size(xp));
  r= sqrt(xp.^2 + yp.^2);
  a=ones(size(xp));
  b=ones(size(xp));
  cur= ones(size(xp));

  [bx1,by1]= rect_2d(a,b,cur,xp,yp);
  br1= sqrt(bx1.^2+by1.^2);

%  diary rect_2d.diary

  disp('[bx1,bx1,br1]')
  disp([bx1,by1,br1])

%  diary off


  figure(1)
  clf
  plot(r,br1)
  hold on
  grid on
%  axis([0,2,0,.4])

  
% ---------------
% Contour Test
  format long
  nargin= 5;

  xp= (-1.:.01:1.000);
  yp= (-1.:.01:1)';
%  yp= (-2.:.01:2)';
  nx= length(xp);  ny= length(yp);

  xpp= ones(length(yp),1)*xp;
  ypp= yp*ones(1,length(xp));
  xpp= reshape(xpp,nx*ny,1);
  ypp= reshape(ypp,nx*ny,1);
  a= 1; b=1; cur= 1;
  
%  cur= ones(nx*ny,1);

  [bx,by]= rect_2d(a,b,cur,xpp,ypp);

  bx= reshape(bx,ny,nx);
  by= reshape(by,ny,nx);

  bt= sqrt(bx.^2+by.^2);

  xpp= reshape(xpp,ny,nx);
  ypp= reshape(ypp,ny,nx);

% calculate field from filament:

  r= sqrt(xpp.^2 + ypp.^2);
  bt1= 2.0e-7./r;
  dbt= bt-bt1;

  figure(1)
  clf
  contour(xpp,ypp,bt)
  axis equal
  hold on
  grid on
%  axis([0,2,0,.4])

% ---------------
% look at individual components around square
  format short
  i=50; j=50;
  disp([i,j,xpp(j,i),ypp(j,i),bx(j,i),by(j,i)])
  i=50; j=101;
  i=50; j=152;
  i=101; j=50;
  i=101; j=101;
  i=101; j=152;
  i=152; j=50;
  i=152; j=101;
  i=152; j=152;
