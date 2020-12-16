  function az= fil_a(x,y,z,c)

% PURPOSE:
% Generates the Z directed vector potential Az from by a straight filament
% centered at the origin [0,0,0] carrying Unit current (1A) in the +Z direction
% Filament limits are z=-c/2 to z=+c/2
% Only component is Az, All other components are zero (Ax, Ay == 0)
% Note: All Units are MKS (meters, Amps, Tesla)
%
% SYNTAX:
%        az= fil_a(x,y,z,c)
%
% INPUT:
%       x,y,z = coordinates of field point for calculation of Az [m]
%       c     = length of filament (scalar or same length as x) [m]
%
%       Note: x,y,z or c can be a scalar or matrix but all matricies must be
%             same size. Any of the components can be a scalar. 
%
% OUTPUT:
%       az == Vector Potential, Az, at point x,y,z [Web/m]
%
% Caution: On filament (x=y=0) solution is infinite between -a/2<=z<=a/2
%          solution gives INF. On axis (r=0), above and below the filament
%          the solution is finite and correctly calculated
%          DO NOT USE FOR SELF INDUCT.
%

% NOTE: see filr_a for filaments with finite radius computation inside rad.??
% Jim Leuer, General Atomics, 1-25-00
% Reference: see d3d book 26 & a_fil.mth & a_fil.xls
% 6/8/01 found that NAN result could actually be an INF result 
%        and oscillation in solution for r=>0 for +z but not for -z
%        modified to only compute az in -z plane and use symmetry to get +z
%        this allows removal of all NAN and INF stuff and is much more efficient
%        however not sure why negitive plane works but + plane has oscillations
%        as r=> sqrt(eps) to zero
% ==========================================================================

% Initialization stuff

  if nargin <= 3
    prog_name= 'fil_a';
    disp('%ERROR fil_a: Must have at least 4 arguments');
    help fil_a
    return
  end

%  warning off % this suppreses warning messages for NAN
  
% ---------------------------------------------------------------------------
% compute magnetics using the least memory and fastest computation

% ---------------
% start computation

% Constants
  x= x.^2+y.^2; % Caution reuse x = r^2 to reduce memory

% ---------------
% fix error seems to be to use up/down symmetry and use -z
% then we have no NAN or INF and no oscillation near r=0; wierd
  z= -abs(z); 
  z= z - 0.5*c;               % Caution reuse of z as z-c/2 lower mem
  az= sqrt(z.^2 + x)-z;       % numerator 
  z= z + c;                   % this is actualy z+c/2
  az= az./(sqrt(z.^2 + x)-z); % devide by numerator

% fix NAN for r=0, Z>a/2 (note negative Z plane z<-a/2 works without fix: wierd) 

% -abs(z) should fix all NAN's so dont needed
  
%  id=find(~isfinite(az)); % assumes all NAN & INF are correctable
%  if ~isempty(id)
%    if length(z)>1  z1= z(id); else z1= z; end % needed for z scalar
%    if length(c)>1  c1= c(id); else c1= c; end % needed for c scalar
%    az(id)= z1./(z1-c1);
%  end

%  warning backtrace % turn warning message back on

%  az= log(az);
  az= 1.0e-07*log(az);             % constant I*mu/(4pi)= 1*4pi*1e-7/(4pi)

% ---------------
  return

% ==========================================================================
% Testing Below All ok on Unix 1-26-01 JAL
% ---------------

% test input parameters

%  if exist('fil_a.diary') == 2 
%     ! rm fil_a.diary
%  end
%  format compact
%  diary fil_a.diary

%  nargin= 4;
  x= [0;1;-1;0;0];
  y= [0;1;-1;0;0];
  z= [0;1;-1;2;-2];
  c= 1;
  az= fil_a1(x,y,z,c)
  
% check out axis r=0  z < -a/2  
  x= 0;  y= 0; z= -2;  c= 1;

% check out axis r=0  z  +a/2  
  x= 0;  y= 0; z= +2;  c= 1;

% problem with "inf" result 6/8/01:
% Actually This shows large error in cylinder around z axis
  x=-1.243556591035144e-09
  y= 7.755451975066663e-10
  c= 0.10895123259799
  z= 0.21790246551551
  z= .8*c;
  
  x= linspace(0,10e+0,100)'*x;
  y= linspace(0,10e+0,100)'*y;
  z= ones(size(x))*z;
%  z= ones(size(x))*z*.1; % point 1 is inside c & gives imaginary solution 
  r= sqrt(x.^2+y.^2);  

  a= fil_a1(x,y,z,c);

  z=-z;
  am= fil_a(x,y,z,c);
  plot(r,a,'b',r,am,'r')  

% LOOKS LIKE -Z COMPUTS CORRECT ANSWER WHILE +Z OSCILLATES???
% FROM SYMMETRY CAN ALWAYS US -|Z| TO GET CORRECT ANSWER => NEEDS CHECKING INTO  
% LOOKS GOOD 6/8/01

% check out axis r=0  z < -a/2  
  x= [0;0];
  y= [0;0];
  z= [2;-2];
  c= 1;
  
%  x= x(1); y= y(1); z= z(1); c= c(1); 

   az=fil_a1(x,y,z,c)

% START CONOUR PLOT: Compare results of B computed using analytical and grad(Az)

   dx= 0.2; dy=0.2 % CAUTION: dont use smaller than dx=0.2 for quiver
   xo= 2; yo= 2;
   dx= xo/10; dy=yo/10;
%   xo=1e-9; yo=1e-10;
%   dx= xo/100; dy=yo/100;
   [x,y] = meshgrid(-xo:dx:xo, -yo:dy:yo);
   z= +0; 
   z=-1; 
%  z=-2; 
  z=100
   c= 1;

   figure(1)
   az=fil_a1(x,y,z,c);
   az= az*1e+7;
   [con,hdl]= contour(x,y,az);
   clabel(con,hdl);
   axis image
   grid on
   
   figure(2)
   clf
   [bx,by,bz]= fil_3d_b(x,y,z,c);
   b= 1e+7*sqrt(bx.^2+by.^2);
   [con,hdl]= contour(x,y,b);
   clabel(con,hdl);
   hold on
   quiver(x,y,bx,by);
   axis image
   grid on

   figure(3)
   clf
   [dadx,dady] = gradient(az,dx,dy);
    bx1= dady;
    by1= -dadx;
    b1= sqrt(bx1.^2+by1.^2);
   [con,hdl]= contour(x,y,b1);
   clabel(con,hdl);
   hold on
   quiver(x,y,bx,by);
   axis image
   grid on

   figure(4)
   id= 2:length(x)-1; % remove edges since not good accuracy at edge
   b1e= (b-b1)./b*1e+6;
   [con,hdl]= contour(x(id,id),y(id,id),b1e(id,id));
   clabel(con,hdl);
   axis image
   title('% error between Analytical B and B from grad(A) (ppm error)')
   grid on
