  function m= mut_fil_fil(xa,xb,fa,fb)

% PURPOSE:
% Calculate the mutual inductance between two sets of filaments A & B.
% Filament is defined by center and deltas x= [xc,yc,zc,dx,dy,dz]. 
% Sending filament set is "A", receiving filament set is "B"
%
% Note: All Units are MKS (meters, Amps, Tesla, Henries ...)
%
% SYNTAX:
%        mut_fil_fil(xa,xb,fa,fb)
%
% INPUT: (a= sending filament, b= receiving filament)
%       xa   = center coordinates. & delta lengths [xa,ya,za,dxa,dya,dza]
%       xb   = same for receiving filaments [xb,yb,zb,dxb,dyb,dzb]
%       fa,fb= Optional: fraction of current in filament for parallel filaments
%
% OUTPUT:
%       m == mutual inductance between filament sets A to B [H]
%
% CAUTION: Will not work for filaments on top of each other (ie self inductance)

% Jim Leuer General Atomics 11-29-01

% =============================================================================

  if nargin <=1
    disp('%ERROR mut_fil_fil: Must have at least 2 arguments');
    help mut_fil_fil
    return
  elseif nargin <= 3
    fa= ones(length(xa(:,1)),1);
    fb= 1;
  end

  
% ----------------------------------------------------------------------
% START Calculation Vector Potential from filament A to filament b
% ----------------------------------------------------------------------
% M_ba = sum_a[sum_b(f_a*A_ba*dl_b*f_b)];

% ------------------------------------
  na= length(xa(:,1));
  la= sqrt(xa(:,4).^2+xa(:,5).^2+xa(:,6).^2);% length of sending filament
  m= 0;
  for ii= 1:na % Loop over all sending filaments A
% note: z_dircos takes direction number or direction cosine as argument
    fil_rot= z_dircos(xa(ii,4:6)); % rotational transform X_local= rot*X_glob
    fil_rot_t= fil_rot';
% ------------------------------------
%  field point infomation
    xp= [xb(:,1)-xa(ii,1)...
         xb(:,2)-xa(ii,2)...
	 xb(:,3)-xa(ii,3)]; % translation global coord to filament A center
    xp= (fil_rot*xp')';     % rotate global to local filament A coordinates
    a_z= fil_a(xp(:,1),xp(:,2),xp(:,3),la(ii)); % vector potential A=>B
    a_xyz= (fil_rot_t(1:3,3)*a_z')'; % local Az to Global Ax,Ay,Az
    a_dl= sum((a_xyz.*xb(:,4:6))')'; % dot product A_ba[x,y,z].[dxb,dyb,dzb]
    m= m+fa(ii)*sum(a_dl.*fb); % sum over all receiving filaments     
  end % end sending filaments A loop

  return

%?????????????????????????????????????????????????????????
% TEST1 coaxial elements sqrt(3) long seperated by sqrt(3)
%  nargin=2;
  fa=1;
  fb=1;
  
  figure(1)
  clf

  na0= 10;
  xpta= linspace(0,1,na0+1)';
  
  xpta= [xpta, xpta+2, xpta+3];  
  plot3(xpta(:,1),xpta(:,2),xpta(:,3),'r')
  hold on
  xfila= pts_to_fil(xpta);
%  arrow(xfila(:,1:3),xfila(:,4:6));
  plot3(xpta(:,1),xpta(:,2),xpta(:,3),'ro')
  xa= fil_to_filc(xfila);
  la= sqrt(xa(:,1).^2 + xa(:,2).^2 + xa(:,2).^2);  
  ia= ones(size(la));
   grid on
   
   nb0= 10;
   xptb= linspace(xpta(end,1)+1,xpta(end,1)+2,nb0+1)';
   xptb= [xptb, xptb+2, xptb+3];  

   plot3(xptb(:,1),xptb(:,2),xptb(:,3),'b');
   xfilb= pts_to_fil(xptb);
%   arrow(xfilb(:,1:3),xfilb(:,4:6));
   xb= fil_to_filc(xfilb);  
   lb= sqrt(xb(:,1).^2 + xb(:,2).^2 + xb(:,2).^2);  
   ib= ones(size(lb));
   plot3(xptb(:,1),xptb(:,2),xptb(:,3),'bx');


format long

  mba=mut_fil_fil(xa,xb)
  mab=mut_fil_fil(xb,xa)
  
  
%ref grover Pg46 Eqn 35

  mg= sqrt(3); lg=sqrt(3); dg= sqrt(3);
  mgrover= 1e-7*((lg+mg+dg)*log(lg+mg+dg)...
                   -(lg+dg)*log(lg+dg)...
		   -(mg+dg)*log(mg+dg)...
		        +dg*log(dg) )

 error= (mba-mgrover)/mgrover


%?????????????????????????????????????????????????????????
% TEST 2 parallel conductors - see grover pg 45 eqn (28)
% sending A at 0,0,0 to 1,0,0
% receiving B at 2,0,1 to 3,0,1

%  figure(2)
%  clf

%  nargin=2;
  fa=1;
  fb=1;

% grover distance parameters for parallel filaments 
  la_g= 1;  % length of A filament
  lb_g= 1;  % length of B filament
  dab_g= 1; % distance between ends of filaments
  gap_g= 1; % perpendicular lenght between filament axis

  na0= 10;
  dcosa= [1,0,0];
  xpta0= [0,0,0]; % origin of filament a
  xpta1= xpta0+dcosa*la_g;
  xpta= [linspace(xpta0(1),xpta1(1),na0+1)'...
         linspace(xpta0(2),xpta1(2),na0+1)'...
         linspace(xpta0(3),xpta1(3),na0+1)'];

  plot3(xpta(:,1),xpta(:,2),xpta(:,3),'r')
  hold on
  xlabel('X'); ylabel('Y'); zlabel('Z')
  xfila= pts_to_fil(xpta);
  plot3(xpta(:,1),xpta(:,2),xpta(:,3),'ro')
  xa= fil_to_filc(xfila);
  la= sqrt(xa(:,1).^2 + xa(:,2).^2 + xa(:,2).^2);  
  ia= ones(size(la));
   grid on

  nb0= 10;
  dcosb= dcosa;
  dcos_perp= [0. 0. 1];
  xptb0= xpta(end,:) + dcosb*dab_g + dcos_perp*gap_g; % origin of filament b
  xptb1= xptb0+dcosb*lb_g;
  xptb= [linspace(xptb0(1),xptb1(1),nb0+1)'...
         linspace(xptb0(2),xptb1(2),nb0+1)'...
         linspace(xptb0(3),xptb1(3),nb0+1)'];
   plot3(xptb(:,1),xptb(:,2),xptb(:,3),'b');
   xfilb= pts_to_fil(xptb);
   xb= fil_to_filc(xfilb);  
   lb= sqrt(xb(:,1).^2 + xb(:,2).^2 + xb(:,2).^2);  
   ib= ones(size(lb));
   plot3(xptb(:,1),xptb(:,2),xptb(:,3),'bx');


format long
%ref grover Pg45 Eqn 28
   a= la_g + lb_g + dab_g;
   b= la_g + dab_g;
   g= lb_g + dab_g;

  mba=mut_fil_fil(xa,xb)
  mab=mut_fil_fil(xb,xa)
    
  mgrover= 1e-7*(a*asinh(a/gap_g)      -b*asinh(b/gap_g)...
                -g*asinh(g/gap_g)  +dab_g*asinh(dab_g/gap_g)...
		-sqrt(a^2+gap_g^2) +sqrt(b^2+gap_g^2)...
		+sqrt(g^2+gap_g^2) -sqrt(dab_g^2+gap_g^2))
 
 error= (mba-mgrover)/mgrover

% Now check snow National Bureau of Standards circular 544 1954 pg 28 # 2.19
 l1= la_g;
 l2= lb_g;
 c=  dab_g +0.5*(l1+l2);
 d=  gap_g;
 x1= c+0.5*(l2+l1);
 x2= c-0.5*(l2+l1);
 x3= c+0.5*(l2-l1);
 x4= c-0.5*(l2-l1);
 msnow= 1e-7*( (abs(x1)*log((sqrt(x1^2+d^2)+abs(x1))/d)-sqrt(x1^2+d^2))...            
	      +(abs(x2)*log((sqrt(x2^2+d^2)+abs(x2))/d)-sqrt(x2^2+d^2))...
              -(abs(x3)*log((sqrt(x3^2+d^2)+abs(x3))/d)-sqrt(x3^2+d^2))...
              -(abs(x4)*log((sqrt(x4^2+d^2)+abs(x4))/d)-sqrt(x4^2+d^2)) )

%?????????????????????????????????????????????????????????
% TEST 3 two rectangular coils per snow NBS Pg 29 eqn 2.20

%  figure(2)
%  clf

%  nargin=2;
  fa=1;
  fb=1;

% grover distance parameters for parallel filaments 
  a= 1;  % length of one side of rectangles
  b= 1; % length of 2nd side of rectangles
  d= 1; % distance between rectangles

  na0= 10;
  xpta0= [0,0,0]
  dcosx= [1,0,0];
  dcosy= [0,1,0];
  dcosz= [0,0,1];
  xpta=  [0,0,0;...
          a,0,0;...
          a,b,0;...
          0,b,0;...
          0,0,0];
  xptb=  [0,0,d;...
          a,0,d;...
          a,b,d;...
          0,b,d;...
          0,0,d];   
  xpta= xpta + ones(length(xpta(:,1)),1)*[0,1,0]; % move in y by 1
  xptb= xptb + ones(length(xptb(:,1)),1)*[0,1,0];
  

  plot3(xpta(:,1),xpta(:,2),xpta(:,3),'r')
  hold on
  xlabel('X'); ylabel('Y'); zlabel('Z')
  xfila= pts_to_fil(xpta);
  na0= 10;
  xfila= fil_regrid(xfila,na0);
  
  plot3([xfila(:,1) xfila(:,4)]',...
        [xfila(:,2) xfila(:,5)]',...
	[xfila(:,3) xfila(:,3)]','ro');
%  plot3(xpta(:,1),xpta(:,2),xpta(:,3),'ro')
  xa= fil_to_filc(xfila);
  la= sqrt(xa(:,1).^2 + xa(:,2).^2 + xa(:,2).^2);  
  ia= ones(size(la));
   grid on

   plot3(xptb(:,1),xptb(:,2),xptb(:,3),'b');
   xfilb= pts_to_fil(xptb);
   nb0= 10;
   xfilb= fil_regrid(xfilb,nb0);
   plot3([xfilb(:,1) xfilb(:,4)]',...
         [xfilb(:,2) xfilb(:,5)]',...
	 [xfilb(:,3) xfilb(:,3)]','bx');
   xb= fil_to_filc(xfilb);  
   lb= sqrt(xb(:,1).^2 + xb(:,2).^2 + xb(:,2).^2);  
   ib= ones(size(lb));
%   plot3(xptb(:,1),xptb(:,2),xptb(:,3),'bx');


format long
%ref snow nbs 544 pg29 eqn 2.20

  mba=mut_fil_fil(xa,xb)
  mab=mut_fil_fil(xb,xa)
    
  msnow= 4e-7*(a*log((a+sqrt(a^2+d^2))/(a+sqrt(a^2+b^2+d^2))*sqrt(b^2+d^2)/d)...
              +b*log((b+sqrt(b^2+d^2))/(b+sqrt(a^2+b^2+d^2))*sqrt(a^2+d^2)/d)...
	      +2*(sqrt(a^2+b^2+d^2)-sqrt(a^2+d^2)-sqrt(b^2+d^2)+d))
  
 error= (mba-msnow)/msnow
