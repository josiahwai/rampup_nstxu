  function l= self_fil(xa,aa,fa,full)

% PURPOSE:
% Calculate the self inductance of a set of filaments A
% Filament is defined by center and deltas x= [xc,yc,zc,dx,dy,dz]. 
% The filament radius, aa is needed to compute self inductance
% This is complementary to mut_fil_fil which computes mutual between
% two filament sets A & B.
%
% Note: All Units are MKS (meters, Amps, Tesla, Henries ...)
%
% SYNTAX:
%        self_fil(xa,fa)
%
% INPUT: (a= sending filament & receiving filament) [default]
%      xa  = center coordinates. & delta lengths [xa,ya,za,dxa,dya,dza]
%      aa  = Radius of filament for self inductance (vector or scalor)
%      fa  = Optional: fraction of current in filament for parallel currents [1]
%      full= 0=does only half matrix (assume symmetry) =1 does full matrix [0]
%
% OUTPUT:
%       l == self inductance of filament A sets [H]
%
% CAUTION: Filament radius aa cannot be zero

% Jim Leuer General Atomics 3-21-01

% =============================================================================

  if nargin <=0
    disp('%ERROR self_fil: Must have at least 1 argument');
    help self_fil
    return
  elseif nargin <= 1
    aamn= 0.01*min(sqrt(xa(:,4).^2 + xa(:,5).^2 + xa(:,6).^2));
    aa=   aamn*ones(length(xa(:,1)),1);	 
    disp(['%CAUTION self_fil setting filament radius to: ',mum2str(aa(1))])
  elseif nargin <= 2
    fa= ones(length(xa(:,1)),1);
    full= 0;
  elseif nargin <= 3
    full= 0;
  end

  if length(fa) ~= length(xa(:,1))
    disp('%Caution: self_fil: fa size changed to size of xa')
    fa= fa(1)*ones(length(xa(:,1)),1);
  end
  
% ----------------------------------------------------------------------
% START Calculation Vector Potential from filament A to filament A
% ----------------------------------------------------------------------
% L = M_aa = sum_a[sum_a(f_a*A_aa*dl_a*f_a)];

% ------------------------------------
  na= length(xa(:,1));
  la= sqrt(xa(:,4).^2+xa(:,5).^2+xa(:,6).^2);% length of sending filament
  l= 0;
% ===========================================================
  if na > 1 % Only do off diagionals if more than 1 element

% ===========================================================
   if full==1 % do all elements of matrix (above and below diagonals)

%   Start Full Matrix Computation of off diagonals
    iao= 1:na; % full index of filaments
    for ii= 1:na % sending filaments
     ia= find(ii~=iao); % remove diagonal element
     fil_rot= z_dircos(xa(ii,4:6)); % rotational transform X_local= rot*X_glob
     fil_rot_t= fil_rot';
     xp= [xa(ia,1)-xa(ii,1)...
          xa(ia,2)-xa(ii,2)...
          xa(ia,3)-xa(ii,3)]; % translation global coord to filament A center
     xp= (fil_rot*xp')';     % rotate global to local filament A coordinates
     a_z= fil_a(xp(:,1),xp(:,2),xp(:,3),la(ii)); % vector potential A=>B
     a_xyz= (fil_rot_t(1:3,3)*a_z')'; % local Az to Global Ax,Ay,Az
     a_dl= sum((a_xyz.*xa(ia,4:6))')'; % dot product A_ba[x,y,z].[dxb,dyb,dzb]
     l= l+fa(ii)*sum(a_dl.*fa(ia)); % sum over all receiving filaments     
    end % for ii=1:na

% ===========================================================
   else % do only elements above diagonals

%   Start Upper Diagonal and assume symmetry to get lower diagonal contribution
    for ii= 1:na-1 % sending filaments
     ia= (ii+1):na; % Upper region without diagonal
     fil_rot= z_dircos(xa(ii,4:6)); % rotational transform X_local= rot*X_glob
     fil_rot_t= fil_rot';
     xp= [xa(ia,1)-xa(ii,1)...
          xa(ia,2)-xa(ii,2)...
          xa(ia,3)-xa(ii,3)]; % translation global coord to filament A center
     xp= (fil_rot*xp')';     % rotate global to local filament A coordinates
     a_z= fil_a(xp(:,1),xp(:,2),xp(:,3),la(ii)); % vector potential A=>B
     a_xyz= (fil_rot_t(1:3,3)*a_z')'; % local Az to Global Ax,Ay,Az
     a_dl= sum((a_xyz.*xa(ia,4:6))')'; % dot product A_ba[x,y,z].[dxb,dyb,dzb]
     l= l+fa(ii)*sum(a_dl.*fa(ia)); % sum over all receiving filaments     
    end % for ii=1:na-1
    l= 2*l; % double from symmetry
   end % if full==1
  
  end % if na>1
% ===========================================================
%  Do diagonal computation
%  Using standard grover, AIP etc
   ll= 2e-7*la.*fa.^2.*(log(2*la./aa)-0.75); %CAUTION assumes aa<<la
   l= l+sum(ll);
%  Using Snow
%   ep= aa./la;
%   ll= 2e-7*la.*fa.^2.*(log((sqrt(1+ep.^2)+ep)./ep)-sqrt(1+ep.^2)+ep+0.25);
%   l= l+sum(ll);
  return

%?????????????????????????????????????????????????????????
% TEST1 coaxial elements sqrt(3) long seperated by sqrt(3)
  nargin=2;
  fa=1;
  aa=.0001;
  lg= sqrt(3); % length of line

  format long
  
  figure(1)
  clf

  na0= 100;
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
   
  na0
  laa= self_fil(xa,aa,fa,1)
  laa= self_fil(xa,aa,fa)

% ref grover pg 35 eqn 7 , Snow(nbs544) eqn 2-17  
  lgrover= 2e-7*lg*(log(2*lg/aa)-.75)
  egrov= (laa-lgrover)/lgrover

  lsnow=   2e-7*(lg*log((sqrt(lg^2+aa^2)+aa)/aa) - sqrt(lg^2+aa^2) + lg/4 + aa)
  esnow= (laa-lsnow)/lsnow


%?????????????????????????????????????????????????????????
% TEST 2 parallel conductors - see grover pg 45 eqn (28)
% sending A at 0,0,0 to 1,0,0
% receiving B at 2,0,1 to 3,0,1

%  figure(2)
%  clf

  nargin=2;
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
% TEST 3 rectangle coil per snow NBS Pg 29 eqn 2.21

%  figure(2)
%  clf

  nargin=2;
  fa=1;
  aa=.0001;

% grover distance parameters for parallel filaments 
  a= 1;  % length of one side of rectangles
  b= 1; % length of 2nd side of rectangles
  d= 1; % distance between rectangles

  na0= 100;
  xpta0= [0,0,0]
  dcosx= [1,0,0];
  dcosy= [0,1,0];
  dcosz= [0,0,1];
  xpta=  [0,0,0;...
          a,0,0;...
          a,b,0;...
          0,b,0;...
          0,0,0];
  xpta= xpta + ones(length(xpta(:,1)),1)*[0,1,0]; % move in y by 1
  

  plot3(xpta(:,1),xpta(:,2),xpta(:,3),'r')
  hold on
  xlabel('X'); ylabel('Y'); zlabel('Z')
  xfila= pts_to_fil(xpta);
  xfila= fil_regrid(xfila,na0);
  
  plot3([xfila(:,1) xfila(:,4)]',...
        [xfila(:,2) xfila(:,5)]',...
	[xfila(:,3) xfila(:,3)]','ro');
%  plot3(xpta(:,1),xpta(:,2),xpta(:,3),'ro')
  xa= fil_to_filc(xfila);
  la= sqrt(xa(:,1).^2 + xa(:,2).^2 + xa(:,2).^2);  
  ia= ones(size(la));
   grid on


format long
%ref snow nbs 544 pg29 eqn 2.21

  laa=self_fil(xa,aa,fa)
    
  lsnow= 4e-7*(...
               (a+b)*log((sqrt(4*(a+b)^2+aa^2)+aa)/aa)...
	       - b  *log((sqrt(b^2+a^2)+b)/a)...
	       - a  *log((sqrt(a^2+b^2)+a)/b)...
	       +2*sqrt(a^2+b^2) + aa/2 - 3/4*(a+b)...
	       -0.5*sqrt(4*(b+a)^2+aa^2) )
  error= (laa-lsnow)/lsnow

% Ref Grover Pg 60 Eqn 58: [where sinh^-1(x)= log(x+sqrt(x^2+1))]
% See derive: rec_ind.mth where both are programmed: definitly different

 lgrover= 4e-7*(a*log(2*a/aa)+b*log(2*b/aa)+2*sqrt(a^2+b^2)...
           -a*log(a/b+sqrt((a/b)^2+1))-b*log(b/a+sqrt((b/a)^2+1))...
           -2*(a+b) + (a+b)/4)
