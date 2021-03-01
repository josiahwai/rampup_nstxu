  function [rp,zp,phi,br,bz,bt]=  mag_fil_pt(xp,yp,zp,x1,y1,z1,x2,y2,z2,cc)
% 
% mag_fil_pt calculates the 3-d magrnetic field from a set of filaments
%
% SYNTAX:
%         [rp,zp,phi,br,bz,bt]=  mag_fil_pt(xp,yp,zp,x1,y1,z1,x2,y2,z2,cc)
% INPUT:
%         xp,yp,zp   = coordinates of field points [m]
%         x1,y1...z2 = coordinates of starting (1) and end (2) of filament [m]
%         cc         = current in filament going from pt 1 to pt 2 [A]
%
%         NOTE: all input can be vectors
%
% OUTPUT:
%         rp,zp      = location of field point [m]
%         phi        = toroidal location of field point [deg 0-2pi]
%         br,bz,bt   = B-field at field point radial,axial,toroidal [T]
%     

% Jim Leuer 9-26-99
% taken from tf_yiirot.m 
%
% ----------------------------------------------------------------------
% initiation stuff
% ----------------------------------------------------------------------
   xc= 0.5*(x1+x2); % filament centers
   yc= 0.5*(y1+y2);
   zc= 0.5*(z1+z2);
   nf= length(x1);
   np= length(xp);
% ----------------------------------------------------------------------
% START CALCULATE FIELD AT FIELD POINTS FROM ALL FILAMENTS
% ----------------------------------------------------------------------

  bx= zeros(size(xp));
  by= zeros(size(xp));
  bz= zeros(size(xp));

% ------------------------------------
  for ii= 1:nf % Loop over all filaments
    dx= x2(ii) - x1(ii);
    dy= y2(ii) - y1(ii);
    dz= z2(ii) - z1(ii);
    dl= norm([dx dy dz]); % length of filament
    dir_cos= [dx/dl; dy/dl; dz/dl];
    rot= z_dircos(dir_cos); % rotational transform X_local= rot*X_glob
    rot_t= rot';
    xcc= [xc(ii); yc(ii); zc(ii)];
    ccc= cc(ii);
% ------------------------------------
%  field point infomation
    xpp= [xp-xc(ii), yp-yc(ii), zp-zc(ii)]; % trans x,y,z
    xpp= (rot*xpp')';                   % local to filament x,y,z [np,3]
    dll= dl*ones(np,1);                 % Length of filament
    rfil= dll*0.001;                    % radius of filament for r=0 protect
%    [b_x b_y b_z]= fil_3d_b(xpp(:,1),xpp(:,2),xpp(:,3),dll); 
    [b_x b_y b_z]= filr_3d_b(xpp(:,1),xpp(:,2),xpp(:,3),dll,rfil); %r=0protect
    b= ccc*rot_t*[b_x, b_y, b_z]';
    b= b';
    bx= bx + b(:,1);      
    by= by + b(:,2);      
    bz= bz + b(:,3);      
  end % filament loop
% ------------------------------------
% find radial and toroidal components
  phi= atan2(yp,xp);          % toroidal angle (-pi<thet<=pi)
  id= find(phi<0);
  phi(id)= phi(id)+2*pi;          % make toroidal angle (0<thet<=2pi)
  ct= cos(phi); st= sin(phi);      
  phi= phi*180/pi;           % change to degrees
  rp= sqrt(yp.^2+xp.^2);  % field point radius
  br= bx.*ct + by.*st;
  bt= by.*ct - bx.*st;

% ----------------------------------------------------------------------

  return

% ----------------------------------------------------------------------
% Test Stuff Below
   xp= [0;1];   yp= [2;3];  zp= [4;5];
   x1= [6;7;8];   y1= [9;10;11]; z1= [12;13;14];
   x2= [16;17;18];   y2= [19;20;21]; z2= [22;23;24];
   cc= [26;27;28];

   [rp,zp,phi,br,bz,bt]=  mag_fil_pt(xp,yp,zp,x1,y1,z1,x2,y2,z2,cc)

% test using tf_rad.m as good results
   fld_x= xp;  fld_y= yp;  fld_z= zp;  
   fld_n= length(fld_x);
   fil_x1= x1;     fil_y1= y1;     fil_z1= z1;     
   fil_x2= x2;     fil_y2= y2;     fil_z2= z2;     
   fil_xc= 0.5*(x1+x2);  fil_yc= 0.5*(y1+y2);  fil_zc= 0.5*(z1+z2);  
   fil_i= cc;
   fil_n= length(x1);

% invoke code at bottom of tf_rad.m

   format long

   disp([fld_r fld_z fld_thet])
   disp([fld_br fld_bz fld_bt])

   [rp,zp,phi,br,bz,bt]=  mag_fil_pt(xp,yp,zp,x1,y1,z1,x2,y2,z2,cc);

   disp([rp zp phi])
   disp([br bz bt])

% exact match => looks good JAL 9-28-99

   
