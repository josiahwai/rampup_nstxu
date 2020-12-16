function [r,z,n] = gscontour22(y,r0,z0,np,ra,za)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   [r,z,n] = gscontour22(y,r0,z0,np,ra,za)
%
%  PURPOSE: Find closed contour in matrix y passing through r0, z0
%
%  INPUTS:  y, matrix to find contours in
%           r0, z0, starting point (units are floating index)
%           np, length of output vectors
%           ra, za, contour must enclose ra,za if included and ~nan
%
%  OUTPUTS: r,z, vectors of length np, beginning with r0, z0,
%           then n-1 points along contour back to r0, z0, then nans
%           If no closed contour exists r(1), z(1) = nan
%           gs_interp2(y,r,z) = gs_interp2(y,r0,z0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  METHOD: The solution to cubic equations avoids complex numbers
%          since not allowed in Simulink matlab function blocks
%          Most work done in parallel and dots connected at the end
%	
	
%  WRITTEN BY:  Anders Welander ON 2016-04-17
%
%  MODIFICATION HISTORY:
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% All variables have fixed sizes, meant to work from matlab function blocks

% Number of valid contour points
n = 0;

% For interpolation
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;

% Number of iterations to determine cubic roots
ni = 5;

% Size of y
[nz, nr] = size(y);

% length of output vectors
if nargin < 4
  np = 8*(nr+nz);
end

% r of axis
if nargin < 5
  ra = nan;
end

% z of axis
if nargin < 6
  za = nan;
end

% Memory allocation for contour
r = nan(np,1);
z = nan(np,1);

% The grid for y
ngg = nz*nr;
rgg = ones(nz,1)*(1:nr);
zgg = (1:nz)'*ones(1,nr);

% A 2x2 denser grid for the contour points
nz2 = 2*nz-1;
nr2 = 2*nr-1;
rgg2 = zeros(nz2,nr2);
rgg2(1:2:nz2,1:2:nr2) = rgg;
rgg2(2:2:nz2,1:2:nr2) = rgg(1:nz-1,:);
rgg2(1:2:nz2,2:2:nr2) = rgg(:,1:nr-1);
rgg2(2:2:nz2,2:2:nr2) = rgg(1:nz-1,1:nr-1);
zgg2 = zeros(nz2,nr2);
zgg2(1:2:nz2,1:2:nr2) = zgg;
zgg2(2:2:nz2,1:2:nr2) = zgg(1:nz-1,:);
zgg2(1:2:nz2,2:2:nr2) = zgg(:,1:nr-1);
zgg2(2:2:nz2,2:2:nr2) = zgg(1:nz-1,1:nr-1);
ngg2 = nz2*nr2;

% Values of y at zg = 1:0.5:nz
yh = zeros(nz2,nr);
yh(1:2:nz2,:) = y;
% yh(2,:) = -0.0625*y(0,:)+0.5625*y(1,:)+0.5625*y(2,:)-0.0625*y(3,:);
% y(0,:) = 3*y(1,:)-3*y(2,:)+y(3,:);
yh(2,:) = 0.375*y(1,:)+0.75*y(2,:)-0.125*y(3,:);
yh(4:2:nz2-3,:) = -0.0625*y(1:nz-3,:)+0.5625*y(2:nz-2,:)+...
                   0.5625*y(3:nz-1,:)-0.0625*y(4:nz,:);
% yh(nz2-1,:) = -0.0625*y(nz-2,:)+0.5625*y(nz-1,:)+0.5625*y(nz,:)-0.0625*y(nz+1,:);
% y(nz+1,:) = 3*y(nz,:)-3*y(nz-1,:)+y(nz-2,:);
yh(nz2-1,:) = -0.125*y(nz-2,:)+0.75*y(nz-1,:)+0.375*y(nz,:);

% Values of y at rg = 1:0.5:nr
yv = zeros(nz,nr2);
yv(:,1:2:nr2) = y;
yv(:,2) = 0.375*y(:,1)+0.75*y(:,2)-0.125*y(:,3);
yv(:,4:2:nr2-3) = -0.0625*y(:,1:nr-3,:)+0.5625*y(:,2:nr-2)+...
                   0.5625*y(:,3:nr-1)-0.0625*y(:,4:nr,:);
yv(:,nr2-1) = -0.125*y(:,nr-2)+0.75*y(:,nr-1)+0.375*y(:,nr);

% Find local min or max on horizontal grid lines
Ah = zeros(nz2,nr);
Bh = zeros(nz2,nr);
Ch = zeros(nz2,nr);
Dh = zeros(nz2,nr);
HA = [nan(nz2,1) [nan(1,nr-3);-yh(2:nz2-1,1:nr-3)+3*yh(2:nz2-1,2:nr-2)-...
  3*yh(2:nz2-1,3:nr-1)+yh(2:nz2-1,4:nr); nan(1,nr-3)] nan(nz2,2)]/2;
HB = [nan(nz2,1) [nan(1,nr-3); 2*yh(2:nz2-1,1:nr-3)-5*yh(2:nz2-1,2:nr-2)+...
  4*yh(2:nz2-1,3:nr-1)-yh(2:nz2-1,4:nr); nan(1,nr-3)] nan(nz2,2)]/2;
HC = [nan(nz2,1) [nan(1,nr-3); -yh(2:nz2-1,1:nr-3)+yh(2:nz2-1,3:nr-1); ...
  nan(1,nr-3)] nan(nz2,2)]/2;
HR1 = zeros(nz2,nr);
HR2 = ones(nz2,nr);
HD = -HB./HA/3;
HS = HD.^2-HC./HA/3;
h = HS >= 0; % SPECIAL CASE HS==0
HS(h) = sqrt(HS(h));
HR1(h) = HD(h)-HS(h);
HR2(h) = HD(h)+HS(h);
HR1(HR1 < 0 | HR1 >= 1) = 0;
HR2(HR2 <= 0 | HR2 > 1) = 1;
H0 = yh;
H1 = yh+HC.*HR1+HB.*HR1.^2+HA.*HR1.^3;
H2 = yh+HC.*HR2+HB.*HR2.^2+HA.*HR2.^3;
H3 = [yh(:,2:nr) nan(nz2,1)];
Hm = yh+HC./2+HB/4+HA/8;

% Find local min or max on vertical grid lines
Av = zeros(nz,nr2);
Bv = zeros(nz,nr2);
Cv = zeros(nz,nr2);
Dv = zeros(nz,nr2);
VA = [nan(1,nr2); [nan(nz-3,1),-yv(1:nz-3,2:nr2-1)+3*yv(2:nz-2,2:nr2-1)-...
  3*yv(3:nz-1,2:nr2-1)+yv(4:nz,2:nr2-1), nan(nz-3,1)]; nan(2,nr2)]/2;
VB = [nan(1,nr2); [nan(nz-3,1),2*yv(1:nz-3,2:nr2-1)-5*yv(2:nz-2,2:nr2-1)+...
  4*yv(3:nz-1,2:nr2-1)-yv(4:nz,2:nr2-1), nan(nz-3,1)]; nan(2,nr2)]/2;
VC = [nan(1,nr2); [nan(nz-3,1),-yv(1:nz-3,2:nr2-1)+yv(3:nz-1,2:nr2-1), ...
  nan(nz-3,1)]; nan(2,nr2)]/2;
VZ1 = zeros(nz,nr2);
VZ2 = ones(nz,nr2);
VD = -VB./VA/3;
VS = VD.^2-VC./VA/3;
v = VS >= 0;
VS(v) = sqrt(VS(v));
VZ1(v) = VD(v)-VS(v);
VZ2(v) = VD(v)+VS(v);
VZ1(VZ1 < 0 | VZ1 >= 1) = 0;
VZ2(VZ2 <= 0 | VZ2 > 1) = 1;
V0 = yv;
V1 = yv+VC.*VZ1+VB.*VZ1.^2+VA.*VZ1.^3;
V2 = yv+VC.*VZ2+VB.*VZ2.^2+VA.*VZ2.^3;
V3 = [yv(2:nz,:); nan(1,nr2)];
Vm = yv+VC./2+VB/4+VA/8;

% Translate exit point from cell to entry point for adjacent cell
oi = [9 8 7 12 11 10 3 2 1 6 5 4];
op = [-1 -1 -1 -nz2 -nz2 -nz2 1 1 1 nz2 nz2 nz2];
dr2 = [0 0 0 -1 -1 -1 0 0 0 1 1 1];
dz2 = [-1 -1 -1 0 0 0 1 1 1 0 0 0];
% exit F(p,i) is entry F(p+op(i),oi(i))

% Half of small quadrangle width
d = 1e-5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Done with prerequisites for contouring %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% The starting point
iz = max(2,min(nz-2,floor(z0))); % grid line below z0
ir = max(2,min(nr-2,floor(r0))); % grid line inside r0
tz = z0-iz;
tr = r0-ir;

% For interpolation to starting point r0, z0
wz0 = [1 tz tz^2 tz^3]*mx;
wr0 = [1 tr tr^2 tr^3]*mx;
wz1 = [0 1 2*tz 3*tz^2]*mx;
wr1 = [0 1 2*tr 3*tr^2]*mx;
y44 = y(iz-1:iz+2,ir-1:ir+2);

% The value to contour and r,z derivatives of y at r0, z0 (floating index units)
yc = wz0*y44*wr0';
yr = wz0*y44*wr1';
yz = wz1*y44*wr0';

% Index to cell containing r0, z0 in denser grid
p1 = floor(2*z0-1) + nz2*floor(2*r0-2);

% Contour points on small quad around r0, z0
s1 = init(y44,yc,1,tr+[-1 1]*d,tz+[-1 1]*d);

% Number of contours coming out from starting point
n1 = sum(double(s1.f));

% The full cell k is the region [rgg(k), rgg(k)+1], [zgg(k), zgg(k)+1]


% Position of contour points on 2x1 denser grid relative to inner edge of full cell

Rc1 = nan(nz2,nr);
h = H0 >= yc & yc > H1 | H0 <= yc & yc < H1;
Rc1(h) = (yc-H0(h))./(H1(h)-H0(h)).*HR1(h);
for i = 1:ni
  Ah(h) = 6*HA(h).*Rc1(h)+2*HB(h);
  Bh(h) = 3*HA(h).*Rc1(h).^2+2*HB(h).*Rc1(h)+HC(h);
  Ch(h) = HA(h).*Rc1(h).^3+HB(h).*Rc1(h).^2+HC(h).*Rc1(h)+H0(h)-yc;
  Dh(h) = -Bh(h)./Ah(h);
  HS(h) = Dh(h).^2-2*Ch(h)./Ah(h);
  
  % Normal fast convergence for smooth variations in y
  g = h & abs(Ah) >= 1e-4 & HS >= 0;
  Rc1(g) = Rc1(g)+Dh(g)-sign(Dh(g)).*sqrt(HS(g));
  
  % Ch+yc instead of H0, Rc1 instead of 0, for linear interpolation
  g = h & (abs(Ah) < 1e-4 | HS < 0) & (H0-yc).*Ch >= 0;
  Rc1(g) = Rc1(g)-Ch(g)./(H1(g)-Ch(g)-yc).*(HR1(g)-Rc1(g));
  
  % Ch+yc instead of H1, Rc1 instead of HR1, for linear interpolation
  g = h & (abs(Ah) < 1e-4 | HS < 0) & (H0-yc).*Ch < 0;
  Rc1(g) = (yc-H0(g))./(Ch(g)+yc-H0(g)).*Rc1(g);
end
% In case some point still goes bad
h = h & (Rc1 < 0 | Rc1 > HR1 | abs(Ch./(H1-H0)) > 1e-3);
if any(h(:))
  Ah(h) = 0;
  Bh(h) = HR1(h);
  Ch(h) = H0(h);
  Dh(h) = H1(h);
  for i = 1:30 % Bisection method
    Rc1(h) = (Ah(h)+Bh(h))/2;
    HS(h) = HA(h).*Rc1(h).^3+HB(h).*Rc1(h).^2+HC(h).*Rc1(h)+H0(h);
    g = h & (Ch-yc).*(HS-yc) >= 0;
    Ah(g) = Rc1(g);
    Ch(g) = HS(g);
    g = h & (Ch-yc).*(HS-yc) < 0;
    Bh(g) = Rc1(g);
    Dh(g) = HS(g);
  end
end

Rc2 = nan(nz2,nr);
h = H1 >= yc & yc > H2 | H1 <= yc & yc < H2;
Rc2(h) = HR1(h)+(yc-H1(h))./(H2(h)-H1(h)).*(HR2(h)-HR1(h));
for i = 1:ni
  Ah(h) = 6*HA(h).*Rc2(h)+2*HB(h);
  Bh(h) = 3*HA(h).*Rc2(h).^2+2*HB(h).*Rc2(h)+HC(h);
  Ch(h) = HA(h).*Rc2(h).^3+HB(h).*Rc2(h).^2+HC(h).*Rc2(h)+H0(h)-yc;
  Dh(h) = -Bh(h)./Ah(h);
  HS(h) = Dh(h).^2-2*Ch(h)./Ah(h);
  
  % Normal fast convergence for smooth variations in y
  g = h & abs(Ah) >= 1e-4 & HS >= 0;
  Rc2(g) = Rc2(g)+Dh(g)-sign(Dh(g)).*sqrt(HS(g));
  
  % Ch+yc instead of H1, Rc2 instead of HR1, for linear interpolation
  g = h & (abs(Ah) < 1e-4 | HS < 0) & (H1-yc).*Ch >= 0;
  Rc2(g) = Rc2(g)-Ch(g)./(H2(g)-Ch(g)-yc).*(HR2(g)-Rc2(g));
  
  % Ch+yc instead of H2, Rc2 instead of HR2, for linear interpolation
  g = h & (abs(Ah) < 1e-4 | HS < 0) & (H1-yc).*Ch < 0;
  Rc2(g) = HR1(g)+(yc-H1(g))./(Ch(g)+yc-H1(g)).*(Rc2(g)-HR1(g));
end
% In case some point still goes bad
h = h & (Rc2 < HR1 | Rc2 > HR2 | abs(Ch./(H2-H1)) > 1e-3);
if any(h(:))
  Ah(h) = HR1(h);
  Bh(h) = HR2(h);
  Ch(h) = H1(h);
  Dh(h) = H2(h);
  for i = 1:30 % Bisection method
    Rc2(h) = (Ah(h)+Bh(h))/2;
    HS(h) = HA(h).*Rc2(h).^3+HB(h).*Rc2(h).^2+HC(h).*Rc2(h)+H0(h);
    g = h & (Ch-yc).*(HS-yc) >= 0;
    Ah(g) = Rc2(g);
    Ch(g) = HS(g);
    g = h & (Ch-yc).*(HS-yc) < 0;
    Bh(g) = Rc2(g);
    Dh(g) = HS(g);
  end
end

Rc3 = nan(nz2,nr);
h = H2 >= yc & yc > H3 | H2 <= yc & yc < H3;
Rc3(h) = HR2(h)+(yc-H2(h))./(H3(h)-H2(h)).*(1-HR2(h));
for i = 1:ni
  Ah(h) = 6*HA(h).*Rc3(h)+2*HB(h);
  Bh(h) = 3*HA(h).*Rc3(h).^2+2*HB(h).*Rc3(h)+HC(h);
  Ch(h) = HA(h).*Rc3(h).^3+HB(h).*Rc3(h).^2+HC(h).*Rc3(h)+H0(h)-yc;
  Dh(h) = -Bh(h)./Ah(h);
  HS(h) = Dh(h).^2-2*Ch(h)./Ah(h);
  
  % Normal fast convergence for smooth variations in y
  g = h & abs(Ah) >= 1e-4 & HS >= 0;
  Rc3(g) = Rc3(g)+Dh(g)-sign(Dh(g)).*sqrt(HS(g));
  
  % Ch+yc instead of H2, Rc3 instead of HR2, for linear interpolation
  g = h & (abs(Ah) < 1e-4 | HS < 0) & (H2-yc).*Ch >= 0;
  Rc3(g) = Rc3(g)-Ch(g)./(H3(g)-Ch(g)-yc).*(1-Rc3(g));
  
  % Ch+yc instead of H3, Rc3 instead of 1, for linear interpolation
  g = h & (abs(Ah) < 1e-4 | HS < 0) & (H2-yc).*Ch < 0;
  Rc3(g) = HR2(g)+(yc-H2(g))./(Ch(g)+yc-H2(g)).*(Rc3(g)-HR2(g));
end
% In case some point still goes bad
h = h & (Rc3 < HR2 | Rc3 > 1 | abs(Ch./(H3-H2)) > 1e-3);
if any(h(:))
  Ah(h) = HR2(h);
  Bh(h) = 1;
  Ch(h) = H2(h);
  Dh(h) = H3(h);
  for i = 1:30 % Bisection method
    Rc3(h) = (Ah(h)+Bh(h))/2;
    HS(h) = HA(h).*Rc3(h).^3+HB(h).*Rc3(h).^2+HC(h).*Rc3(h)+H0(h);
    g = h & (Ch-yc).*(HS-yc) >= 0;
    Ah(g) = Rc3(g);
    Ch(g) = HS(g);
    g = h & (Ch-yc).*(HS-yc) < 0;
    Bh(g) = Rc3(g);
    Dh(g) = HS(g);
  end
end


% Position of contour points on 1x2 denser grid relative to lower edge of full cell

Zc1 = nan(nz,nr2);
v = V0 >= yc & yc > V1 | V0 <= yc & yc < V1;
Zc1(v) = (yc-V0(v))./(V1(v)-V0(v)).*VZ1(v);
for i = 1:ni
  Av(v) = 6*VA(v).*Zc1(v)+2*VB(v);
  Bv(v) = 3*VA(v).*Zc1(v).^2+2*VB(v).*Zc1(v)+VC(v);
  Cv(v) = VA(v).*Zc1(v).^3+VB(v).*Zc1(v).^2+VC(v).*Zc1(v)+V0(v)-yc;
  Dv(v) = -Bv(v)./Av(v);
  VS(v) = Dv(v).^2-2*Cv(v)./Av(v);

  % Normal fast convergence for smooth variations in y
  w = v & abs(Av) >= 1e-4 & VS >= 0;
  Zc1(w) = Zc1(w)+Dv(w)-sign(Dv(w)).*sqrt(VS(w));
  
  % Cv+yc instead of V0, Zc1 instead of 0, for linear interpolation
  w = v & (abs(Av) < 1e-4 | VS < 0) & (V0-yc).*Cv >= 0;
  Zc1(w) = Zc1(w)-Cv(w)./(V1(w)-Cv(w)-yc).*(VZ1(w)-Zc1(w));
  
  % Cv+yc instead of V1, Zc1 instead of VZ1, for linear interpolation
  w = v & (abs(Av) < 1e-4 | VS < 0) & (V0-yc).*Cv < 0;
  Zc1(w) = (yc-V0(w))./(Cv(w)+yc-V0(w)).*Zc1(w);
end
% In case some point still goes bad
v = v & (Zc1 < 0 | Zc1 > VZ1 | abs(Cv./(V1-V0)) > 1e-3);
if any(v(:))
  Av(v) = 0;
  Bv(v) = VZ1(v);
  Cv(v) = V0(v);
  Dv(v) = V1(v);
  for i = 1:30 % Bisection method
    Zc1(v) = (Av(v)+Bv(v))/2;
    VS(v) = VA(v).*Zc1(v).^3+VB(v).*Zc1(v).^2+VC(v).*Zc1(v)+V0(v);
    w = v & (Cv-yc).*(VS-yc) >= 0;
    Av(w) = Zc1(w);
    Cv(w) = VS(w);
    w = v & (Cv-yc).*(VS-yc) < 0;
    Bv(w) = Zc1(w);
    Dv(w) = VS(w);
  end
end

Zc2 = nan(nz,nr2);
v = V1 >= yc & yc > V2 | V1 <= yc & yc < V2;
Zc2(v) = VZ1(v)+(yc-V1(v))./(V2(v)-V1(v)).*(VZ2(v)-VZ1(v));
for i = 1:ni
  Av(v) = 6*VA(v).*Zc2(v)+2*VB(v);
  Bv(v) = 3*VA(v).*Zc2(v).^2+2*VB(v).*Zc2(v)+VC(v);
  Cv(v) = VA(v).*Zc2(v).^3+VB(v).*Zc2(v).^2+VC(v).*Zc2(v)+V0(v)-yc;
  Dv(v) = -Bv(v)./Av(v);
  VS(v) = Dv(v).^2-2*Cv(v)./Av(v);
  
  % Normal fast convergence for smooth variations in y
  w = v & abs(Av) >= 1e-4 & VS >= 0;
  Zc2(w) = Zc2(w)+Dv(w)-sign(Dv(w)).*sqrt(VS(w));
  
  % Cv+yc instead of V1, Zc2 instead of VZ1, for linear interpolation
  w = v & (abs(Av) < 1e-4 | VS < 0) & (V1-yc).*Cv >= 0;
  Zc2(w) = Zc2(w)-Cv(w)./(V2(w)-Cv(w)-yc).*(VZ2(w)-Zc2(w));
  
  % Cv+yc instead of V2, Zc2 instead of VZ2, for linear interpolation
  w = v & (abs(Av) < 1e-4 | VS < 0) & (V1-yc).*Cv < 0;
  Zc2(w) = VZ1(w)+(yc-V1(w))./(Cv(w)+yc-V1(w)).*(Zc2(w)-VZ1(w));
end
% In case some point still goes bad
v = v & (Zc2 < VZ1 | Zc2 > VZ2 | abs(Cv./(V2-V1)) > 1e-3);
if any(v(:))
  Av(v) = VZ1(v);
  Bv(v) = VZ2(v);
  Cv(v) = V1(v);
  Dv(v) = V2(v);
  for i = 1:30 % Bisection method
    Zc2(v) = (Av(v)+Bv(v))/2;
    VS(v) = VA(v).*Zc2(v).^3+VB(v).*Zc2(v).^2+VC(v).*Zc2(v)+V0(v);
    w = v & (Cv-yc).*(VS-yc) >= 0;
    Av(w) = Zc2(w);
    Cv(w) = VS(w);
    w = v & (Cv-yc).*(VS-yc) < 0;
    Bv(w) = Zc2(w);
    Dv(w) = VS(w);
  end
end

Zc3 = nan(nz,nr2);
v = V2 >= yc & yc > V3 | V2 <= yc & yc < V3;
Zc3(v) = VZ2(v)+(yc-V2(v))./(V3(v)-V2(v)).*(1-VZ2(v));
for i = 1:ni
  Av(v) = 6*VA(v).*Zc3(v)+2*VB(v);
  Bv(v) = 3*VA(v).*Zc3(v).^2+2*VB(v).*Zc3(v)+VC(v);
  Cv(v) = VA(v).*Zc3(v).^3+VB(v).*Zc3(v).^2+VC(v).*Zc3(v)+V0(v)-yc;
  Dv(v) = -Bv(v)./Av(v);
  VS(v) = Dv(v).^2-2*Cv(v)./Av(v);
  
  % Normal fast convergence for smooth variations in y
  w = v & abs(Av) >= 1e-4 & VS >= 0;
  Zc3(w) = Zc3(w)+Dv(w)-sign(Dv(w)).*sqrt(VS(w));
  
  % Cv+yc instead of V2, Zc3 instead of VZ2, for linear interpolation
  w = v & (abs(Av) < 1e-4 | VS < 0) & (V2-yc).*Cv >= 0;
  Zc3(w) = Zc3(w)-Cv(w)./(V3(w)-Cv(w)-yc).*(1-Zc3(w));
  
  % Cv+yc instead of V3, Zc3 instead of 1, for linear interpolation
  w = v & (abs(Av) < 1e-4 | VS < 0) & (V2-yc).*Cv < 0;
  Zc3(w) = VZ2(w)+(yc-V2(w))./(Cv(w)+yc-V2(w)).*(Zc3(w)-VZ2(w));
end
% In case some point still goes bad
v = v & (Zc3 < VZ2 | Zc3 > 1 | abs(Cv./(V3-V2)) > 1e-3);
if any(v(:))
  Av(v) = VZ2(v);
  Bv(v) = 1;
  Cv(v) = V2(v);
  Dv(v) = V3(v);
  for i = 1:30 % Bisection method
    Zc3(v) = (Av(v)+Bv(v))/2;
    VS(v) = VA(v).*Zc3(v).^3+VB(v).*Zc3(v).^2+VC(v).*Zc3(v)+V0(v);
    w = v & (Cv-yc).*(VS-yc) >= 0;
    Av(w) = Zc3(w);
    Cv(w) = VS(w);
    w = v & (Cv-yc).*(VS-yc) < 0;
    Bv(w) = Zc3(w);
    Dv(w) = VS(w);
  end
end

% Insert all points in 2x2 denser grid
R1 = nan(nz2,nr2);
R1(:,1:2:nr2) = Rc1;
R1(:,2:2:nr2) = Rc1(:,1:nr-1);
R1b = [R1(2:nz2,:); nan(1,nr2)];
R2 = nan(nz2,nr2);
R2(:,1:2:nr2) = Rc2;
R2(:,2:2:nr2) = Rc2(:,1:nr-1);
R2b = [R2(2:nz2,:); nan(1,nr2)];
R3 = nan(nz2,nr2);
R3(:,1:2:nr2) = Rc3;
R3(:,2:2:nr2) = Rc3(:,1:nr-1);
R3b = [R3(2:nz2,:); nan(1,nr2)];
Z1 = nan(nz2,nr2);
Z1(1:2:nz2,:) = Zc1;
Z1(2:2:nz2,:) = Zc1(1:nz-1,:);
Z1b = [Z1(:,2:nr2) nan(nz2,1)];
Z2 = nan(nz2,nr2);
Z2(1:2:nz2,:) = Zc2;
Z2(2:2:nz2,:) = Zc2(1:nz-1,:);
Z2b = [Z2(:,2:nr2) nan(nz2,1)];
Z3 = nan(nz2,nr2);
Z3(1:2:nz2,:) = Zc3;
Z3(2:2:nz2,:) = Zc3(1:nz-1,:);
Z3b = [Z3(:,2:nr2) nan(nz2,1)];

% For grid lines in 2x2 denser grid
Ri = zeros(nz2,nr2);
Ri(:,2:2:nr2) = 0.5;
Ro = ones(nz2,nr2);
Ro(:,1:2:nr2) = 0.5;
Zl = zeros(nz2,nr2);
Zl(2:2:nz2,:) = 0.5;
Zu = ones(nz2,nr2);
Zu(1:2:nz2,:) = 0.5;

% Three regions may contain contour points along each line surrounding a cell
% These are sorted beginning with the lower outer corner and going clock-wise

% Position of contour points around small cells relative to lower inner edge of full cell
R = [R3(:) R2(:) R1(:) Ri(:) Ri(:) Ri(:) R1b(:) R2b(:) R3b(:) Ro(:) Ro(:) Ro(:)];
Z = [Zl(:) Zl(:) Zl(:) Z1(:) Z2(:) Z3(:) Zu(:) Zu(:) Zu(:) Z3b(:) Z2b(:) Z1b(:)];
% Invalid points are either outside the range for the small cell or contain a nan

% Flags walking from lower right corner clockwise around each of the small cells
F = false(ngg2,12);

% flags for one region at a time, all cells in denser grid
f = false(nz2,nr2);

% G will contain values that depend on which flags are set
G = zeros(ngg2,1);

% Set all flags for regions with contour points around the edges of small cells

f(:,1:2:nr2) = Rc3 > 0 & Rc3 <= 0.5;
f(:,2:2:nr2) = Rc3(:,1:nr-1) > 0.5 & Rc3(:,1:nr-1) <= 1;
F(:,1) = f(:);
G(f(:)) = G(f(:))+1;

f(:,1:2:nr2) = Rc2 > 0 & Rc2 <= 0.5;
f(:,2:2:nr2) = Rc2(:,1:nr-1) > 0.5 & Rc2(:,1:nr-1) <= 1;
F(:,2) = f(:);
G(f(:)) = G(f(:))+2;

f(:,1:2:nr2) = Rc1 > 0 & Rc1 <= 0.5;
f(:,2:2:nr2) = Rc1(:,1:nr-1) > 0.5 & Rc1(:,1:nr-1) <= 1;
F(:,3) = f(:);
G(f(:)) = G(f(:))+4;

f(1:2:nz2,:) = Zc1 >= 0 & Zc1 < 0.5;
f(2:2:nz2,:) = Zc1(1:nz-1,:) >= 0.5 & Zc1(1:nz-1,:) < 1;
F(:,4) = f(:);
G(f(:)) = G(f(:))+8;

f(1:2:nz2,:) = Zc2 >= 0 & Zc2 < 0.5;
f(2:2:nz2,:) = Zc2(1:nz-1,:) >= 0.5 & Zc2(1:nz-1,:) < 1;
F(:,5) = f(:);
G(f(:)) = G(f(:))+16;

f(1:2:nz2,:) = Zc3 >= 0 & Zc3 < 0.5;
f(2:2:nz2,:) = Zc3(1:nz-1,:) >= 0.5 & Zc3(1:nz-1,:) < 1;
F(:,6) = f(:);
G(f(:)) = G(f(:))+32;

f(1:nz2-1,1:2:nr2) = Rc1(2:nz2,:) >= 0 & Rc1(2:nz2,:) < 0.5;
f(1:nz2-1,2:2:nr2) = Rc1(2:nz2,1:nr-1) >= 0.5 & Rc1(2:nz2,1:nr-1) < 1;
F(:,7) = f(:);
G(f(:)) = G(f(:))+64;

f(1:nz2-1,1:2:nr2) = Rc2(2:nz2,:) >= 0 & Rc2(2:nz2,:) < 0.5;
f(1:nz2-1,2:2:nr2) = Rc2(2:nz2,1:nr-1) >= 0.5 & Rc2(2:nz2,1:nr-1) < 1;
F(:,8) = f(:);
G(f(:)) = G(f(:))+128;

f(1:nz2-1,1:2:nr2) = Rc3(2:nz2,:) >= 0 & Rc3(2:nz2,:) < 0.5;
f(1:nz2-1,2:2:nr2) = Rc3(2:nz2,1:nr-1) >= 0.5 & Rc3(2:nz2,1:nr-1) < 1;
F(:,9) = f(:);
G(f(:)) = G(f(:))+256;

f(1:2:nz2,1:nr2-1) = Zc3(:,2:nr2) <= 0.5 & Zc3(:,2:nr2) > 0;
f(2:2:nz2,1:nr2-1) = Zc3(1:nz-1,2:nr2) <= 1 & Zc3(1:nz-1,2:nr2) > 0.5;
F(:,10) = f(:);
G(f(:)) = G(f(:))+512;

f(1:2:nz2,1:nr2-1) = Zc2(:,2:nr2) <= 0.5 & Zc2(:,2:nr2) > 0;
f(2:2:nz2,1:nr2-1) = Zc2(1:nz-1,2:nr2) <= 1 & Zc2(1:nz-1,2:nr2) > 0.5;
F(:,11) = f(:);
G(f(:)) = G(f(:))+1024;

f(1:2:nz2,1:nr2-1) = Zc1(:,2:nr2) <= 0.5 & Zc1(:,2:nr2) > 0;
f(2:2:nz2,1:nr2-1) = Zc1(1:nz-1,2:nr2) <= 1 & Zc1(1:nz-1,2:nr2) > 0.5;
F(:,12) = f(:);
G(f(:)) = G(f(:))+2048;

% The following puts all contour points in Rc(F), Zc(F)
% Rc = rgg2(:)*ones(1,12)+R; Zc = zgg2(:)*ones(1,12)+Z;
% What remains is to connect them in the right order

% Number of contour points surrounding each cell in denser grid
N = sum(double(F)');

% Sparse matrix containing the io vector for all possible G
mio = zeros(2^12,12);
for i = 1:12
  for j = i+1:12
    mio(2^(i-1)+2^(j-1),[i j]) = [j i];
  end
end

% Matrix to connect entry and exit points of the cells
io = zeros(ngg2,12);
io(N==2,:) = mio(G(N==2),:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Connect the dots quickly using the information in io %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prevent rounding error that can occur when r0, z0 is an x-point close to a grid line
i2r0 = round(2*r0-1);
i2z0 = round(2*z0-1);
if n1 > 2 & (abs(2*r0-1-i2r0) < d | abs(2*z0-1-i2z0) < d)
  % Find the contour points on the little quadrangle s1 and use to create points
  rd = nan(1,4);
  zd = nan(1,4);
  n = 0;
  for i = 1:12
    if s1.f(i)
      s1.i = i;
      n = n+1;
      [rd(n), zd(n)] = point(s1);
    end
    if n == 4
      break
    end
  end
  i2r0 = round(2*r0-1);
  i2z0 = round(2*z0-1);
  if abs(2*r0-1-i2r0) < d
    % r0, z0 is very close to inner boundary of cell p
    p = floor(2*z0-1) + nz2*(i2r0-1);
    dr = (2*r0-1-i2r0)/2;
    tz3 = [tz tz nan];
    if rd(3)-rd(1) ~= 0
      tz3(1) = tz + (zd(3)-zd(1))/(rd(3)-rd(1))*dr;
    end
    if rd(4)-rd(2) ~= 0
      tz3(2) = tz + (zd(4)-zd(2))/(rd(4)-rd(2))*dr;
    end
    tz3 = sort(tz3);
    F(p,4:6) = [1 1 0];
    Z(p,4:6) = tz3;
    % Connect 4 & 5
    for i = [3 2 1 12 11 10 9 8]
      if F(p,i)
        io(p,i) = 4;
	io(p,4) = i;
	break
      end
    end
    for i = [7 8 9 10 11 12 1 2]
      if F(p,i)
        io(p,i) = 5;
	io(p,5) = i;
	break
      end
    end
    % Move to the cell inside p
    p = p-nz2;
    % Now r0, z0 is very close to outer boundary of cell p
    F(p,10:12) = [0 1 1];
    Z(p,10:12) = tz3([3 2 1]);
    % Connect 11 & 12
    for i = [1 2 3 4 5 6 7 8]
      if F(p,i)
        io(p,i) = 12;
	io(p,12) = i;
	break
      end
    end
    for i = [9 8 7 6 5 4 3 2]
      if F(p,i)
        io(p,i) = 11;
	io(p,11) = i;
	break
      end
    end
  end
  if abs(z0-round(z0)) < d
    % r0, z0 is very close to lower boundary of cell p
  end
end
f1 = F(p1,:);
% The logic below does not work for tests with r0,z0 near grid line
f1test_possible = abs(2*r0-round(2*r0))>d & abs(2*z0-round(2*z0))>d;
if n1 < N(p1) & f1test_possible
  % Not all f1 connected to starting point and it is possible to find which are
  f1(:) = 0;
  R2 = floor(2*tr)/2+[0 0.5];
  Z2 = floor(2*tz)/2+[0 0.5];
  for i = 1:12 % Check connections from starting point
    if s1.f(i)
      s = init(y44,yc,1,R2,Z2);
      if i < 4 % Go down
	if tz-d > Z2(1)
	  s.z(s.z > Z2(1)) = tz-d;
	else % Point outside p1 because of d shift
	  s.z(s.z == Z2(1)) = tz-d;
	  s.z(s.z > tz-d) = Z2(1);
	end
      elseif i < 7 % Go in
	if tr-d > R2(1)
	  s.r(s.r > R2(1)) = tr-d;
	else
	  s.r(s.r == R2(1)) = tr-d;
	  s.r(s.r > tr-d) = R2(1);
	end
      elseif i < 10 % Go up
	if tz+d < Z2(2)
	  s.z(s.z < Z2(2)) = tz+d;
	else % Point outside p1 because of d shift
	  s.z(s.z == Z2(2)) = tz+d;
	  s.z(s.z < tz+d) = Z2(2);
	end
      else % Go out
	if tr+d < R2(2)
	  s.r(s.r < R2(2)) = tr+d;
	else
	  s.r(s.r == R2(2)) = tr+d;
	  s.r(s.r < tr+d) = R2(2);
	end
      end
      s = update(s);
      s.i = region(s,(s1.r(i)+s1.r(i+1))/2,(s1.z(i)+s1.z(i+1))/2);
      s = getout(s);
      for j = 1:12
        if F(p1,j)
	  if (R(p1,j)-s.r(s.i))*(R(p1,j)-s.r(s.i+1)) <= 0 & ...
	     (Z(p1,j)-s.z(s.i))*(Z(p1,j)-s.z(s.i+1)) <= 0
	    f1(j) = 1;
	  end
	end
      end
    end
  end
end % f1 contains connected points
closed = false;
ongridline = tr == 0 | tr == 0.5 | tz == 0 | tz == 0.5;
if ongridline
  [~,i11] = min((R(p1,:)-tr).^2+(Z(p1,:)-tz).^2);
else
  i11 = 1;
end
for i1 = i11:12
  if f1(i1) % Starting point is connected to grid line point p1,i1
    if ongridline % Starting point on grid line
      r(1) = r0;
      z(1) = z0;
      n = 1; % Number of points that have been archived
    else
      r(1) = r0;
      z(1) = z0;
      r(2) = rgg2(p1)+R(p1,i1); % First point on grid line
      z(2) = zgg2(p1)+Z(p1,i1);
      n = 2;
    end
    p = p1+op(i1); % p now points to the next cell to visit
    i = oi(i1); % Point i in new p is now the same as i1 in p1
    % If the contour passes through an x-point tracing is restored to it to test all paths
    restore.p = 0;
    restore.f = false(1,12);
    restore.n = 0;
    while F(p,i)
      F(p,i) = 0; % Point i, cell p is already archived
      j = io(p,i); % Point i traces to j
      if j == 0 % N(p) must be > 2, find io(p,:)
	ir = rgg2(p);
	iz = zgg2(p);
	if iz < 2 | iz > nz-2 | ir < 2 | ir > nr-2
	  break
	end
	if rgg2(p) == rgg2(p+nz2)
	  R2 = [0 0.5];
	else
	  R2 = [0.5 1];
	end
	if zgg2(p) == zgg2(p+1)
	  Z2 = [0 0.5];
	else
	  Z2 = [0.5 1];
	end
	s = init(y(iz-1:iz+2,ir-1:ir+2),yc,i,R2,Z2);
	s.i = region(s,R(p,i),Z(p,i));
	s = getout(s);
	if s.d < 1e-4 % More than one exit is allowed
	  % Create a restore point so all exits can be tested
	  restore.p = p;
	  restore.f = F(p,:);
	  restore.i = i;
	  restore.n = n;
	  for j = 1:12
	    if restore.f(j)
	      restore.f(j) = 0;
	      break
	    end
	  end
	else % only one exit point and it's in region s.i
	  for j = 1:12
	    if (R(p,j)-s.r(s.i))*(R(p,j)-s.r(s.i+1))<=0 & ...
	       (Z(p,j)-s.z(s.i))*(Z(p,j)-s.z(s.i+1))<=0
	      break
	    end
	  end
	end
      end
      % Archive point j where contour exits cell p
      if n < np
	n = n+1;
	r(n) = rgg2(p)+R(p,j);
	z(n) = zgg2(p)+Z(p,j);
      end
      F(p,j) = 0; % Exit j from cell p has been archived
      p = p+op(j); % p now points to the next cell to visit
      i = oi(j); % Point i in new p is now the same as j in old p
      if p == p1 % If back in the original cell
        closed = f1(i); % then contour is closed if point of entry wasn't visited earlier
	if closed
	  break
	end
      end
      ir = rgg2(p);
      iz = zgg2(p);
      if iz < 2 | iz > nz-2 | ir < 2 | ir > nr-2
	F(p,i) = 0;
      end
      if ~F(p,i) & any(restore.f)
        % restore and try another j
	p = restore.p;
	i = restore.i;
	n = restore.n;
	F(p,i) = 1;
	for j = 1:12
	  if restore.f(j)
	    restore.f(j) = 0;
	    io(p,i) = j;
	    break
	  end
	end        
      end
    end % End of loop: while F(p,i)
    if closed
      if n < np
	n = n+1;
	r(n) = r0;
	z(n) = z0;
      end
      r(n+1:np) = nan;
      z(n+1:np) = nan;
      if isnan(ra+za)
        break % break out of loop: for i1 = i11:12
      else
        if max(r) > ra & min(r) < ra & max(z) > za & min(z) < za
          break % break out of loop: for i1 = i11:12
        else
	  r(1) = nan; % Would be r0, z0 but are nans
          z(1) = nan; % to flag that contour does not enclose ra,za
	end
      end
    else
      r(1) = nan; % Would be r0, z0 but are nans
      z(1) = nan; % to flag that contour not closed
    end
  end
end % End of loop: for i1 = i11:12

function s = init(y,v,i,R,Z)
s.y = y; % The 4x4 y-values needed for cubic hermite interpolation
s.v = v; % y value being traced
s.i = i; % Traced point is between r(i:i+1),z(i:i+1)
s.f = false(1,12); % Will be 1 for intervals with y=v point
s.n = 0; % Will be sum(s.f)
s.R = R; % [min(r), max(r)] for rectangle to trace
s.Z = Z; % [min(z), max(z)] for rectangle to trace
% Corners and extremes going clockwise from lower right 
s.r = R(1)+(R(2)-R(1))*[1 1 0 0 0 0 0 0 1 1 1 1 1];
s.z = Z(1)+(Z(2)-Z(1))*[0 0 0 0 0 1 1 1 1 1 1 0 0];
s.p = zeros(1,13); % Will hold y-values at r,z
s.x = false; % Will be true when traced point exits s.R, s.Z
s = update(s);
s.i = i; % Traced point is between r(i:i+1),z(i:i+1)
s.d = inf; % closest approach of contour to discovered x-point inside s.R, s.Z

function s = update(s)
% Updates r,z with points of local extremes in y, then updates p,f,n,i
f = s.f;
r = s.r;
z = s.z;
y = s.y;
v = s.v;
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % for interpolation
e = mx*[y'*mx'*[1 z(1) z(1)^2 z(1)^3]', ...
        y*mx'*[1 r(4) r(4)^2 r(4)^3]', ...  % e = [d;c;b;a]
	y'*mx'*[1 z(7) z(7)^2 z(7)^3]', ...
	y*mx'*[1 r(10) r(10)^2 r(10)^3]'];
d4 = -e(3,:)./e(4,:)/3;
s4 = d4.^2-e(2,:)./e(4,:)/3;
s4(s4<0) = nan;
s4 = sqrt(s4);
t4 = d4-s4;
if t4(1) > r(4) & t4(1) < r(1)
  r(3) = t4(1);
else
  r(3) = r(4);
end
if t4(2) > z(4) & t4(2) < z(7)
  z(5) = t4(2);
else
  z(5) = z(4);
end
if t4(3) > r(7) & t4(3) < r(10)
  r(8) = t4(3);
else
  r(8) = r(7);
end
if t4(4) > z(1) & t4(4) < z(10)
  z(12) = t4(4);
else
  z(12) = z(1);
end
t4 = d4+s4;
if t4(1) > r(4) & t4(1) < r(1)
  r(2) = t4(1);
else
  r(2) = r(1);
end
if t4(2) > z(4) & t4(2) < z(7)
  z(6) = t4(2);
else
  z(6) = z(7);
end
if t4(3) > r(7) & t4(3) < r(10)
  r(9) = t4(3);
else
  r(9) = r(10);
end
if t4(4) > z(1) & t4(4) < z(10)
  z(11) = t4(4);
else
  z(11) = z(10);
end
p = zeros(1,13);
for j = 1:13
  wz0 = [1 z(j) z(j)^2 z(j)^3]*mx;
  wr0 = [1 r(j) r(j)^2 r(j)^3]*mx;
  p(j) = wz0*y*wr0';
end
n = 0;
for i = 1:12
  j = i; % Normally j = i
  if p(i) == v
    % Unusual situation, contour at corner or extreme
    while r(i) == r(j) & z(i) == z(j)
      % Still same place, back more
      if j > 1
	j = j-1;
      else
	j = 12;
      end
    end
  end
  if p(j) < v & v < p(i+1) | p(j) > v & v > p(i+1)
    n = n+1;
    f(i) = 1;
  else
    f(i) = 0;
  end
end
% Make sure s.i is still correct
i = s.i;
if i > 0 & i < 13
  tr = (s.r(i)+s.r(i+1))/2;
  tz = (s.z(i)+s.z(i+1))/2;
else
  tr = nan;
  tz = nan;
end
if ~isnan(i)
  for i = 1:12
    if (r(i)-tr)*(r(i+1)-tr)<=0 & (z(i)-tz)*(z(i+1)-tz)<=0
      break
    end
  end
end
s.i = i;
s.r = r;
s.z = z;
s.p = p;
s.f = f;
s.n = n;

function s = split(s)
% Split cell, the new r,z will only be corners, i will be updated
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % for interpolation
i = s.i; % index to range containing traced point
f = s.f;
r = s.r;
z = s.z;
v = s.v;
y = s.y;
p = s.p;
ro = r(1); % Corners
ri = r(4);
zl = z(4);
zu = z(7);
r2 = r(i:i+1); % range containing traced point
z2 = z(i:i+1);
for j = [2 3 8 9]
  if i==j&f(j-1) | i+1==j&f(j) | i+1>j&any(f(j-1:j)) | i<j&any(f(j-1:j))
    if r(j) < min(r2)
      ri = max(ri,r(j));
    elseif r(j) < max(r2) % r(i)~=r(i+1) so z(i)=z(i+1)
      wz0 = [1 z(i) z(i)^2 z(i)^3]*mx;
      wr0 = [1 r(j) r(j)^2 r(j)^3]*mx;
      y1 = wz0*y*wr0';
      if y1 < v & p(i) < v | y1 > v & p(i) > v
	% Can use r(j) instead of r(i) and still have v within
	if r(i) < r(i+1)
	  ri = max(ri,r(j));
	else
	  ro = min(ro,r(j));
	end
      else
	% Can use r(j) instead of r(i+1)
	if r(i+1) < r(i)
	  ri = max(ri,r(j));
	else
	  ro = min(ro,r(j));
	end
      end
    elseif r(j) > ri
      ro = min(ro,r(j));
    end
  end
end % End of: for j = [2 3 8 9]
for j = [5 6 11 12]
  if i==j&f(j-1) | i+1==j&f(j) | i+1>j&any(f(j-1:j)) | i<j&any(f(j-1:j))
    if z(j) < min(z2)
      zl = max(zl,z(j));
    elseif z(j) < max(z2) % z(i)~=z(i+1) so r(i)=r(i+1)
      wz0 = [1 z(j) z(j)^2 z(j)^3]*mx;
      wr0 = [1 r(i) r(i)^2 r(i)^3]*mx;
      y1 = wz0*y*wr0';
      if y1 < v & p(i) < v | y1 > v & p(i) > v
	% Can use z(j) instead of z(i) and still have v within
	if z(i) < z(i+1)
	  zl = max(zl,z(j));
	else
	  zu = min(zu,z(j));
	end
      else
	% Can use z(j) instead of z(i+1)
	if z(i+1) < z(i)
	  zl = max(zl,z(j));
	else
	  zu = min(zu,z(j));
	end
      end
    elseif z(j) > zl
      zu = min(zu,z(j));
    end
  end
end % End of: for j = [5 6 11 12]
if ro == r(1) & ri == r(4) & zl == z(4) & zu == z(7)
  % No splitting was done, split in half in both dimensions
  if r(i) == r(i+1)
    if r(i) < ro
      ro = (ri+ro)/2; % Here r(i) == ri
    else
      ri = (ri+ro)/2; % Here r(i) == ro
    end
    half = (z(i)+z(i+1))/2;
    wz0 = [1 half half^2 half^3]*mx;
    wr0 = [1 r(i) r(i)^2 r(i)^3]*mx;
    y1 = wz0*y*wr0';
    if y1 < v & p(i) < v | y1 > v & p(i) > v
      % Can use half instead of z(i) and still have v within
      z2(1) = half;
      if z(i) < z(i+1)
        zl = half;
      else
        zu = half;
      end
    else
      % Can use half instead of z(i+1) and still have v within
      z2(2) = half;
      if z(i+1) < z(i)
        zl = half;
      else
        zu = half;
      end
    end
  else
    % z(i) == z(i+1)
    if z(i) < zu
      zu = (zl+zu)/2; % Here z(i) == zl
    else
      zl = (zl+zu)/2; % Here z(i) == zu
    end
    half = (r(i)+r(i+1))/2;
    wz0 = [1 z(i) z(i)^2 z(i)^3]*mx;
    wr0 = [1 half half^2 half^3]*mx;
    y1 = wz0*y*wr0';
    if y1 < v & p(i) < v | y1 > v & p(i) > v
      % Can use half instead of r(i) and still have v within
      r2(1) = half;
      if r(i) < r(i+1)
        ri = half;
      else
        ro = half;
      end
    else
      % Can use half instead of r(i+1) and still have v within
      r2(2) = half;
      if r(i+1) < r(i)
        ri = half;
      else
        ro = half;
      end
    end    
  end
end
r(r<ri) = ri;
r(r>ro) = ro;
z(z<zl) = zl;
z(z>zu) = zu;
r2(r2<ri) = ri;
r2(r2>ro) = ro;
z2(z2<zl) = zl;
z2(z2>zu) = zu;
for i = 1:12
  if all([r2 z2] == [r(i:i+1) z(i:i+1)])
    break
  end
end
s.r = r;
s.z = z;
s.i = i;

function s = walk(s)
% Find index where contour exits (ASSUMES s.n=2)
for i = 1:12
  if s.f(i) & i~=s.i
    break
  end
end
s.i = i; % this index of exit interval may change to index of entry below
% Make new rectangle for continued tracing if not at edge of s.R, s.Z
if i < 4      % Moving down
  if s.z(1) > s.Z(1)
    s.z = s.Z(1)+(s.z(1)-s.Z(1))*[0 0 0 0 0 1 1 1 1 1 1 0 0];
    s.r([8 9]) = s.r([3 2]); % New top = old bottom
    s.i = 10-i;
  else
    s.x = true;
  end
elseif i < 7  % Moving left
  if s.r(4) > s.R(1)
    s.r = s.R(1)+(s.r(4)-s.R(1))*[1 1 0 0 0 0 0 0 1 1 1 1 1];
    s.z([11 12]) = s.z([6 5]); % New right = old left
    s.i = 16-i;
  else
    s.x = true;
  end
elseif i < 10 % Moving up
  if s.z(7) < s.Z(2)
    s.z = s.z(7)+(s.Z(2)-s.z(7))*[0 0 0 0 0 1 1 1 1 1 1 0 0];
    s.r([2 3]) = s.r([9 8]); % New bottom = old top
    s.i = 10-i;
  else
    s.x = true;
  end
else          % Moving right
  if s.r(1) < s.R(2)
    s.r = s.r(1)+(s.R(2)-s.r(1))*[1 1 0 0 0 0 0 0 1 1 1 1 1];
    s.z([5 6]) = s.z([12 11]); % New left = old right
    s.i = 16-i;
  else
    s.x = true;
  end
end

% Walk until reaching the border s.R, s.Z of rectangle
function s = getout(s,mxiter)
if nargin < 2
  mxiter = 9999;
end
count = 0;
while ~s.x & count < mxiter
  % update sub cell
  s = update(s);
  if s.n > 2
    % Decide whether to split cell
    dr = s.r(1)-s.r(4); % lower outer - lower inner corners
    dz = s.z(7)-s.z(4); % upper inner - lower inner corners
    if dr < 1e-4 | dr < 1e-4
      % Find the x-point
       mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % for interpolation
       tr = (s.r(1)+s.r(4))/2;
       tz = (s.z(7)+s.z(4))/2;
       for k = 1:9
	 wz0 = [1 tz tz^2 tz^3]*mx;
	 wr0 = [1 tr tr^2 tr^3]*mx;
	 wz1 = [0 1 2*tz 3*tz^2]*mx;
	 wr1 = [0 1 2*tr 3*tr^2]*mx;
	 wz2 = [0 0 2 6*tz]*mx;
	 wr2 = [0 0 2 6*tr]*mx;
	 y = wz0*s.y*wr0';
	 yr = wz0*s.y*wr1';
	 yz = wz1*s.y*wr0';
	 yrr = wz0*s.y*wr2';
	 yrz = wz1*s.y*wr1';
	 yzz = wz2*s.y*wr0';
	 cnull = -[yrr yrz; yrz yzz]\[yr; yz];
	 ma = max(abs(cnull));
	 cnull = cnull/max(1,8*ma);
	 tr = tr+cnull(1);
	 tz = tz+cnull(2);
	 if ma < 1e-14
	   break
	 end
       end % End of k loop
       s.d = sqrt(min(abs(2*(s.v-y)/yrr),abs(2*(s.v-y)/yzz)));
       s.x = s.d < 1e-4; % All exits are valid since contour very close to x-point
    end
    if ~s.x
      % split sub cell
      s = split(s);
    end
  else
    % Walk through to new sub cell
    s = walk(s);
  end
  count = count+1;
end

% Find which of the 12 regions along the rectangle contains point tr, tz
function i = region(s,tr,tz)
for i = 1:12
  if (s.r(i)-tr)*(s.r(i+1)-tr)<=0 & (s.z(i)-tz)*(s.z(i+1)-tz)<=0
    break
  end
end

% Return coordinates for the point in region s.i where s.y=s.v
function [tr tz] = point(s)
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % for interpolation
i = s.i; % index to range containing traced point
f = s.f;
r = s.r;
z = s.z;
v = s.v;
y = s.y;
p = s.p;
ro = r(1); % Corners
ri = r(4);
zl = z(4);
zu = z(7);
r2 = r(i:i+1); % range containing traced point
z2 = z(i:i+1);
if any([1 2 3 7 8 9] == i)
  tz = z(i);
  tr = ((p(i+1)-v)*r(i)+(v-p(i))*r(i+1))/(p(i+1)-p(i));
  for j = 1:5
    wz0 = [1 tz tz^2 tz^3]*mx;
    wr0 = [1 tr tr^2 tr^3]*mx;
    wr1 = [0 1 2*tr 3*tr^2]*mx;
    wr2 = [0 0 2 6*tr]*mx;
    A = wz0*y*wr2';
    B = wz0*y*wr1';
    C = wz0*y*wr0'-v;
    D = -B/A;
    S2 = D^2-2*C/A;
    if abs(A) > 1e-4 & S2 >= 0
      tr = tr+D-sign(D)*sqrt(S2);
    else
      tr = tr-C/B;
    end
  end
  if (s.r(i)-tr)*(s.r(i+1)-tr) <= 0
    % The tr is still within range, it should be so all is well
  else
    % tr is either nan or outside range so revert to linear interpolation
    tr = ((p(i+1)-v)*r(i)+(v-p(i))*r(i+1))/(p(i+1)-p(i));
  end
else
  tr = s.r(s.i);
  tz = ((p(i+1)-v)*z(i)+(v-p(i))*z(i+1))/(p(i+1)-p(i));
  for j = 1:5
    wz0 = [1 tz tz^2 tz^3]*mx;
    wz1 = [0 1 2*tz 3*tz^2]*mx;
    wz2 = [0 0 2 6*tz]*mx;
    wr0 = [1 tr tr^2 tr^3]*mx;
    A = wz2*y*wr0';
    B = wz1*y*wr0';
    C = wz0*y*wr0'-v;
    D = -B/A;
    S2 = D^2-2*C/A;
    if abs(A) > 1e-4 & S2 >= 0
      tz = tz+D-sign(D)*sqrt(S2);
    else
      tz = tz-C/B;
    end
  end
  if (s.z(i)-tz)*(s.z(i+1)-tz) <= 0
    % The tz is still within range, it should be so all is well
  else
    % tz is either nan or outside range so revert to linear interpolation
    tz = ((p(i+1)-v)*z(i)+(v-p(i))*z(i+1))/(p(i+1)-p(i));
  end
end
