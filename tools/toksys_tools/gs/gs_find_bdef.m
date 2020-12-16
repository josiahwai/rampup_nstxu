%  USAGE:   gs_find_bdef
%
%  PURPOSE: Find boundary-defining point (x or touch point)
%
%  INPUTS: psizr, the flux at nz vertical positions * nr radial positions
%          ilimgg, flags for grid, -1 = outside vessel, >-1 = inside
%          rgg, zgg, dr, dz, nr, nz (grid variables)
%          For cubic interpolation on the grid:
%            mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2
%            neighbors = reshape((-1:2)'*ones(1,4)+ones(4,1)*(-1:2)*nz,1,16)
%
%  OUTPUTS: rbdef, zbdef, position that defines the boundary (touch or x point)
%           rx1, zx1, psix1, position and flux of most important x-point below axis
%           rx2, zx2, psix2, position and flux of most important x-point below axis
%           rtl, ztl, psitl, position and flux of most important limiter-point
%           wb, iib, weights and indices such that psibry = wb*psizr(iib)
%           drbefdpsi, dzbdefdpsi, weights such that drbdef = drbdefdpsi*dpsizr(iib)
%           Weights and indices for x1, x2, tl are called:
%             wx1, iix1, drx1dpsi, dzx1dpsi, wx2, iix2, drx2dpsi, dzx2dpsi, wtl, iitl
%           ix1, flag for lower x-point defines the boundary
%           ix2, flag for upper x-point defines the boundary
%           itl, flag for plasma touches limiter
%           zbot, ztop, limits for z position of plasma
%           zbotg, ztopg, if zgg(k)==zbotg+dz/2 then zbot is in range zgg(k)+[0 dz]
%	
%  METHOD: 
	
%  NOTES:  
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	4/12/13
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for LOWER NULL
if ~(x1exists & x1inside)
  RBmin_per_row = min((psizr(2:nz,1:nr-1)-psizr(1:nz-1,1:nr-1))'.^2+(psizr(1:nz-1,2:nr)-psizr(1:nz-1,1:nr-1))'.^2);
  [dum, k] = min(RBmin_per_row(zg(1:nz-1) < zmaxis-dz));
  zx1 = zg(k) + dz/2;
  [dum, k] = min((psizr(k+1,1:nr-1)-psizr(k,1:nr-1)).^2+(psizr(k,2:nr)-psizr(k,1:nr-1)).^2);
  rx1 = rg(k) + dr/2;
end

twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
j = 9;
while j > 0 & twopirbrzmax > 1e-10 % Try zooming in on x-point with Newton Rhapson.
  j = j-1;
  % Find indices and weights for grid points around the x-point
  kr0 = min(nr-3,max(1,floor((rx1-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
  kz1 = min(nz-2,max(2,ceil((zx1-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
  k = kr0*nz+kz1;
  iix1 = k+neighbors'; % indexes 16 points around x point
  pp = psizr(iix1);
  tr = (rx1-rgg(k))/dr;
  tz = (zx1-zgg(k))/dz;
  wx1 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  wx1r = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
  wx1z = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  wx1rr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
  wx1zz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  wx1rz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
  psix1_r = wx1r*pp;
  psix1_z = wx1z*pp;
  psix1_rr = wx1rr*pp;
  psix1_rz = wx1rz*pp;
  psix1_zz = wx1zz*pp;
  x1shift_Bx = -inv([psix1_rr psix1_rz; psix1_rz psix1_zz]);
  cx1 = x1shift_Bx*[psix1_r; psix1_z];
  rx1 = rx1+cx1(1);
  zx1 = zx1+cx1(2);
  twopirbrzmax = max(abs([psix1_r psix1_z]));
end
psix1 = wx1*pp;
drx1dpsi = x1shift_Bx(1,1)*wx1r+x1shift_Bx(1,2)*wx1z;
dzx1dpsi = x1shift_Bx(2,1)*wx1r+x1shift_Bx(2,2)*wx1z;

if twopirbrzmax > 1e-10 || zx1 > zmaxis-dz
  x1exists = false; % No x-point found
else
  x1exists = true; % x-point found
end
if x1exists % X-point found. Now check if it is inside the limiter
  if ilimgg(k) == 0 % In this case it is clearly inside
    x1inside = true; % Flag that x1 is a real (and important) x-point	  
  elseif ilimgg(k) == -1 % That means clearly outside
    x1inside = false;
  else
    x1inside = true;
    for j = max(1,ilimgg(k)-5):min(nlim-1,ilimgg(k)+5);
      mlinex = [Rlim(j+1)-Rlim(j) rmaxis-rx1; Zlim(j+1)-Zlim(j) zmaxis-zx1];
      if rcond(mlinex) > 0
	kg = inv(mlinex)*[rmaxis-Rlim(j); zmaxis-Zlim(j)];
	if kg(1)>=0 & kg(1)<=1 & kg(2)>0 & kg(2)<1 % x1 outside limiter
	  x1inside = false;
	end
      end
    end
  end
end

if x1exists & x1inside
  zbot = zx1; % Min z for where plasma can be
else
  psix1 = inf/(psibry-psimag);
  zbot = min(zbbbs(1:nbbbs))-dzbbbs_max;
end
% If zgg(k) > zbotg then zgg(k) *may* be a lower inner corner of cell with point of z > zbot
zbotg = (floor((zbot-zg(1))/dz)+0.5)*dz+zg(1);
psibarx1 = (psix1-psimag)/(psibry-psimag);


% Check for UPPER NULL
if ~(x2exists & x2inside)
  RBmin_per_row = min((psizr(2:nz,1:nr-1)-psizr(1:nz-1,1:nr-1))'.^2+(psizr(1:nz-1,2:nr)-psizr(1:nz-1,1:nr-1))'.^2);
  inz1 = zg(1:nz-1) > zmaxis+dz;
  [dum, k] = min(RBmin_per_row(inz1));	
  k = k + sum(~inz1);
  zx2 = zg(k) + dz/2;
  [dum, k] = min((psizr(k+1,1:nr-1)-psizr(k,1:nr-1)).^2+(psizr(k,2:nr)-psizr(k,1:nr-1)).^2);
  rx2 = rg(k) + dr/2;
end

twopirbrzmax = inf; % Maximum of 2*pi*R*Br, 2*pi*R*Bz at calculated null position
j = 9;
while j > 0 & twopirbrzmax > 1e-10 % Try zooming in on x-point with Newton Rhapson.
  j = j-1;
  % Find indices and weights for grid points around the x-point
  kr0 = min(nr-3,max(1,floor((rx2-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
  kz1 = min(nz-2,max(2,ceil((zx2-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
  k = kr0*nz+kz1;
  iix2 = k+neighbors'; % indexes 16 points around x point
  pp = psizr(iix2);
  tr = (rx2-rgg(k))/dr;
  tz = (zx2-zgg(k))/dz;
  wx2 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  wx2r = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
  wx2z = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  wx2rr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
  wx2zz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  wx2rz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
  psix2_r = wx2r*pp;
  psix2_z = wx2z*pp;
  psix2_rr = wx2rr*pp;
  psix2_rz = wx2rz*pp;
  psix2_zz = wx2zz*pp;
  x2shift_Bx = -inv([psix2_rr psix2_rz; psix2_rz psix2_zz]);
  cx2 = x2shift_Bx*[psix2_r; psix2_z];
  rx2 = rx2+cx2(1);
  zx2 = zx2+cx2(2);
  twopirbrzmax = max(abs([psix2_r psix2_z]));
end
psix2 = wx2*pp;
drx2dpsi = x2shift_Bx(1,1)*wx2r+x2shift_Bx(1,2)*wx2z;
dzx2dpsi = x2shift_Bx(2,1)*wx2r+x2shift_Bx(2,2)*wx2z;

if twopirbrzmax > 1e-10 || zx2 < zmaxis+dz
  x2exists = false; % No x-point found
else
  x2exists = true; % x-point found
end
if x2exists % X-point found. Now check if it is inside the limiter
  if ilimgg(k) == 0 % In this case it is clearly inside
    x2inside = true; % Flag that x2 is a real (and important) x-point
  elseif ilimgg(k) == -1 % That means clearly outside
    x2inside = false;
  else
    x2inside = true;
    for j = max(1,ilimgg(k)-5):min(nlim-1,ilimgg(k)+5);
      mlinex = [Rlim(j+1)-Rlim(j) rmaxis-rx2; Zlim(j+1)-Zlim(j) zmaxis-zx2];
      if rcond(mlinex) > 0
	kg = inv(mlinex)*[rmaxis-Rlim(j); zmaxis-Zlim(j)];
	if kg(1)>=0 & kg(1)<=1 & kg(2)>0 & kg(2)<1 % x2 outside limiter
	  x2inside = false;
	end
      end
    end
  end
end

if x2exists & x2inside
  ztop = zx2; % Max z for where plasma can be
else
  psix2 = inf/(psibry-psimag);
  ztop = max(zbbbs(1:nbbbs))+dzbbbs_max;
end
% If zgg(k) < ztopg then zgg(k) *may* be a lower inner corner of cell with point of z < ztop
ztopg = (floor((ztop-zg(1))/dz)+0.5)*dz+zg(1);
psibarx2 = (psix2-psimag)/(psibry-psimag);


% Find candidate for touch point
psilim = sum(wl.*psizr(iil),2); % Flux at limiter
il = rl>rg(1) & rl<rg(nr) & zl>zbot & zl<ztop; % Flag relevant limiter points
% Find point with lowest psibar among relevant limiter points
dum = inf;
for j = 1:nl
  if il(j) & psilim(j)/(psibry-psimag) < dum
    dum = psilim(j)/(psibry-psimag);
    ktl = j;
  end
end
% Indices to adjacent points
if ktl > 1
  ktlm = ktl-1;
else
  ktlm = nl-1;
end
if ktl < nl
  ktlp = ktl+1;
else
  ktlp = 2;
end
iitl = iil(ktl,:)';
% Is the touch point candidate between ktl and a lower index?
vlim = [rl(ktlm)-rl(ktl) zl(ktlm)-zl(ktl)];
psi1para = vlim*([wlr(ktl,:); wlz(ktl,:)]*psizr(iitl));
psi2para = vlim*([wlr(ktlm,:); wlz(ktlm,:)]*psizr(iil(ktlm,:))');
if psi1para*psi2para > 0 
  % Is the touch point candidate between ktl and a higher index?
  vlim = [rl(ktlp)-rl(ktl) zl(ktlp)-zl(ktl)];
  psi1para = vlim*([wlr(ktl,:); wlz(ktl,:)]*psizr(iitl));
  psi2para = vlim*([wlr(ktlp,:); wlz(ktlp,:)]*psizr(iil(ktlp,:))');
end

if psi1para*psi2para < 0 % Field parallel with limiter between 1 and 2
  
  rtl = rl(ktl);
  ztl = zl(ktl);
  ulim = vlim/norm(vlim);
  j = 9;
  psild_min = realmax;
  psild = realmax;
  while abs(psild) > 1e-9 & j > 0
    kr0 = min(nr-3,max(1,floor((rtl-rg(1))/dr)));
    kz1 = min(nz-3,max(1,floor((ztl-zg(1))/dz)))+1;
    k = kr0*nz+kz1;
    iitl = k+neighbors';
    tr = (rtl-rgg(k))/dr;
    tz = (ztl-zgg(k))/dz;

    wr0 = [1 tr tr^2 tr^3]*mx;
    wz0 = mx'*[1 tz tz^2 tz^3]';
    wr1 = [0 1 2*tr 3*tr^2]*mx/dr*ulim(1);
    wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz*ulim(2);
    wr2 = [0 0 2 6*tr]*mx/dr^2*ulim(1)^2;
    wz2 = mx'*[0 0 2 6*tz]'/dz^2*ulim(2)^2;

    wtld = reshape(wz0*wr1 + wz1*wr0, 1, 16);
    wtlb = reshape(wz0*wr2 + 2*wz1*wr1 + wz2*wr0, 1, 16);

    psild = wtld*psizr(iitl);
    psilb = wtlb*psizr(iitl);
    
    % Do this to get an okay value if these iterations should diverge
    if abs(psild) < psild_min
      rtl_best = rtl;
      ztl_best = ztl;
      psild_min = abs(psild);
    end

    dlim = -psild/psilb;

    rtl = rtl + dlim*ulim(1);
    ztl = ztl + dlim*ulim(2);
    
    j = j-1;
  end
  rtl = rtl_best;
  ztl = ztl_best;
   
else % Field not parallel to limiter => this is corner sticking into plasma
  rtl = rl(ktl);
  ztl = zl(ktl);
  ulim = [0 0];
end

tr = (rtl-rgg(iitl(6)))/dr;
tz = (ztl-zgg(iitl(6)))/dz;
wr0 = [1 tr tr^2 tr^3]*mx;
wz0 = mx'*[1 tz tz^2 tz^3]';
wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
wr2 = [0 0 2 6*tr]*mx/dr^2;
wz2 = mx'*[0 0 2 6*tz]'/dz^2;
wtl = reshape(wz0*wr0,1,16);
wtlr = reshape(wz0*wr1,1,16);
wtlz = reshape(wz1*wr0,1,16);
wtlrr = reshape(wz0*wr2,1,16);
wtlrz = reshape(wz1*wr1,1,16);
wtlzz = reshape(wz2*wr0,1,16);
wtld = wtlr*ulim(1)+wtlz*ulim(2);
wtlb = reshape(wz0*wr2*ulim(1)^2+2*wz1*wr1*ulim(1)*ulim(2)+wz2*wr0*ulim(2)^2,1,16);
psitlbis = wtlb*psizr(iitl);
if abs(psitlbis) < 1e-9
  psitlbis = 1e-9;
end

psitl = wtl*psizr(iitl);
psibartl = (psitl-psimag)/(psibry-psimag);

% Which point defines the boundary?
if psibartl < min(psibarx1,psibarx2)
  itl = true;
  ix1 = false;
  ix2 = false;
  rbdef = rtl;
  zbdef = ztl;
  wb = wtl;
  wbr = wtlr;
  wbz = wtlz;
  wbrr = wtlrr;
  wbrz = wtlrz;
  wbzz = wtlzz;
  iib = iitl;
  drbdefdpsi = -ulim(1)/psitlbis*wtld;
  dzbdefdpsi = -ulim(2)/psitlbis*wtld;
elseif psibarx1 < psibarx2
  itl = false;
  ix1 = true;
  ix2 = false;
  rbdef = rx1;
  zbdef = zx1;
  wb = wx1;
  wbr = wx1r;
  wbz = wx1r;
  wbrr = wx1rr;
  wbrz = wx1rz;
  wbzz = wx1zz;
  iib = iix1;
  drbdefdpsi = drx1dpsi;
  dzbdefdpsi = dzx1dpsi;
else
  itl = false;
  ix1 = false;
  ix2 = true;
  rbdef = rx2;
  zbdef = zx2;
  wb = wx2;
  wbr = wx2r;
  wbz = wx2r;
  wbrr = wx2rr;
  wbrz = wx2rz;
  wbzz = wx2zz;
  iib = iix2;
  drbdefdpsi = drx2dpsi;
  dzbdefdpsi = dzx2dpsi;
end
