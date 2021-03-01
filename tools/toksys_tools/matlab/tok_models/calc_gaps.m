function [gapresp,gaps,gapspec] = calc_gaps(gapspec,tok_data_struct,resp,eq,idoplot)
%
%  USAGE: [gapresp,gaps,gapspec] = calc_gaps(gapspec,tok_data_struct,resp,eq,idoplot)
%
%  PURPOSE: Calculate gaps and responses thereof at specified locations
%           A gap is defined as the distance between the plasma boundary 
%           and the wall along a line through a specified "gap location", 
%           with direction of measurement defined by gradient of flux for a 
%           nominal equilibrium at that prescribed "gap location".
%           (This nominal equilibrium is also used to generate the 
%           linearized plasma response model.)
%
%  INPUTS: gapspec: gap specification on the form [r z gr gz]
%                   r,z is gap location and gr,gz nominal gradient
%                   If nominal gradient is 0,0 then it is taken from eq
%                   This creates a gapspec with eq as nominal equilibrium
%                   If gapspec is empty, a default inner, upper, outer, lower gap are set
%  tok_data_struct: toksys description of tokamak
%             resp: (optional) the output from a response model (gspert or rzrig)
%               eq: equilibrium structure
%          idoplot: flag to plot limiter, boundary, gap location, and gap vector
%
%  OUTPUTS: gapresp: response structure with fields:
%                    dgapdis: gap response to conductors [m/A]
%                    dgapdip,*dli,*dbetap, etc: exogenous responses
%              gaps: distances from boundary to wall through r,z along gr,gz [m]
%           gapspec: Complete gap specification on form [r z gr gz].

%
%  RESTRICTIONS: gapspec points must be within the grid

%
%  METHOD: 
%	
%  VERSION @(#)calc_gaps.m	1.6 02/26/15
%
%  WRITTEN BY:  Anders Welander  ON	5/18/11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4
  help calc_gaps
  disp('*********************************************************')
  disp('Supply at least 4 arguments! (gapspec, resp can be empty)')
  disp('*********************************************************')
  return
end

if ~exist('idoplot','var')
  idoplot = 0;
end

if exist('tok_data_struct','var') & isstruct(tok_data_struct)
  if size(tok_data_struct.limdata,2) == 2
    limdata = tok_data_struct.limdata';
  else
    limdata = tok_data_struct.limdata;
  end
  rlim = limdata(2,:);
  zlim = limdata(1,:);
  if min(rlim) < 0
    rlim = limdata(1,:);
    zlim = limdata(2,:);
  end
else
  rlim = eq.xlim;
  zlim = eq.ylim;
end
if rlim(1) ~= rlim(end) | zlim(1) ~= zlim(end)
  rlim(end+1) = rlim(1); zlim(end+1) = zlim(1);
end
nlim = length(rlim);

if isempty(gapspec) % Use default gapspec with inner, upper, outer, lower gaps
   k = find(zlim(1:end-1).*zlim(2:end) <= 0);
  [x2, j] = max(rlim(k));
  [x1, j] = min(rlim(k));
  xc = (x1+x2)/2; % Center of limiter
  r = interp1(zlim(k(j)+[0 1]),rlim(k(j)+[0 1]),0);
  gapspec(1,:) = [r 0 1 0];
  % Define default upper gap
  k = find((rlim(1:end-1)-xc).*(rlim(2:end)-xc) <= 0 & zlim(1:end-1)>0);
  [dum, j] = min(zlim(k));
  z = interp1(rlim(k(j)+[0 1]),zlim(k(j)+[0 1]),xc);
  gapspec(2,:) = [xc z 0 -1];
  % Define default outer gap
  k = find(zlim(1:end-1).*zlim(2:end) <= 0);
  [dum, j] = max(rlim(k));
  r = interp1(zlim(k(j)+[0 1]),rlim(k(j)+[0 1]),0);
  gapspec(3,:) = [r 0 -1 0];
  % Define default lower gap
  k = find((rlim(1:end-1)-xc).*(rlim(2:end)-xc) <= 0 & zlim(1:end-1)<0);
  [dum, j] = max(zlim(k));
  z = interp1(rlim(k(j)+[0 1]),zlim(k(j)+[0 1]),xc);
  gapspec(4,:) = [xc z 0 1];
end

ngap = size(gapspec,1);
if size(gapspec,2) < 4
  gapspec(1,4) = 0;
end 
gapresp = []; gaps = [];

nr = eq.nw;
nz = eq.nh;
dr = eq.dr;
dz = eq.dz;
rg = eq.rg;
zg = eq.zg;
rgg = ones(nz,1)*rg(:)';
zgg = zg(:)*ones(1,nr);
psizr = eq.psizr;
psibry = eq.psibry;
psimag = eq.psimag;
rbbbs = eq.rbbbs;
zbbbs = eq.zbbbs;
nbbbs = eq.nbbbs;
rmaxis = eq.rmaxis;
zmaxis = eq.zmaxis;

% Find weights in grid points to calculate value at a point using cubic Hermite spline (ref. wikipedia)
mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % Used for interpolation
ir4 = (-1:2)'*nz; % For finding flux on horizontal grid lines
iz4 = (-1:2)';    % For finding flux on vertical grid lines
neighbors = reshape(iz4*ones(1,4)+ones(4,1)*ir4',1,16);

% Determine gr,gz for any 0,0 entries
[y, yr, yz] = gs_interp2(rg, zg, psizr, gapspec(:,1), gapspec(:,2));
ifix = gapspec(:,3)==0 & gapspec(:,4)==0;
gapspec(ifix,3:4) = [yr(ifix) yz(ifix)];

% Determine the limiter points that the gap vectors intersect and also normalize gr,gz
for j = 1:ngap
  R = gapspec(j,1);
  Z = gapspec(j,2);
  gr = gapspec(j,3);
  gz = gapspec(j,4);
  g = norm([gr gz]);
  gr = gr/g;
  gz = gz/g;
  gapspec(j,3:4) = [gr gz]; % Now they are normalized in the specification
  d2min = inf; % Will be minimum distance^2 along gr,gz from r,z to limiter
  for k = 1:nlim-1
    m = [gr rlim(k+1)-rlim(k); gz zlim(k+1)-zlim(k)];
    if rank(m)==2
      kk = m\[R-rlim(k);Z-zlim(k)];
      if kk(2)>=0 & kk(2)<=1
	Rt = rlim(k)+kk(2)*(rlim(k+1)-rlim(k));
	Zt = zlim(k)+kk(2)*(zlim(k+1)-zlim(k));
	d2 = (R-Rt)^2+(Z-Zt)^2;
	if d2<d2min
	  rgl(j) = Rt;
	  zgl(j) = Zt;
	  d2min = d2;
	end
      end
    end
  end
end

% Make hires boundary
rbhires = interp1(1:nbbbs,rbbbs(1:nbbbs),1:.1:nbbbs);
zbhires = interp1(1:nbbbs,zbbbs(1:nbbbs),1:.1:nbbbs);

% Now calculate gaps for eq
for j = 1:ngap
  gr = gapspec(j,3);
  gz = gapspec(j,4);
  [d2, k] = min((rbhires-gapspec(j,1)).^2+(zbhires-gapspec(j,2)).^2);
  d = gr*(rbhires(k)-gapspec(j,1))+gz*(zbhires(k)-gapspec(j,2));
  R = gapspec(j,1)+d*gr;
  Z = gapspec(j,2)+d*gz;
  dstep = 1e-3;    
  psiprev = inf;
  while abs(dstep) > 1e-6
    R = R+gr*dstep;
    Z = Z+gz*dstep;
    kr0 = min(nr-3,max(1,floor((R-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
    kz1 = min(nz-2,max(2,ceil((Z-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
    k = kr0*nz+kz1;
    ii = k+neighbors; % ii indexes 16 points around R,Z
    tr = (R-rgg(k))/dr; tz = (Z-zgg(k))/dz;
    w = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
    psi = psizr(ii)*w';
    if (psi>psibry & dstep>0) | (psi<psibry & dstep<0) | (abs(psi-psibry)>abs(psiprev-psibry))
      dstep = -dstep/2;
    end
    psiprev = psi;
  end
  rs(j) = R;
  zs(j) = Z;
  gaps(j) = norm([rs(j)-rgl(j) zs(j)-zgl(j)]);
end
[y, yr, yz] = gs_interp2(rg, zg, (psizr-psimag)/(psibry-psimag), rs, zs);
for i = 1:ngap
  gaps(i) = (yr(i)*(rgl(i)-rs(i))+yz(i)*(zgl(i)-zs(i)))/norm([yr(i) yz(i)]);
end

if idoplot
  plot(gapspec(:,1),gapspec(:,2),'go','markerSize',8)
  hold on
  plot(rlim,zlim,'k','linew',2)
  plot(rbbbs(1:nbbbs),zbbbs(1:nbbbs),'b')
  plot(rgl,zgl,'rx','markerSize',12,'linew',2)
  for j = 1:ngap
    plot([rgl(j) rs(j)],[zgl(j) zs(j)],'m')
    gr = gapspec(j,3); gz = gapspec(j,4);
    a = angle(gr+i*gz);
    plot(rs(j)-[0 2*dr]*cos(a-pi/12),zs(j)-[0 2*dr]*sin(a-pi/12),'m')
    plot(rs(j)-[0 2*dr]*cos(a+pi/12),zs(j)-[0 2*dr]*sin(a+pi/12),'m')
    a = angle((rs(j)-rmaxis)+i*(zs(j)-zmaxis));
    rho = sqrt((rs(j)-rmaxis)^2+(zs(j)-zmaxis)^2);
    %text(rmaxis+0.8*rho*cos(a),zmaxis+0.8*rho*sin(a),num2str(j),'fonts',12,'fontw','bold','color','g','hori','center')
  end
  axis('image')
  legend('Gap locations');
  hold off
end

if isempty(resp)
  return
end

% Find the point that defines the boundary
% Begin by finding touch point candidate
rlhires = interp1(1:nlim,rlim,1:.1:nlim);
zlhires = interp1(1:nlim,zlim,1:.1:nlim);
psilim = interp2(rg,zg,psizr,rlhires,zlhires,'spline');
for j = 1:length(rlhires)
  gaphires(j) = min(sqrt((rlhires(j)-rbhires).^2+(zlhires(j)-zbhires).^2));
end
[mingap, ilim] = min(gaphires);
if ilim == 1 | ilim==nlim
  k = [nlim-1 1 2];
else
  k = ilim+[-1 0 1];
end
d1 = min(-sqrt((rlhires(k(1))-rlhires(k(2)))^2+(zlhires(k(1))-zlhires(k(2)))^2),-1e-19);
d2 = max(+sqrt((rlhires(k(3))-rlhires(k(2)))^2+(zlhires(k(3))-zlhires(k(2)))^2),+1e-19);
if (psilim(k(1))-psilim(k(2)))*(psilim(k(3))-psilim(k(2)))>0
  warning off; p = polyfit([d1 0 d2],psilim(k)-psilim(k(2)),2); d = -p(2)/2/p(1); warning on
else
  d = 0;
end
rlim = spline([d1 0 d2],rlhires(k),d); zlim = spline([d1 0 d2],zlhires(k),d);
% Now find x-point candidate
for j = 1:nbbbs-1
  k = 1+floor((rbbbs(j)-rg(1))/dr)*nz+floor((zbbbs(j)-zg(1))/dz);
  ii = k+neighbors; % ii indexes 16 points around R,Z
  tr = (rbbbs(j)-rgg(k))/dr; tz = (zbbbs(j)-zgg(k))/dz;
  wr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
  wz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  bzbbbs(j) = +psizr(ii)*wz'/rbbbs(j);
  brbbbs(j) = -psizr(ii)*wr'/rbbbs(j);
end
bpbbbs = sqrt(brbbbs.^2+bzbbbs.^2);
[bpmin, ix] = min(bpbbbs);
rx = rbbbs(ix); zx = zbbbs(ix);
kr0 = min(nr-3,max(1,floor((rx-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
kz1 = min(nz-2,max(2,ceil((zx-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
k = kr0*nz+kz1;
iix = k+neighbors; % iix indexes 16 points around x point
pp = psizr(iix); j = j-1;
tr = (rx-rgg(k))/dr; tz = (zx-zgg(k))/dz;
wx = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
wxr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
wxz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
wxrr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
wxzz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
wxrz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
xshift_Bx = -inv([pp*wxrr' pp*wxrz';pp*wxrz' pp*wxzz']);
cx = xshift_Bx*[pp*wxr'; pp*wxz']; rx = rx+cx(1); zx = zx+cx(2);
% If the candidate x-point is further from a null than min(gaphires) then limited plasma
if norm(cx)>min(gaphires)
  ilimited = 1;
  r0 = rlim;
  z0 = zlim;
else
  r0 = rbbbs(ix);
  z0 = zbbbs(ix);
end
kr0 = min(nr-3,max(1,floor((r0-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
kz1 = min(nz-2,max(2,ceil((z0-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
k = kr0*nz+kz1;
ii0 = k+neighbors; % ii0 indexes 16 points around limiting point
tr = (r0-rgg(k))/dr; tz = (z0-zgg(k))/dz;
w0 = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
m0 = zeros(16,nz*nr); % will hold mutuals to the 16 points around r0,z0
for jj = 1:16
  z_idx = mod(ii0(jj)-1,nz)+1; % z index of output element ii(jj)
  r_idx = ceil(ii0(jj)/nz);    % r index of output element ii(jj)
  for jr=1:nr
    k1 = (r_idx-1)*nz;
    k2 = (jr-1)*nz;
    for jz=1:nz
      m0(jj,k2+jz) = mpp(k1+abs(z_idx-jz)+1,jr);
    end
  end
end

% Parse the resp structure
if isfield(resp,'dcphidis') % true for gspert
elseif isfield(resp,'drdis') & isfield(resp,'dcdr') % true for rzrig
  resp.dcphidis = resp.dcdr(:)*resp.drdis+resp.dcdz(:)*resp.dzdis;
  resp.units.dcphidis = 'A/A';
  resp.dcphidip = resp.dcdr(:)*resp.drdip+resp.dcdz(:)*resp.dzdip;
  resp.units.dcphidip = 'A/A';
  resp.dcphidli = resp.dcdr(:)*resp.drdli;
  resp.units.dcphidli = 'A/1';
  resp.dcphidbetap = resp.dcdr(:)*resp.drdbetap+resp.dcdz(:)*resp.dzdbetap;
  resp.units.dcphidbetap = 'A/1';
else
  wait('error calc_gap: can not interpret structure resp')
  return
end
% Since exogenous variables may vary, use field names to extract all that apply
f = fieldnames(resp); % f holds names of fields in resp
ind_resp_in_f = []; % will become indices in f containing dcphid*
for j = 1:length(f)
  if strfind(char(f(j)),'dcphid') == 1
    ind_resp_in_f(end+1) = j;
  end
end

% Now calculate responses using dcphid* for plasma response and tok_data_struct for mutuals
for j = 1:ngap
  gr = gapspec(j,3); gz = gapspec(j,4);
  s = sign(gr*(rs(j)-rgl(j))+gz*(zs(j)-zgl(j)));
  % Get mutuals from grid to the 16 points that surround this point rs,zs on separatrix.
  kr0 = min(nr-3,max(1,floor((rs(j)-rg(1))/dr))); % r index 0-start, allowed values: 1:nr-3
  kz1 = min(nz-2,max(2,ceil((zs(j)-zg(1))/dz))); % z index 1-start, allowed values: 2:nz-2
  k = kr0*nz+kz1;
  ii = k+neighbors; % ii indexes 16 points around gap location
  m = zeros(16,nz*nr); % will hold mutuals to the 16 points around rs,zs on separatrix
  for jj = 1:16
    z_idx = mod(ii(jj)-1,nz)+1; % z index of output element ii(jj)
    r_idx = ceil(ii(jj)/nz);    % r index of output element ii(jj)
    for jr=1:nr
      k1 = (r_idx-1)*nz;
      k2 = (jr-1)*nz;
      for jz=1:nz
	m(jj,k2+jz) = mpp(k1+abs(z_idx-jz)+1,jr);
      end
    end
  end
  % We now have dpsizr(ii)/dis = m*resp.dcphidis+[mpc(ii,:) mpv(ii,:)], etc
  tr = (rs(j)-rgg(k))/dr; tz = (zs(j)-zgg(k))/dz;
  w = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  wr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
  wz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
  dpsidgv = psizr(ii)*(gr*wr+gz*wz)'; % d(psi)/d(gapvector)
  % We now have dpsi(rs,zs)/dis = w*(m*dpsizr(ii))/dis, etc
  for jj = 1:length(ind_resp_in_f)
    respobj = char(f(ind_resp_in_f(jj)));
    if strmatch(respobj,'dcphidis')
      gapresp.dgapdis(j,:) = -((w*m-w0*m0)*resp.dcphidis + w*[mpc(ii,:) mpv(ii,:)]-w0*[mpc(ii0,:) mpv(ii0,:)])/dpsidgv;
    elseif eval(['prod(size(resp.' respobj '))']) == nr*nz % respobj may be on form (nz,nr) so collapse
      eval(['gapresp.dgapd' respobj(7:end) '(j,:) = -(w*m-w0*m0)*resp.' respobj '(:)/dpsidgv;'])
    else
      eval(['gapresp.dgapd' respobj(7:end) '(j,:) = -(w*m-w0*m0)*resp.' respobj '/dpsidgv;'])
    end
  end
  gapresp.dgapdpsi(j,1) = -1/dpsidgv;
end
for jj = 1:length(ind_resp_in_f)
  respobj = char(f(ind_resp_in_f(jj)));
  eval(['gapresp.descriptions.dgapd' respobj(7:end) ' = ''gap response to ' respobj(7:end) ''';'])
  s = eval(['resp.units.dcphid' respobj(7:end)]);
  eval(['gapresp.units.dgapd' respobj(7:end) ' = [''m/A * '' s]; '])
end
gapresp.descriptions.dgapdpsi = 'dgap/dpsi, where dpsi = (dpsi at gap) - (dpsi at boundary defining point)';
gapresp.units.dgapdpsi = 'm/Wb';
%  gapresp.descriptions.dpsidis = 'flux response to is (conductor currents) [Wb/A]';
%  gapresp.descriptions.dbrdis = 'Br response to is (conductor currents) [T/A]';
%  gapresp.descriptions.dbzdis = 'Bz response to is (conductor currents) [T/A]';


