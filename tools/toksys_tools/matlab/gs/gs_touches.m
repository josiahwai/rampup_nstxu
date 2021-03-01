%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_touches
%
%  PURPOSE: Find points where plasma may touch xor initiate
%
%  INPUTS:  psizr, flux on grid
%
%  OUTPUTS:  touches, information about possible touch points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander ON 2015-02-22
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 0; % Number of possible touch points
psilim =  sum(wl.*psizr(iil),2); % Flux at limiter
psidlim = sum(wld.*psizr(iil),2);
psiblim = sum(wlb.*psizr(iil),2);
psitlim = sum(wlt.*psizr(iil),2);
limsq = psiblim.^2-2*psitlim.*psidlim;
sqrtlimsq = sqrt(limsq);
xl1 = (-psiblim-sqrtlimsq)./psitlim; % 2 Solutions for local flux extremum
xl2 = (-psiblim+sqrtlimsq)./psitlim;
kl1 = limsq >= 0 & xl1 >= 0 & xl1 <= dl; % True if extremum between lim points
kl2 = limsq >= 0 & xl2 >= 0 & xl2 <= dl;
kl3 = diff(psilim).*diff(psilim([nl-1 1:nl-1])) < 0;
for i = 1:nl-1
  if psidlim(i)*psidlim(i+1) < 0 & ~kl1(i) & ~kl2(i)
    % A point of psidlim = 0 is indeed in this range
    xl1(i) = psidlim(i)/(psidlim(i)-psidlim(i+1))*dl(i);
    kl1(i) = 1;
  end
  if kl1(i) | kl2(i)
    iitl = iil(i,:)';
    if kl1(i)
      flim = xl1(i)/dl(i);
    else
      flim = xl2(i)/dl(i);
    end
    rtl = rl(i)+flim*drl(i);
    ztl = zl(i)+flim*dzl(i);
    k = iitl(6);
    pp = psizr(iitl);
    tr = (rtl-rgg(k))/dr;
    tz = (ztl-zgg(k))/dz;
    m = 0;
    ur = drl(i)/dl(i);
    uz = dzl(i)/dl(i);
    dlim = dl(i);
    while abs(dlim) > 1e-9*dl(i) & m < 9
      m = m+1;
      wr0 = [1 tr tr^2 tr^3]*mx;
      wz0 = mx'*[1 tz tz^2 tz^3]';
      wr1 = [0 1 2*tr 3*tr^2]*mx/dr;
      wz1 = mx'*[0 1 2*tz 3*tz^2]'/dz;
      wr2 = [0 0 2 6*tr]*mx/dr^2;
      wz2 = mx'*[0 0 2 6*tz]'/dz^2;
      wr3 = [0 0 0 6]*mx/dr^3;
      wz3 = mx'*[0 0 0 6]'/dz^3;
      wtld = reshape(wz0*wr1*ur + wz1*wr0*uz, 1, 16);
      wtldin = reshape(wz0*wr1*turnin(1,2)*uz + wz1*wr0*turnin(2,1)*ur, 1, 16);
      wtlb = reshape(wz0*wr2*ur^2 + 2*wz1*wr1*ur*uz + wz2*wr0*uz^2, 1, 16);
      wtlt = reshape(wz0*wr3*ur^3 + 3*wz1*wr2*ur^2*uz + ...
                                    3*wz2*wr1*ur*uz^2 + wz3*wr0*uz^3, 1, 16);  
      psild = wtld*pp;
      psildin = wtldin*pp;
      psilb = wtlb*pp;
      psilt = wtlt*pp;
      if abs(psilt) > 1e-9
	dum = psilb^2-2*-psilt*psild;
	sqrtdum = sqrt(dum);
	dlim1 = (psilb-sqrtdum)/psilt;
	dlim2 = (psilb+sqrtdum)/psilt;
	if dum >= 0
	  if abs(dlim1) < abs(dlim2)
	    dlim = dlim1;
	  else
	    dlim = dlim2;
	  end
	else
	  dlim = 0;
	end
      elseif psilb ~= 0
        dlim = -psild/psilb;
      end
      tr = tr + dlim*ur/dr;
      tz = tz + dlim*uz/dz;
      flim = flim + dlim/dl(i);
    end
    rtl = rgg(k) + tr*dr;
    ztl = zgg(k) + tz*dz;
    wtl = reshape(wz0*wr0,1,16);
    psitl = wtl*psizr(iitl);
    drdpsi = -ur/psilb*wtld;
    dzdpsi = -uz/psilb*wtld;
    if abs(flim) < 1.1
      n = n+1;
      touches.count = n;
      touches.r(n) = rtl;
      touches.z(n) = ztl;
      touches.psi(n) = psitl;
      touches.w(n,:) = wtl;
      touches.ii(n,:) = iitl';
      touches.drdpsi(n,:) = drdpsi;
      touches.dzdpsi(n,:) = dzdpsi;
      touches.limiter(n) = i;
      touches.dpsidinward(n) = psildin;
      touches.bdefpossible(n) = psildin*psilb < 0;
      touches.inicursign(n) = -sign(psildin)*(psildin*psilb >= 0);
      touches.placursign(n) = sign(psildin)*(psildin*psilb < 0);
    end
  elseif kl3(i) & concavel(i)
    n = n+1;
    touches.count = n;
    touches.psi(n) = psilim(i);
    touches.r(n) = rl(i);
    touches.z(n) = zl(i);
    touches.w(n,:) = wl(i,:);
    touches.ii(n,:) = iil(i,:);
    touches.drdpsi(n,1:16) = 0;
    touches.dzdpsi(n,1:16) = 0;
    touches.limiter(n) = i;
    touches.dpsidinward(n) = 1;
    touches.bdefpossible(n) = logical(1);
    touches.inicursign(n) = 1;
    touches.placursign(n) = 1;
  end
end
ds = [];
ds.count = 'Number of points';
ds.r = 'Radius of possible touch or breakdown point';
ds.z = 'Height of possible touch or breakdown point';
ds.psi = 'The flux at the points';
ds.w = 'Weights for calculating psi = w*psizr(ii)';
ds.ii = 'Indices for calculating psi = w*psizr(ii)';
ds.drdpsi = 'How point moves radially in response to dpsizr(ii)';
ds.dzdpsi = 'How point moves vertically in response to dpsizr(ii)';
ds.limiter = 'index (i) into rl, zl, point is in range rl(i:i+1),zl(i:i+1)';
ds.dpsidinward = 'Flux gradient toward the inside of the machine';
ds.bdefpossible = 'If true plasma may touch, if false plasma may initiate';
ds.inicursign = 'Sign of total plasma current for initiating plasma';
ds.placursign = 'Sign of total plasma current for existing plasma';
touches.descriptions = ds;

return

if exist('plotit','var') & plotit
  clf
  hold on
  plot(rg(1)+ones(2,1)*linspace(-0.5,nr-0.5,nr+1)*dr,...
     [zg(1)-dz/2; zg(nz)+dz/2]*ones(1,nr+1),'--','color',[.9 .9 .9]);
  plot([rg(1)-dr/2; rg(nr)+dr/2]*ones(1,nz+1),...
     zg(1)+ones(2,1)*linspace(-0.5,nz-0.5,nz+1)*dz,'--','color',[.9 .9 .9]);
  for i = 1:nz
    text(rg(nr)+dr/2,zg(i),num2str(i),'vert','mid','color',[.8 .8 .8])
  end
  for j = 1:nr
    text(rg(j),zg(nz)+dz/2,num2str(j),'hori','cen','color',[.8 .8 .8])
  end
  plot(rl,zl,'k','linew',5)
  contour(rg,zg,psizr,touches.psi)
  k = touches.placursign > 0;
  plot(touches.r(k),touches.z(k),'rx','markers',24,'linew',3)
  k = touches.placursign < 0;
  plot(touches.r(k),touches.z(k),'bx','markers',24,'linew',3)
  k = touches.inicursign > 0;
  plot(touches.r(k),touches.z(k),'ro','markers',18,'linew',3)
  k = touches.inicursign < 0;
  plot(touches.r(k),touches.z(k),'bo','markers',18,'linew',3)
  drl = turnin*[diff(rl); diff(zl)];
  plot([1;1]*(rl(1:nl-1)+rl(2:nl))/2+[0;1]*drl(1,:), ...
       [1;1]*(zl(1:nl-1)+zl(2:nl))/2+[0;1]*drl(2,:))
  axis image
  zoom on
  %plot(rl(concavel),zl(concavel),'bo','linew',5)
  %plot(rl(convexl),zl(convexl),'ro','linew',5)
end
