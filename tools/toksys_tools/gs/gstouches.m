function touches = gstouches(c,psizr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   touches = gstouches(c,psizr);
%
%  PURPOSE: Find limiter points where plasma may touch xor initiate
%
%  INPUTS:  c, configuration data made by gsconfigure
%           psizr, flux on the grid c.rg, c.zg
%
%  OUTPUTS: touches, information about possible touch points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
%
%  WRITTEN BY:  Anders Welander ON 2016-11-08
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2; % Used for interpolation

rl = c.rl;
zl = c.zl;
nl = c.nl;
il = c.il;
wl = c.wl;
wld1 = c.wld1;
wld2 = c.wld2;
drl = c.drl;
dzl = c.dzl;
dl = c.dl;
irl = c.irl;
izl = c.izl;
trl = c.trl;
tzl = c.tzl;
kl = c.kl;
turnin = c.turnin;
concavel = c.concavel;

touches.count = 0;
touches.psi = nan(1,nl-1);
touches.r = nan(1,nl-1);
touches.z = nan(1,nl-1);
touches.w = nan(nl-1,16);
touches.ii = ones(nl-1,16);
touches.drdpsi = nan(nl-1,16);
touches.dzdpsi = nan(nl-1,16);
touches.limiter = nan(1,nl-1);
touches.dpsidinward = false(1,nl-1);
touches.bdefpossible = false(1,nl-1);
touches.inicursign = zeros(1,nl-1);
touches.placursign = zeros(1,nl-1);
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

yl1 = sum(wl.*psizr(il),2);
yd1 = sum(wld1.*psizr(il),2);
yd2 = sum(wld2.*psizr(il),2);

% Sign of derivative along limiter changes at possible touch points
kl = yd1.*yd2 < 0;

% Corners that may touch plasma
kc = yd1.*yd1([2:end 1]) < 0 & ~kl;

for i = 1:nl-1
  if kl(i)
    touches.count = touches.count+1;
    touches.limiter(touches.count) = i;
  end
end

% Process possible touch points on limiter surface
for i = 1:touches.count
  j = touches.limiter(i);
  x1 = 0;
  x2 = 1;
  ur = drl(j)/dl(j);
  uz = dzl(j)/dl(j);
  for k = 1:9
    x = x1+(x2-x1)*yd1(j)/(yd1(j)-yd2(j));
    tr = trl(j)+x*drl(j);
    tz = tzl(j)+x*dzl(j);
    wr0 = [1 tr tr^2 tr^3]*mx;
    wz0 = [1 tz tz^2 tz^3]*mx;
    wr1 = [0 1 2*tr 3*tr^2]*mx;
    wz1 = [0 1 2*tz 3*tz^2]*mx;
    yr = wz0*psizr(izl(j)-1:izl(j)+2,irl(j)-1:irl(j)+2)*wr1';
    yz = wz1*psizr(izl(j)-1:izl(j)+2,irl(j)-1:irl(j)+2)*wr0';
    yd = ur*yr+uz*yz;
    if yd*yd1(j) > 0
      x1 = x;
      yd1(j) = yd;
    elseif yd*yd2(j) > 0
      x2 = x;
      yd2(j) = yd;
    else
      break % Found exact point, break to avoid risk of dividing by 0
    end
  end
  ydin = turnin(1,2)*uz*yr+turnin(2,1)*ur*yz;
  wr2 = [0 0 2 6*tr]*mx;
  wz2 = [0 0 2 6*tz]*mx;
  psi = wz0*psizr(izl(j)-1:izl(j)+2,irl(j)-1:irl(j)+2)*wr0';
  yrr = wz0*psizr(izl(j)-1:izl(j)+2,irl(j)-1:irl(j)+2)*wr2';
  yrz = wz1*psizr(izl(j)-1:izl(j)+2,irl(j)-1:irl(j)+2)*wr1';
  yzz = wz2*psizr(izl(j)-1:izl(j)+2,irl(j)-1:irl(j)+2)*wr0';
  ylb = ur^2*yrr+2*ur*uz*yrz+uz^2*yzz;
  wtl = reshape(wz0'*wr0,1,16);
  wtld = reshape(wz0'*wr1*ur+wz1'*wr0*uz,1,16);
  drdpsi = -ur/ylb*wtld;
  dzdpsi = -uz/ylb*wtld;
  touches.r(i) = irl(j)+tr;
  touches.z(i) = izl(j)+tz;
  touches.psi(i) = psi;
  touches.w(i,:) = wtl;
  touches.ii(i,:) = il(j,:);
  touches.drdpsi(i,:) = drdpsi;
  touches.dzdpsi(i,:) = dzdpsi;
  touches.limiter(i) = j;
  touches.dpsidinward(i) = ydin;
  touches.bdefpossible(i) = ydin*ylb < 0;
  touches.inicursign(i) = -sign(ydin)*(ydin*ylb >= 0);
  touches.placursign(i) = sign(ydin)*(ydin*ylb < 0);
end

n = touches.count;
for i = 1:nl-1
  if kc(i) & concavel(i+1)
    touches.count = touches.count+1;
    touches.limiter(touches.count) = i+1;
  end
end

for i = n+1:touches.count
  j = touches.limiter(i);
  touches.r(i) = rl(j);
  touches.z(i) = zl(j);
  touches.psi(i) = yl1(j);
  touches.w(i,:) = wl(j,:);
  touches.ii(i,:) = il(j,:);
  touches.drdpsi(i,:) = 0;
  touches.dzdpsi(i,:) = 0;
  touches.limiter(i) = j;
  touches.dpsidinward(i) = 1;
  touches.bdefpossible(i) = 1;
  touches.inicursign(i) = 1;
  touches.placursign(i) = 1;
end

