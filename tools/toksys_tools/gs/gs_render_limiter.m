%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE:   gs_render_limiter
%
%  PURPOSE: Calculate how pixels are affected by grid elements of plasma
%
%  INPUTS:  xcam, ycam, zcam, nxpix, nzpix, rl, zl, Plg
%
%  OUTPUTS: Pcg, power on ccd from plasma on grid
%
%  METHOD:   
	
%  VERSION %W% %G%
%
%  WRITTEN BY:  Anders Welander  ON	6/30/14
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check Plg (power of radiation from limiter elements caused by grid element)
if ~exist('Plg','var')
  disp('warning gs_render_limiter: Plg is missing')
  disp('Plg = power of radiation from limiter elements caused by grid elements')
  disp('Calculating Plg with gs_limiter_radiation...')
  gs_limiter_radiation
  disp('Done')
end
if size(Plg,1) ~= nl || size(Plg,2) ~= ngg
  disp('warning gs_render_limiter: Plg is the wrong size')
  disp('Plg = power of radiation from limiter elements caused by grid elements')
  disp('Calculating new Plg with gs_limiter_radiation...')
  gs_limiter_radiation
  disp('Done')
end

% Limiter center & size
Rl = (max(rl)+min(rl))/2;
Zl = (max(zl)+min(zl))/2;
al = (max(rl)-min(rl))/2;
bl = (max(zl)-min(zl))/2;

% Camera point
if ~exist('xcam','var') || isempty(xcam)
  xcam = Rl;
end
if ~exist('ycam','var') || isempty(ycam)
  ycam = -al;
end
if ~exist('zcam','var') || isempty(zcam)
  zcam = Zl;
end

% Pixels in plane y = 0
if ~exist('nxpix','var') | isempty(nxpix)
  nxpix = 64;
end
if ~exist('nzpix','var') | isempty(nzpix)
  nzpix = 128;
end
xpix = Rl + linspace(-al,al,nxpix);
zpix = Zl + linspace(-bl,bl,nzpix);

% For fast rendering
Pcg = zeros(nzpix*nxpix,ngg);

% Fraction of each pixel area within limiter
Apix = polycellint(xpix, zpix, rl, zl);
fpix = Apix/max(Apix(:));

% Trace
drl = diff(rl);
dzl = diff(zl);
rl1 = rl(1:nl-1);
c00 = rl1.*rl1;
c01 = rl1.*drl;
c11 = drl.*drl;
dydt = -ycam;
dzlc = zl(1:nl-1) - zcam;
for i = 1:nzpix
  dzdt = zpix(i)-zcam;
  for j = 1:nxpix
    dxdt = xpix(j)-xcam;
    % r = rl + drl*f
    % z = zl + dzl*f
    % x = xcam + dxdt*t
    % y = ycam + dydt*t
    % z = zcam + dzdt*t
    % zcam + dzdt*t = zl + dzl*f
    % t = (zl + dzl*f - zcam)/dzdt = (dzlc + dzl*f)/dzdt
    % x = xcam + dxdt*(dzlc + dzl*f)/dzdt = a0 + a1*f
    a0 = xcam + dxdt/dzdt*dzlc;
    a1 = dxdt/dzdt*dzl;
    % y = ycam + dydt*(dzlc + dzl*f)/dzdt = b0 + b1*f
    b0 = ycam + dydt/dzdt*dzlc;
    b1 = dydt/dzdt*dzl;
    % (a0 + a1*f)^2 + (b0 + b1*f)^2 = (rl + drl*f)^2
    al = 2*(a1.*a1 + b1.*b1 - c11);
    bl = 2*(a0.*a1 + b0.*b1 - c01);
    cl = 2*(a0.*a0 + b0.*b0 - c00);
    sq = sqrt(bl.*bl - al.*cl);
    f1 = -(bl-sq)./al;
    f2 = -(bl+sq)./al;
    t1 = (dzlc + dzl.*f1)/dzdt;
    t2 = (dzlc + dzl.*f2)/dzdt;
    t = inf;
    for k = 1:nl-1
      if isreal(f1(k)) && f1(k) > 0 && f1(k) < 1 && t1(k) > 0 && t1(k) < t
        t = t1(k);
	f = f1(k);
	n = k;
      end
      if isreal(f2(k)) && f2(k) > 0 && f2(k) < 1 && t2(k) > 0 && t2(k) < t
        t = t2(k);
	f = f2(k);
	n = k;
      end
    end
    xt = xcam + dxdt*t;
    yt = ycam + dydt*t;
    zt = zcam + dzdt*t;
    rt = sqrt(xt^2+yt^2);
    t3 = [dxdt dydt dzdt];
    l3 = [dzl(n)*xt/rt dzl(n)*yt/rt -drl(n)];
    proj = abs(t3*l3'/norm(t3)/norm(l3));
    k = i+(j-1)*nzpix;
    Pcg(k,:) = proj*(f*Plg(n+1,:) + (1-f)*Plg(n,:));
  end
end



return

b = zeros(nzpix,nxpix,3);
im = Pcg*pcurrt(:);
a = reshape(im,nzpix,nxpix)/max(im);
ar = a*0.9;
ag = a*0.8;
ab = a*1.0;
b(:,:,1) = ar.*fpix + 1-fpix;
b(:,:,2) = ag.*fpix + 1-fpix;
b(:,:,3) = ab.*fpix + 1-fpix;

clf
hold on
image(xpix,zpix,b)
plot(rl,zl,'g','linew',4)
axis image
h = get(gcf,'Children');
PlotBoxAspectRatio = get(get(gcf,'Children'),'PlotBoxAspectRatio')
axis_read = axis
axis_write = ...
  [mean(axis_read(1:2))+PlotBoxAspectRatio(2)^2*[axis_read(1:2)-mean(axis_read(1:2))], ...
   mean(axis_read(3:4))+PlotBoxAspectRatio(1)^2*[axis_read(3:4)-mean(axis_read(3:4))]]
axis(axis_write)
PlotBoxAspectRatio = get(h,'PlotBoxAspectRatio')
  set(gcf,'ResizeFcn', [...
    'drawnow' 10 ...
    'axis image' 10 ...
    'drawnow' 10 ...
    'PlotBoxAspectRatio = get(get(1,''Children''),''PlotBoxAspectRatio'')' 10 ...
    'axis_read = axis;' 10 ...
    'axis_write = ' ...   
    '[mean(axis_read(1:2))+PlotBoxAspectRatio(2)^2*' ...
    '[axis_read(1:2)-mean(axis_read(1:2))], ' ...
    'mean(axis_read(3:4))+PlotBoxAspectRatio(1)^2*' ...
    '[axis_read(3:4)-mean(axis_read(3:4))]]' 10 ...
    'axis(h,axis_write)'])

return

tic
% This HD is pretty useless
hd.scalex = 1;
hd.scalez = 1;
hd.xpix = linspace(xpix(1),xpix(end),nxpix*hd.scalex);
hd.zpix = linspace(zpix(1),zpix(end),nzpix*hd.scalez);
hd.ar = gs_interp2(xpix, zpix, ar, hd.xpix, hd.zpix');
hd.ag = gs_interp2(xpix, zpix, ag, hd.xpix, hd.zpix');
hd.ab = gs_interp2(xpix, zpix, ab, hd.xpix, hd.zpix');
hd.ar(hd.ar<0) = 0;
hd.ag(hd.ag<0) = 0;
hd.ab(hd.ab<0) = 0;
hd.ar(hd.ar>1) = 1;
hd.ag(hd.ag>1) = 1;
hd.ab(hd.ab>1) = 1;
hd.b        = hd.ar;
hd.b(:,:,2) = hd.ag;
hd.b(:,:,3) = hd.ab;
clf
hold on
image(hd.xpix,hd.zpix,hd.b)
plot(rl,zl,'g','linew',4)
%axis image
toc


