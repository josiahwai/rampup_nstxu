function [rx,zx,psix,t,br,bz,eqs] = getxpt(shot,t1,t2,efit_source,tokamak,xtarget,idoplot)
%
%  USAGE: [rx,zx,psix,t,br,bz,eqs] = getxpt(shot,t1,t2,efit_source,tokamak,xtarget,idoplot)
%
%  PURPOSE: Get coordinates for 1 x-point below and 1 above the magnetic axis.
%           The points where flux is closest to the boundary flux are chosen.
%           Also return the flux at these x-points and at the boundary.
%           In addition return Br, Bz at xtarget points (if specified)
%
%  INPUTS:    shot: shot number
%           t1, t2: Specification of the time samples (t)
%                   t = all efit times, if t1,t2 are empty or not supplied
%                   t1 <= t <= t2, if length(t1)==1 & length(t2)==1
%                   t = t1,        if length(t1)>=1 & length(t2)==0
%      efit_source: Structure with equilibria or name of MDS tree or directory
%                   with gfiles or corsica flat files (see help read_eq),
%                   default EFIT01 for shot<900000, EFITRT1 for shot>=900000
%          tokamak: default D3D
%          xtarget: specification on the form [r z] where r and z are column vectors
%                   of target x-point locations, for calculation of Br, Bz.
%          idoplot: flag to plot x-point trajectories and limiter,
%                   bit 0 plots x's, bit 1 plots limiter, bit 2 adds text,
%                   bit 3 adds color coding of times, bit 4 adds colorbar,
%                   idoplot=31 plots it all, default is 0.
%
%  OUTPUTS:     rx: radial position of x-points [lower; upper] [m]
%               zx: vertical position of x-points [lower; upper] [m]
%             psix: fluxes at [lower x; upper x; boundary] [Wb]
%                t: time [sec]
%               br: radial magnetic field at xtarget points
%               bz: vertical magnetic field at xtarget points
%              eqs: The structure with equilibria from efit_source

%
%  RESTRICTIONS: 

%
%  METHOD: 
%	
%  VERSION @(#)getxpt.m	1.3 02/01/12
%
%  WRITTEN BY:  Anders Welander  ON	1/31/12
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  xlim = 0; ylim = 0; % Killing a matlab bug
  
  if nargin<1
    help getxpt
    disp('****************************')
    disp('Supply at least shot number!')
    disp('****************************')
    return
  end
  
  % Defaults
  if ~exist('t1','var')
    t1 = [];
  end
  if ~exist('t2','var')
    t2 = [];
  end
  if ~exist('efit_source','var') | isempty(efit_source)
    if shot<900000
      efit_source = 'EFIT01';
    else
      efit_source = 'EFITRT1';
    end
  end
  if ~exist('tokamak','var') | isempty(tokamak)
    tokamak = 'D3D';
  end
  if ~exist('idoplot','var')
    idoplot = 0;
  end
  
  % Read the efits
  if isstruct(efit_source)
    eqs = efit_source;
  else
    eqs = read_eq(shot,[],efit_source,tokamak);
  end
  
  % Select time samples
  if     length(t1) == 0 & length(t2) == 0
    iselect = 1:length(eqs.time);
  elseif length(t1) == 1 & length(t2) == 0
    [dum, j] = min(abs(eqs.time-t1)); iselect = j;
  elseif length(t1) == 1 & length(t2) == 1
    iselect = find(t1<=eqs.time & eqs.time<=t2);
  elseif length(t1) > 1
    iselect = find(ismember(eqs.time,t1));
  end
  t = eqs.time(iselect);
  
  % Look at first equilibrium for limiter information
  if length(eqs.time) == 1
    eq = eqs;
  else
    eq = eqs.gdata(1);
  end
  struct_to_ws(eq);

  % Find weights in grid points to calculate value at a point using cubic Hermite spline (ref. wikipedia)
  mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;
  neighbors = reshape([-1-nh -1 -1+nh -1+2*nh;-nh 0 nh 2*nh;1-nh 1 1+nh 1+2*nh;2-nh 2 2+nh 2+2*nh],1,16);
  rgg = ones(nh,1)*rg';
  zgg = zg*ones(1,nw);
    
  % Make sure the outputs are defined
  rx = zeros(2,length(iselect));
  zx = zeros(2,length(iselect));
  psix = zeros(3,length(iselect));
  br = zeros(0,length(iselect));
  bz = zeros(0,length(iselect));
  
  % Begin processing the equilibria
  ibadbbbs = [];
  for jeq = 1:length(iselect)
    % Get this eq
    if length(eqs.time) == 1
      eq = eqs;
    else
      eq = eqs.gdata(iselect(jeq));
    end
    struct_to_ws(eq);
    
    if nbbbs > 0
      % Find rx, zx
      yr = [zeros(nh,1) psizr(:,3:end)-psizr(:,1:end-2) zeros(nh,1)]/2/dr;
      yz = [zeros(1,nw);psizr(3:end,:)-psizr(1:end-2,:);zeros(1,nw)]/2/dz;
      iNULL = [];
      for j = 3:nw-2
        for k = 3:nh-2
	  if (yr(k,j-1)<0 | yr(k,j+1)<0 | yr(k-1,j)<0 | yr(k+1,j)<0) & ...
	     (yr(k,j-1)>0 | yr(k,j+1)>0 | yr(k-1,j)>0 | yr(k+1,j)>0) & ...
	     (yz(k,j-1)<0 | yz(k,j+1)<0 | yz(k-1,j)<0 | yz(k+1,j)<0) & ...
	     (yz(k,j-1)>0 | yz(k,j+1)>0 | yz(k-1,j)>0 | yz(k+1,j)>0)
	    iNULL(end+1) = k+(j-1)*nh;
	  end
        end
      end
      psix(3,jeq) = psibry;
      for ix = 1:2
	if ix == 1
	  ii = find(zgg(iNULL) < zmaxis);  % Search for lower x-point
	  [psid,k] = min(abs(psizr(iNULL(ii))-psibry));
	else
	  ii = find(zgg(iNULL) > zmaxis);  % Search for upper x-point
	  [psid,k] = min(abs(psizr(iNULL(ii))-psibry));
	end
	k = iNULL(ii(k));
	r = rgg(k); z = zgg(k);
	cx=0;
	c_previous = zg(end)-zg(1); cx = c_previous*.99; j = 9;
	% Try zooming in on x-point with Newton Rhapson.
	while j>0 & norm(cx)<norm(c_previous)
	  % Find indices and weights for grid points around the x-point
	  kr0 = min(nw-3,max(1,floor((r-rg(1))/dr))); % r index 0-start, allowed values: 1:nw-3
	  kz1 = min(nh-2,max(2,ceil((z-zg(1))/dz))); % z index 1-start, allowed values: 2:nh-2
	  k = kr0*nh+kz1;
	  iix = k+neighbors; % iix indexes 16 points around x point
	  pp = psizr(iix); j = j-1; c_previous = cx;
	  tr = (r-rgg(k))/dr; tz = (z-zgg(k))/dz;
	  wx = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	  wxr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	  wxz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	  wxrr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 0 2 6*tr]/dr^2*mx,1,16);
	  wxzz = reshape(([0 0 2 6*tz]/dz^2*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	  wxrz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	  xshift_Bx = -inv([pp*wxrr' pp*wxrz';pp*wxrz' pp*wxzz']);
	  cx = xshift_Bx*[pp*wxr'; pp*wxz'];
	  r = r+cx(1); z = z+cx(2);
	end
	% Test that we converged on a NULL and that it wasn't the axis
	if norm(cx)<dr & (pp*wx'-psimag)/(psibry-psimag)>0.5
	  rx(ix,jeq) = r;  zx(ix,jeq) = z;
	  psix(ix,jeq) = pp*wx';
	else
	  rx(ix,jeq) = NaN; zx(ix,jeq) = NaN;
	  psix(ix,jeq) = NaN;
	end
      end
      % Find br, bz at xtarget points
      for ix = 1:size(xtarget,1)
        r = xtarget(ix,1); z = xtarget(ix,2);
	kr0 = min(nw-3,max(1,floor((r-rg(1))/dr))); % r index 0-start, allowed values: 1:nw-3
	kz1 = min(nh-2,max(2,ceil((z-zg(1))/dz))); % z index 1-start, allowed values: 2:nh-2
	k = kr0*nh+kz1;
	iix = k+neighbors; % iix indexes 16 points around x point
	pp = psizr(iix); j = j-1; c_previous = cx;
	tr = (r-rgg(k))/dr; tz = (z-zgg(k))/dz;
	wxr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
	wxz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	br(ix,jeq) = -pp*wxz'/r/2/pi;
	bz(ix,jeq) = +pp*wxr'/r/2/pi;
      end
    else
      ibadbbbs(end+1) = jeq;
    end
  end
  igood = setdiff(1:length(t),ibadbbbs);
  t = t(igood); rx = rx(:,igood); zx = zx(:,igood);
  br = br(:,igood); bz = bz(:,igood);
      
  if idoplot
    clf
    hold on
    plot(rx',zx','x-')
    title(['x-point trajectories, ' num2str(t(1)) ' to ' num2str(t(end)) ' sec'])
    if bitand(idoplot,2)
      plot(xlim,ylim,'k','linew',2)
      axis('image')
    end
    if bitand(idoplot,4) % display time next to first and last x-point in trajectories
      k = find(~isnan(rx(1,:)));
      if length(k) > 2
	text(rx(1,k(1)),zx(1,k(1)),['t = ' num2str(t(k(1)))])
	text(rx(1,k(end)),zx(1,k(end)),['t = ' num2str(t(k(end)))])
      end
      k = find(~isnan(rx(2,:)));
      if length(k) > 2
	text(rx(2,k(1)),zx(2,k(1)),['t = ' num2str(t(k(1)))])
	text(rx(2,k(end)),zx(2,k(end)),['t = ' num2str(t(k(end)))])
      end
    end
    if bitand(idoplot,8) % plot the x's with different colors
      plot(rx',zx','k-')
      rgb = interp1(linspace(t(1),t(end),64),jet,t);
      for j = 1:length(t)
	plot(rx(:,j),zx(:,j),'x','color',rgb(j,:))
      end
    end
    if bitand(idoplot,16) % Add a colorbar that shows colors of times
      j = 9;
      it1 = ceil(t(1)/10^j);
      it2 = floor(t(end)/10^j);
      while it2<it1
	j = j-1;
        it1 = ceil(t(1)/10^j);
        it2 = floor(t(end)/10^j);
      end
      tc = (it1:it2)*10^j;
      h = colorbar;
      v = interp1(t,linspace(1,64,length(t)),tc);
      set(h,'YTick',v);
      a = '';
      for j = 1:length(v)
        s = num2str(tc(j));
	a(j,1:9+length(s)) = [' t = ' s ' sec'];
      end
      set(h,'YTickLabel',a)
    end
    hold off
  end
  

