function [gaps,t, gapspec] = getgaps(shot,gapspec,t1,t2,efit_source,tokamak,idoplot)
%
%  USAGE: [gaps,t, gapspec] = getgaps(shot,gapspec,t1,t2,efit_source,tokamak,idoplot)
%
%  PURPOSE: Get gaps defined by gapspec
%           A gap is defined as the distance between the plasma boundary 
%           and the wall along a line through a specified "gap location", 
%           with direction of measurement defined by gradient of flux for a 
%           nominal equilibrium at that prescribed "gap location".
%
%  INPUTS:    shot: shot number
%          gapspec: gap specification on the form [r z gr gz]
%                   r,z is gap location and gr,gz nominal gradient
%                   default gaps are: inner, upper, outer, lower
%                   see also calc_gaps
%           t1, t2: Specification of the time samples (t)
%                   t = all efit times, if t1,t2 are empty or not supplied
%                   t1 <= t <= t2, if length(t1)==1 & length(t2)==1
%                   t = t1,        if length(t1)>=1 & length(t2)==0
%      efit_source: string with name of efit-tree OR equil_data
%                   default efit01 for shot<900000, efitrt1 for shot>=900000
%          tokamak: default d3d
%          idoplot: flag to plot limiter, gap locations, and gap vectors
%
%  OUTPUTS:   gaps: distances from boundary to wall through r,z along gr,gz [m]
%                t: time [sec]
%          gapspec: gapspec that was used (supplied in call or default)

%
%  RESTRICTIONS: gapspec points must be within the grid

%
%  METHOD: 
%	
%  VERSION @(#)getgaps.m	1.2 03/19/12
%
%  WRITTEN BY:  Anders Welander  ON	5/18/11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  xlim = 0; ylim = 0; % Killing a matlab bug
  
  if nargin<1
    help getgaps
    disp('****************************')
    disp('Supply at least shot number!')
    disp('****************************')
    return
  end
  
  % Defaults
  if ~exist('gapspec','var')
    gapspec = []; % Empty gapspec means default to inner, upper, outer, lower
  end
  if ~isempty(gapspec) & size(gapspec,2) < 4
    wait('error in getgaps: gapspec has to contain 4 columns with gap locations and vectors')
    return
  end 
  if ~exist('t1','var')
    t1 = [];
  end
  if ~exist('t2','var')
    t2 = [];
  end
  if ~exist('efit_source','var')
    if shot<900000
      efit_source = 'efit01';
    else
      efit_source = 'efitrt1';
    end
  end
  if ~exist('tokamak','var')
    tokamak = 'd3d';
  end
  if ~exist('idoplot','var')
    idoplot = 0;
  end
  
  if isstruct(efit_source)
    equilibria = efit_source;
  else
    equilibria = read_eq(shot,[],efit_source,tokamak); % Read the efits
  end
  
  % Select time samples
  if     length(t1) == 0 & length(t2) == 0
    iselect = 1:length(equilibria.time);
  elseif length(t1) == 1 & length(t2) == 0
    [dum, j] = min(abs(equilibria.time-t1)); iselect = j;
  elseif length(t1) == 1 & length(t2) == 1
    iselect = find(t1<=equilibria.time & equilibria.time<=t2);
  elseif length(t1) > 1
    iselect = find(ismember(equilibria.time,t1));
  end
  t = equilibria.time(iselect);
  
  % Look at first equilibrium for limiter information
  if length(equilibria.time) == 1
    eq = equilibria;
  else
    eq = equilibria.gdata(1);
  end
  struct_to_ws(eq);

  % Find weights in grid points to calculate value at a point using cubic Hermite spline (ref. wikipedia)
  mx = [0 2 0 0;-1 0 1 0;2 -5 4 -1;-1 3 -3 1]/2;
  neighbors = reshape([-1-nh -1 -1+nh -1+2*nh;-nh 0 nh 2*nh;1-nh 1 1+nh 1+2*nh;2-nh 2 2+nh 2+2*nh],1,16);
  rgg = ones(nh,1)*rg';
  zgg = zg*ones(1,nw);
  
  % Make default_gapspec with inner, upper, outer, lower gaps
  % Define default inner gap
  k = find(ylim(1:end-1).*ylim(2:end) <= 0);
  [x2, j] = max(xlim(k));
  [x1, j] = min(xlim(k));
  xc = (x1+x2)/2; % Center of limiter
  r = interp1(ylim(k(j)+[0 1]),xlim(k(j)+[0 1]),0);
  default_gapspec(1,:) = [r 0 1 0];
  % Define default upper gap
  k = find((xlim(1:end-1)-xc).*(xlim(2:end)-xc) <= 0 & ylim(1:end-1)>0);
  [dum, j] = min(ylim(k));
  z = interp1(xlim(k(j)+[0 1]),ylim(k(j)+[0 1]),xc);
  default_gapspec(2,:) = [xc z 0 -1];
  % Define default outer gap
  k = find(ylim(1:end-1).*ylim(2:end) <= 0);
  [dum, j] = max(xlim(k));
  r = interp1(ylim(k(j)+[0 1]),xlim(k(j)+[0 1]),0);
  default_gapspec(3,:) = [r 0 -1 0];
  % Define default lower gap
  k = find((xlim(1:end-1)-xc).*(xlim(2:end)-xc) <= 0 & ylim(1:end-1)<0);
  [dum, j] = max(ylim(k));
  z = interp1(xlim(k(j)+[0 1]),ylim(k(j)+[0 1]),xc);
  default_gapspec(4,:) = [xc z 0 1];
  
  if isempty(gapspec)
    gapspec = default_gapspec;
  end
  
  % Make sure the outputs are defined even if length(t) or ngap = 0
  ngap = size(gapspec,1);
  gaps = zeros(ngap,length(t));
  
  % Determine what limiter points the gap vectors intersect and also normalize gr,gz
  nlim = length(xlim);
  for j = 1:ngap
    R = gapspec(j,1);   Z = gapspec(j,2);
    gr = gapspec(j,3); gz = gapspec(j,4);
    g = norm([gr gz]); gr = gr/g; gz = gz/g;
    gapspec(j,3:4) = [gr gz]; % Now they are normalized in the specification
    d2min = inf; % Will be minimum distance^2 along gr,gz from r,z to limiter
    for k = 1:nlim-1
      m = [gr xlim(k+1)-xlim(k); gz ylim(k+1)-ylim(k)];
      if rank(m)==2
	kk = m\[R-xlim(k);Z-ylim(k)];
	if kk(2)>=0 & kk(2)<=1
	  Rt = xlim(k)+kk(2)*(xlim(k+1)-xlim(k));
	  Zt = ylim(k)+kk(2)*(ylim(k+1)-ylim(k));
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
  % Now rgl, zgl holds gap-limiter points
  
  if idoplot
    plot(gapspec(:,1),gapspec(:,2),'go','markerSize',8)
    hold on
    plot(xlim,ylim,'k','linew',2)
    plot(rgl,zgl,'rx','markerSize',12,'linew',2)
    for j = 1:ngap
      gr = gapspec(j,3); gz = gapspec(j,4);
      rs(j) = rgl(j)+5*dr*gr; zs(j) = zgl(j)+5*dr*gz;
      plot([rgl(j) rs(j)],[zgl(j) zs(j)],'m')
      a = angle(gr+i*gz);
      plot(rs(j)-[0 2*dr]*cos(a-pi/12),zs(j)-[0 2*dr]*sin(a-pi/12),'m')
      plot(rs(j)-[0 2*dr]*cos(a+pi/12),zs(j)-[0 2*dr]*sin(a+pi/12),'m')
      a = angle((rs(j)-rmaxis)+i*(zs(j)-zmaxis));
      rho = sqrt((rs(j)-rmaxis)^2+(zs(j)-zmaxis)^2);
      text(rmaxis+0.8*rho*cos(a),zmaxis+0.8*rho*sin(a),num2str(j),'fonts',12,'fontw','bold','color','g','hori','center')
    end
    axis('image')
    legend('Gap locations');
    hold off
  end
  
  % Begin processing the equilibria
  ibadbbbs = [];
  for jeq = 1:length(iselect)
    % Get this eq
    if length(equilibria.time) == 1
      eq = equilibria;
    else
      eq = equilibria.gdata(iselect(jeq));
    end
    struct_to_ws(eq);
    
    if nbbbs > 0
      % Find where gap vectors intersect bbbs
      for j = 1:ngap
	R = gapspec(j,1);   Z = gapspec(j,2);
	gr = gapspec(j,3); gz = gapspec(j,4);
	d2min = inf; % Will be minimum distance^2 along gr,gz from r,z to limiter
	rgb(j) = R; zgb(j) = Z; % dummy values
	for k = 1:nbbbs-1
	  m = [gr rbbbs(k+1)-rbbbs(k); gz zbbbs(k+1)-zbbbs(k)];
	  if rank(m)==2
	    kk = m\[R-rbbbs(k);Z-zbbbs(k)];
	    if kk(2)>=0 & kk(2)<1
	      Rt = rbbbs(k)+kk(2)*(rbbbs(k+1)-rbbbs(k));
	      Zt = zbbbs(k)+kk(2)*(zbbbs(k+1)-zbbbs(k));
	      d2 = (R-Rt)^2+(Z-Zt)^2;
	      if d2<d2min
		rgb(j) = Rt;
		zgb(j) = Zt;
		d2min = d2;
	      end
	    end
	  end
	end
	if 01
	R = rgb(j); Z = zgb(j); % Boundary point found by intersection
	dstep = 1e-3; dtot = 0; % dstep is initial step size in correcting rgb,zgb. dtot= step length 
	psiprev = inf;
	% Refine the boundary points in rgb, zgb by searching for exact point where psi = psibry
	while abs(dstep) > 1e-6 & dtot<dr
	  R = R+gr*dstep; Z = Z+gz*dstep; dtot = dtot+dstep;
	  k = 1+floor((R-rg(1))/dr)*nh+floor((Z-zg(1))/dz);
	  ii = k+neighbors; % ii indexes 16 points around R,Z
	  tr = (R-rgg(k))/dr; tz = (Z-zgg(k))/dz;
	  w = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
	  psi = psizr(ii)*w';
	  if (psi>psibry & dstep>0) | (psi<psibry & dstep<0) | (abs(psi-psibry)>abs(psiprev-psibry))
	    dstep = -dstep/2;
	  end
	  psiprev = psi;
	end
	rgb(j) = R; zgb(j) = Z;
	end
	gaps(j,jeq) = sqrt((rgb(j)-rgl(j))^2+(zgb(j)-zgl(j))^2);
      end
    else
      ibadbbbs(end+1) = jeq;
    end
  end
  igood = setdiff(1:length(t),ibadbbbs);
  t = t(igood); gaps = gaps(:,igood);
  
