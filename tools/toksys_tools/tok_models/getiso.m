function [iso,t] = getiso(shot,rziso,t1,t2,efit_source,tokamak,idoplot)
%
%  USAGE: [iso,t] = getiso(shot,rziso,t1,t2,efit_source,tokamak,idoplot)
%
%  PURPOSE: Get iso fluxes and fields at rziso points
%
%  INPUTS:    shot: shot number
%            rziso: R, Z of isoflux points arranged as [R(:) Z(:)] (unit: meters)
%           t1, t2: Specification of the time samples (t)
%                   t = all efit times, if t1,t2 are empty or not supplied
%                   t1 <= t <= t2, if length(t1)==1 & length(t2)==1
%                   t = t1,        if length(t1)>=1 & length(t2)==0
%      efit_source: string with name of efit-tree OR equil_data
%                   default efit01 for shot<900000, efitrt1 for shot>=900000
%          tokamak: default d3d
%          idoplot: flag to plot limiter and rziso points
%
%  OUTPUTS:    iso: structure with fields:
%                   psi: flux at rziso points
%                   br, bz: Br, Bz at rziso points
%                t: time [sec]

%
%  RESTRICTIONS: rziso points must be within the grid

%
%  METHOD: 
%	
%  VERSION @(#)getiso.m	1.3 09/20/11
%
%  WRITTEN BY:  Anders Welander  ON	5/20/11
%
%  MODIFICATION HISTORY:				
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  xlim = 0; ylim = 0; % Killing a matlab bug
  
  if nargin<2
    help getiso
    disp('*********************************************')
    disp('Supply at least shot number and rziso points!')
    disp('*********************************************')
    return
  end
  
  % Defaults
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
  
  % Make sure the outputs are defined even if length(t) or niso = 0
  niso = size(rziso,1);
  iso = [];
  
  if idoplot
    plot(rziso(:,1),rziso(:,2),'rx','markerSize',8)
    hold on
    plot(xlim,ylim,'k','linew',2)
    for j = 1:niso
      a = angle((rziso(j,1)-rmaxis)+i*(rziso(j,2)-zmaxis));
      rho = sqrt((rziso(j,1)-rmaxis)^2+(rziso(j,2)-zmaxis)^2);
      text(rmaxis+0.8*rho*cos(a),zmaxis+0.8*rho*sin(a),num2str(j),'fonts',12,'fontw','bold','color','g','hori','center')
    end
    axis('image')
    legend('isoflux points');
    hold off
  end
  
  % Begin processing the equilibria
  for jeq = 1:length(iselect)
    % Get this eq
    if length(equilibria.time) == 1
      eq = equilibria;
    else
      eq = equilibria.gdata(iselect(jeq));
    end
    struct_to_ws(eq);
    for j = 1:niso
      kr0 = min(nw-3,max(1,floor((rziso(j,1)-rg(1))/dr))); % r index 0-start, allowed values: 1:nw-3
      kz1 = min(nh-2,max(2,ceil((rziso(j,2)-zg(1))/dz))); % z index 1-start, allowed values: 2:nh-2
      k = kr0*nh+kz1;
      ii = k+neighbors; % ii indexes 16 points around rziso(j,:)
      tr = (rziso(j,1)-rgg(k))/dr; tz = (rziso(j,2)-zgg(k))/dz;
      w = reshape(([1 tz tz^2 tz^3]*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      wr = reshape(([1 tz tz^2 tz^3]*mx)'*[0 1 2*tr 3*tr^2]/dr*mx,1,16);
      wz = reshape(([0 1 2*tz 3*tz^2]/dz*mx)'*[1 tr tr^2 tr^3]*mx,1,16);
      iso.psi(j,jeq) = psizr(ii)*w';
      iso.br(j,jeq) = -psizr(ii)*wz'/rziso(j,1);
      iso.bz(j,jeq) = +psizr(ii)*wr'/rziso(j,1);
    end
  end
  
